import logger from "./logger";
import { Molecule, BaseMolecule, User, StashedMolecule } from "./Entities/entities";
import { simpleflake } from 'simpleflakes';
import Errors, { ApiError, ErrorType } from "./Errors";
import Express from 'express';
import { TokenPayload, ReadedFile } from "./types";
import JsonWebToken from 'jsonwebtoken';
import { KEYS, UPLOAD_ROOT_DIR } from "./constants";
import { Database } from "./Entities/CouchHelper";
import { unlink } from "fs";
import MoleculeOrganizer from "./MoleculeOrganizer";
import SearchWorker from "./search_worker";
import Mailer from "./Mailer/Mailer";
import fs, { promises as FsPromise } from 'fs';
import path from 'path';

export function isDebugMode() {
  return logger.level === "debug" || logger.level === "silly";
}

export function isMolecule(e: BaseMolecule) : e is Molecule {
  return 'approved_by' in e;
}

export function generateSnowflake() {
  return simpleflake(undefined, undefined, Date.UTC(2020, 0, 1)).toString(10);
}

export function sendError(error: ApiError, res: Express.Response) {
  res.status(Number(error.message)).json(error.data);
}

export function cleanMulterFiles(req: Express.Request) {
  if (req.files) {
    if (Array.isArray(req.files)) {
      for (const file of req.files) {
        unlink(file.path, () => {});
      }
    }
    else {
      for (const files of Object.values(req.files)) {
        for (const file of files) {
          unlink(file.path, () => {});
        }
      }
    }
  }
  if (req.file) {
    unlink(req.file.path, () => {});
  }
}

/**
 * Create a closure to catch an ApiError that occurs in a Promise.
 */
export function errorCatcher(res: Express.Response, req?: Express.Request) {
  return function(err: any) {

    if (req) {
      cleanMulterFiles(req);
    }
  
    if (res.headersSent) {
      return;
    }

    if (err instanceof ApiError) {
      return sendError(err, res);
    }

    if (err instanceof Error) {
      logger.error("During request handling, the following error occurred: " + err + "\n" + err.stack);
    }
    else if (typeof err === 'string') {
      logger.error("Unknown error: " + err);
    }
    else if (Array.isArray(err)) {
      logger.error("Unknown error: " + err.join(', '));
    }
    else if (err) {
      logger.error("Unknown error: " + JSON.stringify(Object.getOwnPropertyDescriptors(err), null, 2));
    }
    else {
      logger.error("Unknown error (undefined)");
    }
    
    return sendError(Errors.make(ErrorType.Server), res);
  };
}

export function sanitize(obj: any) {
  const props = [] as string[];

  for (const prop in obj) {
    if (prop.startsWith('_')) {
      props.push(prop);
    }
  }

  for (const prop of props) {
    delete obj[prop];
  }

  return obj;
}

export function signToken(payload: TokenPayload, id: string) {
  return new Promise((resolve, reject) => {
    // Signe le token
    JsonWebToken.sign(
      payload, // Données custom
      { key: KEYS.PRIVATE, passphrase: "" }, // Clé RSA privée
      { 
        algorithm: 'RS256', 
        expiresIn: "720d", // 2 years durability
        issuer: "MArtini Database Server 1", 
        jwtid: id, // ID généré avec snowflake
      }, 
      (err, encoded) => { // Quand le token est généré (ou non), accepte/rejette la promesse
        if (err) reject(err);
        else resolve(encoded);
      }
    );
  }) as Promise<string>;
}

export async function validateToken(token: string) {
  const payload: any = await new Promise((resolve, reject) => {
    JsonWebToken.verify(
      token, 
      { key: KEYS.PUBLIC, passphrase: "" }, 
      { algorithms: ['RS256'] }, 
      (err, payload: any) => {
        if (err) {
          reject(err);
          return;
        }

        resolve(payload);
      }
    )
  });

  return getUserFromToken(payload.jti);
}

export function getUserFromToken(jti: string) {
  // Get the token from string and call done(null, is_revoked)
  return Database.token.get(jti as string)
    .then(() => Database.user.fromToken(jti as string));
}

export function methodNotAllowed(allow: string | string[]) {
  return (_: any, res: Express.Response) => {
      res.setHeader('Allow', typeof allow === 'string' ? allow : allow.join(', '));
      Errors.throw(ErrorType.InvalidMethod);
  };
}

export function getNameAndPathOfUploadedFile(name: string) : [string, string] {
  const file_name = name.includes('/') ? name.split('/').pop()! : name;
  const file_path = UPLOAD_ROOT_DIR + file_name;

  return [file_name, file_path];
}

export function escapeRegExp(string: string) {
  return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); // $& means the whole matched string
}

export function sleep(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

export function withRegex(text: string, is_regex: boolean, flags = "i", strict_match = false) {
  let search_text = is_regex ? text : escapeRegExp(text);

  if (strict_match) {
    search_text = "^" + search_text + "$";
  }

  return { $regex: `(?${flags})${search_text}` };
}

export async function deleteMolecule(id: string, user: User, stashed = false, checked_attached = true) {
  // Delete a stashed molecule
  
  const getChilds = async (id: string, stashed = false) : Promise<BaseMolecule[]> => {
    const db = stashed ? Database.stashed : Database.molecule
    const molDoc = await db.get(id)
    const tree = await db.moleculeTreeOf(molDoc.tree_id)
    if(tree) {
      const childs = tree.getChilds(molDoc)
      return childs
    }
    
    return []
  }

  const dels = []
  if (stashed) {
    if (user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    try {
      const mol = await Database.stashed.get(id)
      const childMols = await getChilds(mol.id, true)

      // Delete attached ZIP
      await MoleculeOrganizer.remove(mol.files);
      const delResp = await Database.stashed.delete(mol);
      dels.push(delResp.id)
      for (const child of childMols){
        const molDoc = await Database.stashed.get(child.id)
        await MoleculeOrganizer.remove(child.files); 
        const delResp = await Database.stashed.delete(molDoc);
        dels.push(delResp.id)
      }

    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }
    return;
  }

  // Delete a published molecule
  try {
    const mol = await Database.molecule.get(id);
    const childMols = await getChilds(mol.id)
    if (user.role !== "admin" && user.id !== mol.owner) {
      return Errors.throw(ErrorType.Forbidden);
    }

    await MoleculeOrganizer.remove(mol.files)

    const delResp = await Database.molecule.delete(mol)
    dels.push(delResp.id)

    for (const mol of childMols){
      const molDoc = await Database.molecule.get(mol.id)
      await MoleculeOrganizer.remove(mol.files)
      const delResp = await Database.molecule.delete(molDoc)
      dels.push(delResp.id)
    }

    return dels

    /*if (checked_attached) {
      // Recherche les sous-versions attachées à cette molecule
      const versions_attached = await Database.molecule.find({
        limit: 99999,
        selector: {
          parent: mol.id,
          tree_id: mol.tree_id,
        },
      });
  
      // Met à jour les liens de parenté
      if (versions_attached.length) {
        if (mol.parent) {
          // If the deleted molecule has a parent: Every molecule will be attached to the new parent
          for (const m of versions_attached) {
            m.parent = mol.parent;
          }
        }
        else {
          // If the deleted molecule doesn't have a parent: The first molecule will be the new parent for everyone
          const first = versions_attached[0];
          for (const other of versions_attached.slice(1)) {
            other.parent = first.id;
          }
          first.parent = null;
        }
  
        // Save everyone
        await Promise.all(versions_attached.map(v => Database.molecule.save(v)));
        SearchWorker.clearCache();
      }
    }

    // Delete attached ZIP
    await MoleculeOrganizer.remove(mol.files);
    await Database.molecule.delete(mol); */
  } catch (e) {
    return Errors.throw(ErrorType.ElementNotFound);
  }
}

export async function informAdminFromAskCreation(new_user: User) {
  const admins = await Database.user.find({
    selector: { role: "admin" },
    limit: 99999
  });

  const promises: any[] = [];
  for (const admin of admins) {
    // Send a mail
    promises.push(Mailer.send({
      to: admin.email,
      subject: "MArtini Database - New account request",
    }, "mail_ask", {
      new_user: {
        name: new_user.name
      },
      name: admin.name
    }));
  }

  await Promise.all(promises);
}

export async function informAdminFromNewMolecule(new_molecule: StashedMolecule, submitter: User) {
  const admins = await Database.user.find({
    selector: { role: "admin" },
    limit: 99999
  });

  const promises: any[] = [];
  for (const admin of admins) {
    // Send a mail
    promises.push(Mailer.send({
      to: admin.email,
      subject: "MArtini Database - New molecule submitted",
    }, "mail_molecule_submitted", {
      submitter,
      name: admin.name,
      molecule: new_molecule
    }));
  }

  await Promise.all(promises);
}

export async function informAdminContact(content: string, sender: string) {
  const admins = await Database.user.find({
    selector: { role: "admin" },
    limit: 99999
  });

  const promises: any[] = [];
  for (const admin of admins) {
    // Send a mail
    promises.push(Mailer.send({
      to: admin.email,
      subject: "MArtini Database - New question asked from contact page",
    }, "mail_contact", {
      content,
      name: admin.name,
      email: sender
    }));
  }

  await Promise.all(promises);
}

// Example: 4.700000e-01    4.990000e+00
export function vanDerWaalsRadius(lennar_johns_1: number, lennar_johns_2: number) {
  const T = 0.191, S = 0.230, R = 0.254;


}

export async function dumpStdFromDir(dir: string) {
  // Can for .stdout and .stderr files
  const stderr: any = {};
  const stdout: any = {};
  
  for (const file of await FsPromise.readdir(dir)) {
    if (file.endsWith('.stdout')) {
      const name = file.slice(0, file.length - '.stdout'.length);
      stdout[name] = await FsPromise.readFile(dir + '/' + file, 'utf-8');
    }
    else if (file.endsWith('.stderr')) {
      const name = file.slice(0, file.length - '.stderr'.length);
      stderr[name] = await FsPromise.readFile(dir + '/' + file, 'utf-8');
    }
  }

  return { stdout, stderr };
}

/**
 * Get the basename of a file without the extension.
 */
export function basenameWithoutExt(src: string) {
  const basename = path.basename(src);
  const last_dot = basename.lastIndexOf('.');

  if (last_dot !== -1) {
    return basename.slice(0, last_dot);
  }
  return basename;
}

/**
 * Create a type that contain values from an array (known at compile-time).
 * 
 * Usage:
 * ```ts
 * const MY_VALUES = ['v1', 'v2', 'v3'] as const;
 * 
 * type AvailableValues = ArrayValues<typeof MY_VALUES>;
 * // AvailableValues = 'v1' | 'v2' | 'v3'
 * ```
 */
export type ArrayValues<T extends ReadonlyArray<unknown>> = T extends ReadonlyArray<infer ArrayValues> ? ArrayValues : never;

export function fileExists(path: string) {
  return FsPromise.access(path, fs.constants.F_OK).then(() => true).catch(() => false);
}

/**
 * Formate un objet Date en chaîne de caractères potable.
 * Pour comprendre les significations des lettres du schéma, se référer à : http://php.net/manual/fr/function.date.php
 * @param schema string Schéma de la chaîne. Supporte Y, m, d, g, H, i, s, n, N, v, z, w
 * @param date Date Date depuis laquelle effectuer le formatage
 * @returns string La chaîne formatée
 */
export function dateFormatter(schema: string, date = new Date()) : string {
  function getDayOfTheYear(now: Date): number {
    const start = new Date(now.getFullYear(), 0, 0);
    const diff = now.getTime() - start.getTime();
    const oneDay = 1000 * 60 * 60 * 24;
    const day = Math.floor(diff / oneDay);

    return day - 1; // Retourne de 0 à 364/365
  }

  const Y = date.getFullYear();
  const N = date.getDay() === 0 ? 7 : date.getDay();
  const n = date.getMonth() + 1;
  const m = (n < 10 ? "0" : "") + String(n);
  const d = ((date.getDate()) < 10 ? "0" : "") + String(date.getDate());
  const L = Y % 4 === 0 ? 1 : 0;

  const i = ((date.getMinutes()) < 10 ? "0" : "") + String(date.getMinutes());
  const H = ((date.getHours()) < 10 ? "0" : "") + String(date.getHours());
  const g = date.getHours();
  const s = ((date.getSeconds()) < 10 ? "0" : "") + String(date.getSeconds());

  const replacements: any = {
    Y, m, d, i, H, g, s, n, N, L, v: date.getMilliseconds(), z: getDayOfTheYear, w: date.getDay()
  };

  let str = "";

  // Construit la chaîne de caractères
  for (const char of schema) {
    if (char in replacements) {
      if (typeof replacements[char] === 'string') {
        str += replacements[char];
      }
      else if (typeof replacements[char] === 'number') {
        str += String(replacements[char]);
      }
      else {
        str += String(replacements[char](date));
      }
    }
    else {
      str += char;
    }
  }

  return str;
}

export async function getFormattedFile(file: string) : Promise<ReadedFile> {
   
  const name = path.basename(file);
  const type = detectType(file.slice(file.indexOf('.') + 1));

  return {
    name,
    type,
    content: await FsPromise.readFile(file, 'utf-8'),
  };
}

function detectType(ext: string) {
  switch (ext) {
    case 'itp': return 'chemical/x-include-topology';
    case 'top': return 'chemical/x-topology';
    case 'pdb': return 'chemical/x-pdb';
  }
  return '';
}

export function rtrim(x: string, characters: string) {
  var start = 0;
  var end = x.length - 1;
  while (characters.indexOf(x[end]) >= 0) {
    end -= 1;
  }
  return x.substr(0, end + 1);
}

