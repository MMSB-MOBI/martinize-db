import logger from "./logger";
import { Molecule, BaseMolecule, User } from "./Entities/entities";
import { simpleflake } from 'simpleflakes';
import Errors, { ApiError, ErrorType } from "./Errors";
import Express from 'express';
import { TokenPayload } from "./types";
import JsonWebToken from 'jsonwebtoken';
import { KEYS, UPLOAD_ROOT_DIR } from "./constants";
import { Database } from "./Entities/CouchHelper";
import { unlink } from "fs";
import MoleculeOrganizer from "./MoleculeOrganizer";
import SearchWorker from "./search_worker";
import Mailer from "./Mailer/Mailer";

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
    else if (err) {
      logger.error("Unknown error: ", JSON.stringify(Object.getOwnPropertyDescriptors(err)));
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
        issuer: "Martinize Database Server 1", 
        jwtid: id, // ID généré avec snowflake
      }, 
      (err, encoded) => { // Quand le token est généré (ou non), accepte/rejette la promesse
        if (err) reject(err);
        else resolve(encoded);
      }
    );
  }) as Promise<string>;
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
  if (stashed) {
    if (user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    try {
      const mol = await Database.stashed.get(id);

      // Delete attached ZIP
      await MoleculeOrganizer.remove(mol.files);
      await Database.stashed.delete(mol);
    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }
    return;
  }

  // Delete a published molecule
  try {
    const mol = await Database.molecule.get(id);

    if (user.role !== "admin" && user.id !== mol.owner) {
      return Errors.throw(ErrorType.Forbidden);
    }

    if (checked_attached) {
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
    await Database.molecule.delete(mol);
  } catch (e) {
    return Errors.throw(ErrorType.ElementNotFound);
  }
}

export async function informAdminFromAskCreation(new_user: User) {
  const admins = await Database.user.find({
    selector: { role: "admin" },
    limit: 99999
  });

  for (const admin of admins) {
    // Send a mail
    await Mailer.send({
      to: admin.email,
      subject: "MArtinize Database - New account request",
    }, "mail_ask", {
      new_user: {
        name: new_user.name
      },
      name: admin.name
    });
  }
}
