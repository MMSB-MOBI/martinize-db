import logger from "./logger";
import { Molecule, BaseMolecule } from "./Entities/entities";
import { simpleflake } from 'simpleflakes';
import Errors, { ApiError, ErrorType } from "./Errors";
import Express from 'express';
import { TokenPayload } from "./types";
import JsonWebToken from 'jsonwebtoken';
import { KEYS, UPLOAD_ROOT_DIR } from "./constants";
import { Database } from "./Entities/CouchHelper";
import { unlink } from "fs";
import MoleculeOrganizer from "./MoleculeOrganizer";

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
      logger.error("During request handling, the following error occured: " + err + "\n" + err.stack);
    }
    else {
      logger.error("Unknown error: ", err);
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

export async function verifyAndCompleteMolecule(molecule: BaseMolecule, edit = false, stashed = false) {
  // Must check every thing: if field exists
  if ([
    molecule.id,
    molecule.name,
    molecule.alias,
    molecule.formula,
    molecule.version,
    molecule.category,
    molecule.command_line,
    molecule.martinize_version,
    molecule.force_field
  ].some(e => typeof e !== 'string')) {
    throw { error: 'Required parameter is missing' };
  }

  // If parent is defined, check if it exists
  let parent: Molecule | undefined;
  if (molecule.parent && typeof molecule.parent === 'string') {
    try {
      parent = await Database.molecule.get(molecule.parent);
    } catch (e) { 
      throw { error: 'Failed to get parent', detail: e };
    }
  }
  else {
    molecule.parent = null;
  }

  if (parent) {
    molecule.tree_id = parent.tree_id;
    molecule.name = parent.name;
    molecule.alias = parent.alias;
    molecule.formula = parent.formula;
    molecule.category = parent.category;
  }
  
  if (edit) {
    if (!molecule.tree_id) {
      throw { error: 'When you edit a molecule, it should have a tree ID' };
    }
  }
  else {
    molecule.tree_id = generateSnowflake();
  }

  if (edit && !stashed) {
    (molecule as Molecule).last_update = new Date().toISOString();
  }

  // TODO CHECK EVERYTHING ELSE except hash & files parameters, 
  // because file isn't necessary saved
  if (molecule.files) {
    // Check if file exists
    const exists = await MoleculeOrganizer.exists(molecule.files);
    if (!exists) {
      throw { error: "Specified file not found. Has it been deleted ?" };
    }

    molecule.hash = await MoleculeOrganizer.hash(molecule.files);
  }


  return molecule;
}

export function escapeRegExp(string: string) {
  return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); // $& means the whole matched string
}

export function withRegex(text: string, is_regex: boolean, flags = "i") {
  const search_text = is_regex ? text : escapeRegExp(text);
  return { $regex: `(?${flags})${search_text}` };
}
