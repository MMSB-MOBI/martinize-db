import logger from "./logger";
import { Molecule, BaseMolecule } from "./Entities/entities";
import { simpleflake } from 'simpleflakes';
import Errors, { ApiError, ErrorType } from "./Errors";
import Express from 'express';
import { TokenPayload } from "./types";
import JsonWebToken from 'jsonwebtoken';
import { KEYS, UPLOAD_ROOT_DIR } from "./constants";

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

/**
 * Create a closure to catch an ApiError that occurs in a Promise.
 */
export function errorCatcher(res: Express.Response) {
  return function(err: any) {
    if (res.headersSent) {
      return;
    }

    if (err instanceof ApiError) {
      return sendError(err, res);
    }

    logger.error("Unknown error: ", err);
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
