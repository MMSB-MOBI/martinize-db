import logger from "./logger";
import { Molecule, BaseMolecule } from "./Entities/entities";
import { simpleflake } from 'simpleflakes';
import { ApiError } from "./Errors";
import Express from 'express';

export function isDebugMode() {
  return logger.level === "debug" || logger.level === "silly";
}

export function isMolecule(e: BaseMolecule) : e is Molecule {
  return 'approved_by' in e;
}

export function generateSnowflake() {
  return simpleflake().toString(10);
}

export function sendError(error: ApiError, res: Express.Response) {
  res.status(Number(error.message)).json(error.data);
}
