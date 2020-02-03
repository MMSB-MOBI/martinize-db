import nano = require("nano");
import { Molecule } from "./entities";

export class MoleculeDatabase {
  constructor(protected _db: nano.DocumentScope<Molecule>) {}

  

  get db() {
    return this._db;
  }
}
