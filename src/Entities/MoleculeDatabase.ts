import nano from "nano";
import { Molecule } from "./entities";

export default class MoleculeDatabase {
  constructor(protected _db: nano.DocumentScope<Molecule>) {}

  find(query: nano.MangoQuery) {
    return this._db.find(query);
  }

  get db() {
    return this._db;
  }
}
