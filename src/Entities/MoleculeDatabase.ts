import nano from "nano";
import { Molecule } from "./entities";

export default class MoleculeDatabase {
  constructor(protected _db: nano.DocumentScope<Molecule>) {}

  async find(query: nano.MangoQuery) {
    const res = await this._db.find(query);
    return res.docs;
  }

  get db() {
    return this._db;
  }
}
