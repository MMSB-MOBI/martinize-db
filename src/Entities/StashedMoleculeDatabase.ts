import nano from "nano";
import { StashedMolecule } from "./entities";

export default class MoleculeDatabase {
  constructor(protected _db: nano.DocumentScope<StashedMolecule>) {}

  async find(query: nano.MangoQuery) {
    const res = await this._db.find(query);
    return res.docs;
  }

  get db() {
    return this._db;
  }
}
