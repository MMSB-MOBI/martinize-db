import { Molecule } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import MoleculeVersionTree from "./MoleculeVersionTree";
import nano = require("nano");
import SearchWorker from "../search_worker";

export default class MoleculeDatabase extends AbstractDatabase<Molecule> {
  /**
   * Construct the `MoleculeVersionTree` of {tree_id}
   *
   * If any molecule have this tree id, returns `undefined`
   */
  async moleculeTreeOf(tree_id: string) {
    const mols = await this.find({ limit: 999999, selector: { tree_id } });

    if (mols.length) {
      return new MoleculeVersionTree(mols);
    }
  }

  /** Bulk create */
  bulkCreate(docs: Molecule[]) {
    SearchWorker.clearCache();
    return super.bulkCreate(docs);
  }

  /** Bulk update */
  bulkUpdate(docs: Molecule[]) {
    SearchWorker.clearCache();
    return super.bulkUpdate(docs);
  }

  /** Bulk delete */
  bulkDelete(documents: Molecule[]) {
    SearchWorker.clearCache();
    return super.bulkDelete(documents);
  }

  /** Bulk create, update or delete */
  bulk(docs: nano.BulkModifyDocsWrapper) {
    SearchWorker.clearCache();
    return super.bulk(docs);
  }

  save(element: Molecule) {
    SearchWorker.clearCache();
    return super.save(element);
  }

  delete(element: Molecule) {
    SearchWorker.clearCache();
    return super.delete(element);
  }
}
