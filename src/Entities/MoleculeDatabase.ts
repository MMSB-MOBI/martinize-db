import { Molecule } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import CompleteMoleculeVersionTree from "./MoleculeVersionTree";
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
      return new CompleteMoleculeVersionTree(mols);
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

  async getVersionsFromAlias(alias: string) {
    const mols = await this.find({ limit: 999999, selector: { alias } });
    return mols
  }

  async stats(){
    SearchWorker.clearCache();
    const mols = await this.all()
    let byCategories: any = {}
    let byForceField: any = {}
    for (const mol of mols){
      if(!(mol.category in byCategories)) byCategories[mol.category] = []
      if(!(mol.force_field in byForceField)) byForceField[mol.force_field] = []
      byCategories[mol.category].push(mol)
      byForceField[mol.force_field].push(mol)
    }
    return({byCategories, byForceField, all:mols})
  }
}
