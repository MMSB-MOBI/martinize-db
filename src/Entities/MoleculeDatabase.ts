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
    const aliased: {[alias: string]: Molecule[]} = {}
    let byCategories: any = {}
    let byForceField: any = {}
    for (const mol of mols){
      if(!(mol.alias in aliased)) aliased[mol.alias] = []
      aliased[mol.alias].push(mol)
    }     
    for (const alias in aliased){
      const mols = aliased[alias]
      const ffs = new Set(mols.map(mol => mol.force_field))
      const cats = new Set(mols.map(mol => mol.category).flat())
      if(alias === "ALA"){
        console.log("FFS", ffs)
        console.log("CATS", cats)
      }
      for (const ff of ffs){
        if(!(ff in byForceField)) byForceField[ff] = 0
        byForceField[ff] += 1
      }
      for (const cat of cats){
        if(!(cat in byCategories)) byCategories[cat] = 0
        byCategories[cat] += 1
      }
    }
    return({byCategories, byForceField, all:aliased})
  }
}
