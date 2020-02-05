import { Molecule } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import MoleculeVersionTree from "./MoleculeVersionTree";

export default class MoleculeDatabase extends AbstractDatabase<Molecule> {
  /**
   * Construct the `MoleculeVersionTree` of {tree_id}
   *
   * If any molecule have this tree id, returns `undefined`
   */
  async moleculeTreeOf(tree_id: string) {
    const mols = await this.find({ limit: 999, selector: { tree_id } });

    if (mols.length) {
      return new MoleculeVersionTree(mols);
    }
  }
}
