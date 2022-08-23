import nano from "nano";
import { StashedMolecule } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import CompleteMoleculeVersionTree from "./MoleculeVersionTree";

export default class MoleculeDatabase extends AbstractDatabase<StashedMolecule> {

    async moleculeTreeOf(tree_id: string) {
        const mols = await this.find({ limit: 999999, selector: { tree_id } });
    
        if (mols.length) {
          return new CompleteMoleculeVersionTree(mols);
        }
      }

}
