import MoleculeOrganizer, { MoleculeSave } from './MoleculeOrganizer';
import * as fs from 'fs';

import { MoleculeChecker } from './routes/molecule/MoleculeChecker';
import { Molecule, BaseMolecule } from './Entities/entities';
import nano = require('nano');
import { Database } from './Entities/CouchHelper';
import { SimuRequest } from './routes/molecule/MoleculeChecker';

export const test = async (req: SimuRequest) => {
    //try {

    try {
        const checker = new MoleculeChecker(req);
    

        let response: nano.DocumentInsertResponse;
        let molecule: BaseMolecule;


        molecule = await checker.check();
        response = await Database.molecule.save(molecule as Molecule);
    
        } catch (e) {
        console.log("une erreur s'est produite");
        console.log(e);
        return ;
        }
   // }

        /*
        let save = await MoleculeOrganizer.save(
            itpfiles, 
            pdbfiles, 
            {originalname: "test.top", path: "/home/achopin/Documents/database/martini-molecule-repository/martini2_lipids/LPC/APC/test.top", size: 140}, 
            [],
            "martini22"
        );
        console.log("ca a fonctionne");
    } catch (e) {
        console.log(e);
    }*/
}


let itpfiles = JSON.parse(fs.readFileSync("/home/achopin/Documents/martinize/martinize-db/src/testitp2.json", 'utf8'));
let pdbfiles = JSON.parse(fs.readFileSync("/home/achopin/Documents/martinize/martinize-db/src/testpdb.json", 'utf8'));
let topfiles = JSON.parse(fs.readFileSync("/home/achopin/Documents/martinize/martinize-db/src/testtop.json", 'utf8'));

let r : SimuRequest = {
    full_user: {
      id: '320054308425769119',
      role: "admin"
    },
    body: {
      name: '1',
      alias: '1',
      smiles: '',
      version: '3.0',
      category: 'GO:0001',
      command_line: '',
      comments: '',
      create_way: 'martinize_2',
      force_field: 'elnedyn22p',
      validation: '',
      citation: '',
      parent: "323729664604577506",
    },
    files: {
      itp: itpfiles,
      pdb: pdbfiles,
      top: topfiles,
      map: []
    }
  };

let t = test(r);