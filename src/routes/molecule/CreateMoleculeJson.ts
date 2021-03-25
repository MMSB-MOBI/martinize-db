const fs = require('fs');
const dree = require('dree');

import { MoleculeChecker } from './MoleculeChecker';
import { Molecule, BaseMolecule } from '../../Entities/entities';
//import nano = require('nano');
import { Database } from '../../Entities/CouchHelper';
import { GoTerms } from '../../types';
import { ErrorType } from '../../Errors';
//import logger from './logger';

/**
 * One itp file version written in the Json file
 */
interface VersionItp {
  number: string,
  itp: SimuFile,
  force_field: string
};

/**
 * All the informations on the molecule stored in the Json file used to simulate a request
 */
interface InfosJson {
  versions: VersionItp[],
  /* Name of the molecule */
  name: string,
  alias: string,
  category: keyof typeof GoTerms,
  create_way: string,
  force_field: string,
  directory: string,
  top: {version: string, infos: SimuFile}[],
  map: SimuFile[],
  gro: SimuFile
};

/**
 * Simulate the files informations contained in the request sent to the server
 */
export interface SimuFile {
  originalname: string;
  path: string;
  size: number
}

/**
 * Simulate the informations contained in a request sent to the server to add a molecule to the database
 */
export interface SimuRequest{
  full_user: {
    id: string,
    role: "admin"
  };
  body: {
    name: string,
    alias: string,
    smiles: string,
    version: string,
    category: keyof typeof GoTerms,
    command_line: string,
    comments: string,
    create_way: string,
    force_field: string,
    validation: string,
    citation: string,
    parent: string | null
  };
  files: {
    itp: SimuFile[],
    pdb: SimuFile[],
    top: SimuFile[] | [],
    map: SimuFile[] | []
  };
}


/**
 * Create a molecule object from its data contained in a Json file and send it to the database
 * @param jsonFile 
 */
export const CreateMoleculeFromJson = async (jsonFile: string) => {

  let infos : InfosJson = JSON.parse(fs.readFileSync(jsonFile));
  let parentMol: string | null = null;

  let verCourante: VersionItp = infos.versions[0];

  let martiniVer = '';
  if (verCourante.force_field == 'v2.0' || verCourante.force_field == 'v2.1' || verCourante.force_field == 'v2.2') {
    martiniVer = 'martini22';
  }

  
  let req : SimuRequest = {
    full_user: {
      id: '320054308425769119',
      role: "admin"
    },
    body: {
      name: infos.name,
      alias: infos.alias,
      smiles: '',
      version: '01',
      category: infos.category,
      command_line: '',
      comments: '',
      create_way: infos.create_way,
      force_field: martiniVer,
      validation: '',
      citation: '',
      parent: parentMol
    },
    files: {
      itp: [verCourante.itp,],
      pdb: [infos.gro,],
      top: [infos.top[0].infos,],
      map: infos.map
    }
  }

  console.log(req);

  const checker = new MoleculeChecker(req);
  try {
    await checker.checkName(req.body.name, '');
  } catch (e) {
    if (e.code === ErrorType.NameAlreadyExists) {
      parentMol = e.data.id;
      console.log("There is already a molecule by the name "+infos.name+" in the database. If the other versions numbers are not already in the database, the molecules will be added as other versions.\n");
    }
  }
  

  try {
    let molecule = await checker.check();
    parentMol = molecule.id;
    let response = await Database.molecule.save(molecule as Molecule);
  } catch (e) {
    if (e.code != ErrorType.NameAlreadyExists) {
      console.log(e.data.message);
    }
  }

  //console.log(parentMol);

  for (let i=1; i < infos.versions.length; i++) {

    let ver : VersionItp = infos.versions[i];

    if (ver.force_field == 'v2.0' || ver.force_field == 'v2.1' || ver.force_field == 'v2.2') {
      martiniVer = 'martini22';
    }

    req = {
      full_user: {
        id: '320054308425769119',
        role: "admin"
      },
      body: {
        name: infos.name,
        alias: infos.alias,
        smiles: '',
        version: ver.number,
        category: infos.category,
        command_line: '',
        comments: '',
        create_way: infos.create_way,
        force_field: martiniVer,
        validation: '',
        citation: '',
        parent: parentMol
      },
      files: {
        itp: [ver.itp,],
        pdb: [infos.gro,],
        top: [infos.top[i].infos,],
        map: infos.map
      }
    }

    try {
      const checker = new MoleculeChecker(req);
      let molecule = await checker.check();
      let response = await Database.molecule.save(molecule as Molecule);
    } catch (e) {
      console.log("Error with the "+ver.number+" version : "+ e.data.message);
    }

  };
}


const options = {
  depth:10,
  extensions: ['json'],
};
const tree = dree.scan('/home/achopin/Documents/database/martini-molecule-repository/martini2_lipids_test/Glycosphingolipids', options, (element : any) => {
  CreateMoleculeFromJson(element.path);
});

