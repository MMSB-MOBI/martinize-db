const fs = require('fs');
const dree = require('dree');

import { MoleculeChecker } from './MoleculeChecker';
import { Molecule, BaseMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import { GoTerms } from '../../types';
import Errors, { ErrorType } from '../../Errors';
import logger from '../../logger';
import { resolve } from 'path';
import { CONNECTED_USER_CLI } from '../../cli/user_cli';
import { Excel } from '../../cli/molecule_cli';





/**
 * One itp file version written in the Json file
 */
 interface VersionItp {
  number: string,
  itp: SimuFile,
  force_field: string,
  protonation? : string
};

/**
 * All the informations on the molecule stored in the Json file used to simulate a request
 */
export interface InfosJson {
  versions: VersionItp[],
  /* Name of the molecule */
  name: string,
  alias: string,
  category: keyof typeof GoTerms,
  create_way: string,
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
    role: string
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
export const CreateMoleculeFromJson = async (batch : InfosJson[]) => {

  await batch.reduce(async (memo, i) => {
    await memo;
    logger.info('Inserting '+i.name);
    await CreateMoleculeFromJsonAux(i);

  }, Promise.resolve());
}

/*
const options = {
  depth:10,
  extensions: ['json'],
};
const tree = dree.scan('/home/achopin/Documents/database/martini-molecule-repository/martini2_lipids_test/Glycosphingolipids', options, (element : any) => {
  CreateMoleculeFromJson(element.path);
});

*/

const CreateMoleculeFromJsonAux = async (infos : InfosJson) => {

  let parentMol: string | null = null;

    for (let i=0; i < infos.versions.length; i++) {

      let ver : VersionItp = infos.versions[i];

      let martiniVer = '';
      if (ver.force_field == 'v2.0' || ver.force_field == 'v2.1' || ver.force_field == 'v2.2') {
        martiniVer = 'martini22';
      }

      let req : SimuRequest = {
        full_user: {
          id: CONNECTED_USER_CLI.id,
          role: CONNECTED_USER_CLI.role
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


      if (req.full_user.role != 'admin') {
        return Errors.throw(ErrorType.Unallowed);
      }

      if (!infos.top[i].infos.originalname.match(ver.number)){
        return Errors.throw(ErrorType.MissingTopFiles);
      }


      if (ver.number == '01') {
        const checker = new MoleculeChecker(req);
        try {
          await checker.checkName(req.body.name, '');
        } catch (e) {
          if (e.code === ErrorType.NameAlreadyExists) {
            logger.warn("There is already a molecule by the name "+infos.name+" in the database.");
            break;
          }
          else {
            return e;
          }
        }
      }

      try {
        const checker = new MoleculeChecker(req);
        let molecule = await checker.check();
        if (ver.number == '01') {
          parentMol = molecule.id;
        }
        let response = await Database.molecule.save(molecule as Molecule);
        Excel.text += infos.name+',,,,X\n';
      } catch (e) {
          if (e.code !== ErrorType.NameAlreadyExists) {
            logger.warn("Error with the "+ver.number+" version of \""+infos.name+"\" : "+ e.data.message);
            if (e.code == ErrorType.InvalidMoleculeFiles) {
              Excel.text += infos.name+',,X,,\n';
            }
            if (req.body.force_field == '') {
              Excel.text += infos.name+',,,X,\n';
            }
            break;
          };
      }

    };
    return new Promise((resolve, reject) => {resolve(true)});
}

