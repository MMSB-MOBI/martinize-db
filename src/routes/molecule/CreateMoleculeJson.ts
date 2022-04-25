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
import { info } from 'console';
import { create } from 'domain';





/**
 * One itp file version written in the Json file
 */
 interface VersionItp {
  number: string,
  itp: SimuFile,
  force_field: string,
  protonation? : string,
  comments: string,
  citation: string,
};

/**
 * All the informations on the molecule stored in the Json file used to simulate a request
 */
export interface InfosJson {
  versions: VersionItp[],
  /* Name of the molecule */
  name: string,
  alias: string,
  category: string[],
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
    category: string[],
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

interface InsertionRecap {
  inserted : string[],
  not_inserted : {[reason: string]: string[]}
}


/**
 * Read a bacth of molecule objects contained in a Json object and send their infos to the database to insert them sequentially
 * @param batch - a list of molecules informations
 */
export const CreateMoleculeFromJson = async (batch : InfosJson[]) : Promise< InsertionRecap > => {

  const recap : InsertionRecap = {'inserted': [], 'not_inserted' : {'other': []}}
  for (const info of batch) {
    logger.info(`Try to insert ${info.name}`)
    try {
      await CreateMoleculeFromJsonAux(info)
      console.log("Inserted")
      recap.inserted.push(info.directory)
    } catch(e) {
      console.log("Not inserted")
      if(e.data && e.data.message){
        if(!(e.data.message in recap.not_inserted)) recap.not_inserted[e.data.message] = []
        recap.not_inserted[e.data.message].push(info.directory)
      }
      else {
        recap.not_inserted['other'].push(info.directory)
      }
      
    }
  }
  return new Promise((res, rej) => res(recap))
}


/**
 * Auxiliary function reading the infos from one molecule of the batch and inserting them in the databse
 * @param infos - The informations on the moelcule
 * @returns a Promise with resolve true :)
 */
const CreateMoleculeFromJsonAux = async (infos : InfosJson) => {

  let parentMol: string | null = null;

    // for each itp versions (versions of the molecule) we insert the infos in the database
    for (let i=0; i < infos.versions.length; i++) {

      let ver : VersionItp = infos.versions[i];

      const martiniVer = ver.force_field;

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
          comments: ver.comments,
          create_way: infos.create_way,
          force_field: martiniVer,
          validation: '',
          citation: ver.citation,
          parent: parentMol
        },
        files: {
          itp: [ver.itp,],
          pdb: [infos.gro,],
          top: [infos.top[i].infos,],
          map: infos.map
        }
      }

      // Error if the user connected is not an admin
      if (req.full_user.role != 'admin') {
        return Errors.throw(ErrorType.Unallowed);
      }

      // If there are no top file correspnding to the itp
      if (infos.versions[i].number !== infos.top[i].version){
        return Errors.throw(ErrorType.MissingTopFiles);
      }

      const checker = new MoleculeChecker(req);
      await checker.checkName(req.body.name, '');
      const molecule = await checker.check()
      await Database.molecule.save(molecule as Molecule);

    };
    return new Promise((resolve, reject) => {resolve(true)});
}