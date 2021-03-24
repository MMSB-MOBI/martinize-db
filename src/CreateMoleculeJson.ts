const fs = require('fs');
const dree = require('dree');

import { MoleculeChecker } from './routes/molecule/MoleculeChecker';
import { Molecule, BaseMolecule } from './Entities/entities';
import nano = require('nano');
import { Database } from './Entities/CouchHelper';
import { SimuRequest } from './routes/molecule/MoleculeChecker';
import { SimuFile } from './MoleculeOrganizer';
import { GoTerms } from './types';
import { ErrorType } from './Errors';

interface Version {
  number: string,
  itp: SimuFile,
  force_field: string
};

interface InfosJson {
  versions: Version[],
  name: string,
  alias: string,
  category: keyof typeof GoTerms,
  create_way: string,
  directory: string,
  top: {version: string, infos: SimuFile}[],
  gro: SimuFile
};

export const CreateMoleculeFromJson = async (jsonFile: string) => {

  let infos : InfosJson = JSON.parse(fs.readFileSync(jsonFile));
  let parentMol: string | null = null;

  let verCourante: Version = infos['versions'][0];
  
  let req : SimuRequest = {
    full_user: {
      id: '320054308425769119',
      role: "admin"
    },
    body: {
      name: infos['name'],
      alias: infos['alias'],
      smiles: '',
      version: '01',
      category: infos['category'],
      command_line: '',
      comments: '',
      create_way: infos['create_way'],
      force_field: 'martini22',                     // TODO Changer pour la bonne valeur
      validation: '',
      citation: '',
      parent: parentMol
    },
    files: {
      itp: [verCourante.itp,],
      pdb: [infos['gro'],],
      top: [infos.top[0]['infos'],],
      map: []
    }
  }

  const checker = new MoleculeChecker(req);
  try {
    await checker.checkName(req.body.name, '');
  } catch (e) {
    if (e.code === ErrorType.NameAlreadyExists) {
      parentMol = e.data.id;
      console.log(req);
      console.log(e.data);
      console.log(e.data.id);
    }
  }

  try {
    let molecule = await checker.check();
    parentMol = molecule.id;
  } catch (e) {
    console.log(e);
  }

  //console.log(parentMol);

  for (let i=1; i < infos['versions'].length; i++) {

    let ver : Version = infos['versions'][i];

    console.log(parentMol);

    req = {
      full_user: {
        id: '320054308425769119',
        role: "admin"
      },
      body: {
        name: infos['name'],
        alias: infos['alias'],
        smiles: '',
        version: ver.number,
        category: infos['category'],
        command_line: '',
        comments: '',
        create_way: infos['create_way'],
        force_field: 'martini22',                     // TODO Changer pour la bonne valeur
        validation: '',
        citation: '',
        parent: parentMol
      },
      files: {
        itp: [ver.itp,],
        pdb: [infos['gro'],],
        top: [infos.top[i]['infos'],],
        map: []
      }
    }

    try {
      let molecule = await checker.check();
      console.log(molecule);
    } catch (e) {
      console.log(e);
    }

  };
  //let response = await Database.molecule.save(molecule as Molecule);
  





  /*
  try {

    //console.log('lecture json');
    //let infos = JSON.parse(fs.readFileSync(jsonFile));
    //let parent : null | string;

    if(infos['versions'].length > 1) {

      console.log('entree if plusieurs versions');

      infos['versions'].forEach(async (version: Version) => {
        try{
          if(version.number === '01') {

            console.log('si version 01 :');


            let topversion = {originalname: '', path: '', size: 0};

            infos['top'].forEach((topfile: { version:string, infos: SimuFile}) => {
              if(topfile['version'] === '01'){
                topversion = topfile['infos'];
              }
            });

            let req : SimuRequest = {
              full_user: {
                id: '320054308425769119',
                role: "admin"
              },
              body: {
                name: infos['name'],
                alias: infos['alias'],
                smiles: '',
                version: '01',
                category: infos['category'],
                command_line: '',
                comments: '',
                create_way: infos['create_way'],
                force_field: 'martini22',                     // TODO Changer pour la bonne valeur
                validation: '',
                citation: '',
                parent: null
              },
              files: {
                itp: [version.itp,],
                pdb: [infos['gro'],],
                top: [topversion,],
                map: []
              }
            }

            console.log('requete 01 :', req);
            const checker = new MoleculeChecker(req);

            //var response: nano.DocumentInsertResponse;
            //var molecule: BaseMolecule;



            molecule = await checker.check();
            console.log('molecule : '+ molecule);
            parent = await molecule.parent;
            response = await Database.molecule.save(molecule as Molecule);

          }
          else {

            console.log('si autre version');

            let topversion = {originalname: '', path: '', size: 0};

            infos['top'].forEach((topfile: { version:string, infos: SimuFile}) => {
              if(topfile['version'] === version.number){
                topversion = topfile['infos'];
              }
            });

            let req : SimuRequest = {
              full_user: {
                id: '320054308425769119',
                role: "admin"
              },
              body: {
                name: infos['name'],
                alias: infos['alias'],
                smiles: '',
                version: version.number,
                category: infos['category'],
                command_line: '',
                comments: '',
                create_way: infos['create_way'],
                force_field: 'martini22',                     // TODO Changer pour la bonne valeur
                validation: '',
                citation: '',
                parent: parent
              },
              files: {
                itp: [version.itp,],
                pdb: [infos['gro'],],
                top: [topversion,],
                map: []
              }
            } 

            console.log('requete autre :', req);
            const checker = new MoleculeChecker(req);

            //var response: nano.DocumentInsertResponse;
            //var molecule: BaseMolecule;

            molecule = await checker.check();
            response = await Database.molecule.save(molecule as Molecule);
          }
        } catch (e) {
          throw(e);
        }

      });
    }
    else {
      //console.log('si une seule version');

      let req : SimuRequest = {
        full_user: {
          id: '320054308425769119',
          role: "admin"
        },
        body: {
          name: infos['name'],
          alias: infos['alias'],
          smiles: '',
          version: '01',
          category: infos['category'],
          command_line: '',
          comments: '',
          create_way: infos['create_way'],
          force_field: 'martini22',                     // TODO Changer pour la bonne valeur
          validation: '',
          citation: '',
          parent: null
        },
        files: {
          itp: [infos['versions'][0]['itp'],],
          pdb: [infos['gro'],],
          top: [infos['top'][0]['infos'],],
          map: []
        }
      }

      //console.log('requete autre :', req);
      const checker = new MoleculeChecker(req);

      //var response: nano.DocumentInsertResponse;
      //var molecule: BaseMolecule;

      molecule = await checker.check();
      response = await Database.molecule.save(molecule as Molecule);
    }

  } catch (e) {
  console.log("une erreur s'est produite");
  console.log(e);
  return ;
  }
  */
}


const options = {
  depth:10,
  extensions: ['json'],
};
const tree = dree.scan('/home/achopin/Documents/database/martini-molecule-repository/martini2_lipids_test/Glycosphingolipids', options, (element : any) => {
  console.log(element.name);
  CreateMoleculeFromJson(element.path);
});

