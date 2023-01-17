import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { BaseMolecule } from '../../Entities/entities';
import JSZip from 'jszip';
import { MangoQuery } from 'nano';

// Get a pdb from a file ID
const GetForcefieldInformationAPI = Router();


//If format isnt provided give the last update of this model 
GetForcefieldInformationAPI.get('/:forcefield', (req, res) => {
  (async () => {
    console.log("GetForcefieldInformationAPI.get() ", req.params)
    const selectruc: MangoQuery = { selector: { force_field: req.params.forcefield } }
    const molcouch = await Database.molecule.find(selectruc)

    // File does not exists
    if (molcouch.length === 0) {
      res.send({ "Error": "forcefield not found" });
    } 
    let reponse: any = {}

    
    for (let i in Object.keys(molcouch)) {
      const alias: string = molcouch[i].alias
      if (Object.keys(reponse).includes(alias)) {
        reponse[alias]["version"].push(molcouch[i].version)
      }
      else {
        reponse[alias] = {
          "name": molcouch[i].name,
          "citation" : molcouch[i].citation,
          "category": molcouch[i].category,
          "version": [molcouch[i].version],
        }
      }

    }

    //transform le dico
    let reponse2 : any = {}
    let c = 0
    for (let i of Object.keys(reponse)) {
      reponse2[c] = {
        "alias" : i,
        "name" : reponse[i].name,
        "citation" : reponse[i].citation,
        "category": reponse[i].category,
        "version" : reponse[i].version
      }
      c++
    }

    res.send(reponse2);
  })().catch(errorCatcher(res));
});

export default GetForcefieldInformationAPI;