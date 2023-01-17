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

    console.log(molcouch)
    let reponse: any = {}

    for (let i in Object.keys(molcouch)) {
      const alias: string = molcouch[i].alias
      if (Object.keys(reponse).includes(alias)) {
        reponse[alias]["version"].push(molcouch[i].version)
      }
      else {
        reponse[alias] = {
          "name": molcouch[i].name,
          "version": [molcouch[i].version],
        }
      }

    }
    res.send(reponse);
  })().catch(errorCatcher(res));
});

export default GetForcefieldInformationAPI;