import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { BaseMolecule } from '../../Entities/entities';
import JSZip from 'jszip';
import { MangoQuery } from 'nano';
import {decodeCategory } from "../../routes/molecule/parser/parser_files";

// Get a pdb from a file ID
const GetInformationAPI = Router();


//If format isnt provided give the last update of this model 
GetInformationAPI.get('/:forcefield', (req, res) => {
  (async () => {
    console.log("GetInformationAPI.get() ", req.params)
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
          "forcefield": molcouch[i].force_field,
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
        "forcefield": reponse[i].forcefield,
        "category": decodeCategory(reponse[i].category[0]),
        "version" : reponse[i].version
      }
      c++
    }

    res.send(reponse2);
  })().catch(errorCatcher(res));
});

export default GetInformationAPI;