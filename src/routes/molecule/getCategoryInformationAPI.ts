import { Router } from 'express';
import { errorCatcher } from '../../helpers'; 
import { Database } from '../../Entities/CouchHelper'; 
import { MangoQuery } from 'nano';
import {decodeCategory } from "./parser/parser_files";

// Get a pdb from a file ID
const GetCategoryInformationAPI = Router();


//If format isnt provided give the last update of this model 
GetCategoryInformationAPI.get('/:category', (req, res) => {
  (async () => {
    console.log("GetcategoryInformationAPI.get() ", req.params)
    const selectruc: MangoQuery = { selector: { force_field: req.params.category } }
    const molcouch = await Database.molecule.find(selectruc)

    // File does not exists
    if (molcouch.length === 0) {
      res.send({ "Error": "category not found" });
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
          "forcefield": molcouch[i].force_field,
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
        "category": decodeCategory(reponse[i].category[0]),
        "forcefield": reponse[i].forcefield,
        "version" : reponse[i].version
      }
      c++
    }

    res.send(reponse2);
  })().catch(errorCatcher(res));
});

export default GetCategoryInformationAPI;