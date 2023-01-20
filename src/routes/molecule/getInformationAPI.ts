import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import { MangoQuery } from 'nano';
import { decodeCategory } from "./parser/parser_files";
import { AvailableForceFields, FORCE_FIELDS } from "../../routes/types";
import settings from '../settings';
// Get a pdb from a file ID
const GetInformationAPI = Router();

function getInformation(field: string, value: string) {
  const selectruc: MangoQuery = { selector: { [field]: value } }
  return Database.molecule.find(selectruc)
    .then(molcouch => {
      if (molcouch.length === 0) {
        return { error: "No molecule found for this forcefield." }
         
      }
      let response: any = {};
      for (let i in Object.keys(molcouch)) {
        const alias = molcouch[i].alias;
        if (Object.keys(response).includes(alias)) {
          response[alias]["version"].push(molcouch[i].version);
        }
        else {
          response[alias] = {
            "name": molcouch[i].name,
            "citation": molcouch[i].citation,
            "forcefield": molcouch[i].force_field,
            "category": molcouch[i].category,
            "version": [molcouch[i].version],
          };
        }
      }
      let response2: any = {};
      let c = 0;
      for (let i of Object.keys(response)) {
        response2[c] = {
          "alias": i,
          "name": response[i].name,
          "citation": response[i].citation,
          "forcefield": response[i].forcefield,
          "category": decodeCategory(response[i].category[0]),
          "version": response[i].version
        };
        c++;
      }
      return response2;
    });
}

//If format isnt provided give the last update of this model 
GetInformationAPI.get('/:field', (req, res) => {
  (async () => {
    const field = req.params.field as AvailableForceFields
    if (FORCE_FIELDS.includes(field)) {
      console.log( field)
      const response = await getInformation("force_field", field);
      res.send(response);

    }
    else {
      res.status(400).send({ error: "Invalid field" });
    }

    // else #CHECK if it's a ctagory field {
    //   const response = await getInformation("category", field);
    //   res.send(response);
    // }

  })().catch(errorCatcher(res));
});

export default GetInformationAPI;