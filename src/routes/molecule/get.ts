import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed, escapeRegExp } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import Errors, { ErrorType } from '../../Errors';
import { Molecule } from '../../Entities/entities';

const GetMoleculeRouter = Router();

GetMoleculeRouter.get('/', (req, res) => {
  (async () => {
    const { alias, version, tag } = req.query;
    
    if (!alias) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const selector: any = { alias: { $regex: '(?i)^' + escapeRegExp(alias as string) + '$' } };
    let search_specific_version = false;

    if (version) {
      selector.id = version;
      search_specific_version = true;
    }
    if (tag) {
      selector.version = tag;
      search_specific_version = true;
    }

    const query: nano.MangoQuery = { selector };

    const desired_molecule = await Database.molecule.find(query);

    if (!desired_molecule.length) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const molecule = desired_molecule[0];
    const tree_id = molecule.tree_id;

    // Find all versions
    const versions = await Database.molecule.find({ selector: { tree_id }, limit: 9999999 });

    let good_version: Molecule = molecule;
    
    if (!search_specific_version) {
      // Find the last version of the molecule
      good_version =  versions.reduce((prev, actual) => {
        const prev_date = new Date(prev.created_at);
        const actual_date = new Date(actual.created_at);
  
        if (actual_date.getTime() > prev_date.getTime()) {
          // Keep the most recent
          return actual;
        }
        return prev;
      });
    }

    res.json({
      molecule: sanitize(good_version),
      versions: versions.map(sanitize),
    });
  })().catch(errorCatcher(res));
});

GetMoleculeRouter.all('/', methodNotAllowed('GET'))

export default GetMoleculeRouter;
