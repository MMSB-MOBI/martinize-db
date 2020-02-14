import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed, escapeRegExp } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import Errors, { ErrorType } from '../../Errors';

const GetMoleculeRouter = Router();

GetMoleculeRouter.get('/', (req, res) => {
  (async () => {
    const alias = req.query.alias;
    
    if (!alias) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const query: nano.MangoQuery = { 
      selector: { alias: { $regex: '(?i)^' + escapeRegExp(alias) + '$' } }
    };

    const desired_molecule = await Database.molecule.find(query);

    if (!desired_molecule.length) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const molecule = desired_molecule[0];
    const tree_id = molecule.tree_id;

    // Find all versions
    const versions = await Database.molecule.find({ selector: { tree_id }, limit: 9999999 });

    // Find the last version of the molecule
    const last_version = versions.reduce((prev, actual) => {
      const prev_date = new Date(prev.created_at);
      const actual_date = new Date(actual.created_at);

      if (actual_date.getTime() > prev_date.getTime()) {
        // Keep the most recent
        return actual;
      }
      return prev;
    });

    res.json({
      molecule: sanitize(last_version),
      versions: versions.map(sanitize),
    });
  })().catch(errorCatcher(res));
});

GetMoleculeRouter.all('/', methodNotAllowed('GET'))

export default GetMoleculeRouter;
