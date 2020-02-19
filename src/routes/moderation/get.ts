import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed, escapeRegExp } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import Errors, { ErrorType } from '../../Errors';
import { StashedMolecule } from '../../Entities/entities';

const GetStashedRouter = Router();

GetStashedRouter.get('/', (req, res) => {
  (async () => {
    const { version } = req.query;

    const selector: any = {};

    if (version) {
      selector.id = version;
    }

    const query: nano.MangoQuery = { selector };

    const desired_molecule = await Database.stashed.find(query);

    if (!desired_molecule.length) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const molecule = desired_molecule[0];

    // Find parent if any
    const parent = molecule.parent ? await Database.molecule.find({ selector: { parent: molecule.parent }, limit: 1 }) : undefined;

    let good_version: StashedMolecule = molecule;
  
    res.json({
      molecule: sanitize(good_version),
      parent: parent ? sanitize(parent) : undefined
    });
  })().catch(errorCatcher(res));
});

GetStashedRouter.all('/', methodNotAllowed('GET'))

export default GetStashedRouter;
