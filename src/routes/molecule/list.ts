import { Router } from 'express';
import { errorCatcher, sanitize } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';

const ListMoleculeRouter = Router();

ListMoleculeRouter.get('/', (req, res) => {
  (async () => {
    // TODO MAKE SEARCH
    // For now, it only returns every molecule

    const all_accepted_mols = await Database.molecule.all();

    res.json(all_accepted_mols.map(m => sanitize(m)));
  })().catch(errorCatcher(res));
});

export default ListMoleculeRouter;
