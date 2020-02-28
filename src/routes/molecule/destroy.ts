import { Router } from 'express';
import { errorCatcher, methodNotAllowed, deleteMolecule } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import SearchWorker from '../../search_worker';

const DestroyMoleculeRouter = Router();

DestroyMoleculeRouter.delete('/:id', (req, res) => {
  (async () => {
    // TODO MAKE SEARCH
    // For now, it only returns every molecule
    const id = req.params.id;

    if (!id || typeof id !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const user = req.full_user!;

    if (!user) {
      return Errors.throw(ErrorType.Forbidden);
    }

    await deleteMolecule(id, user, false, true);

    res.send();
  })().catch(errorCatcher(res));
});

DestroyMoleculeRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyMoleculeRouter;
