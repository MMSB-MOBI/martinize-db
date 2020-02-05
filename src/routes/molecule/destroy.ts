import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';

const DestroyMoleculeRouter = Router();

DestroyMoleculeRouter.delete('/:id', (req, res) => {
  (async () => {
    // TODO MAKE SEARCH
    // For now, it only returns every molecule
    const id = req.params.id;

    if (!id || typeof id !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    try {
      const mol = await Database.molecule.get(id);

      // Delete attached ZIP
      await MoleculeOrganizer.remove(mol.files);
      await Database.molecule.delete(mol);
    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    res.send();
  })().catch(errorCatcher(res));
});

DestroyMoleculeRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyMoleculeRouter;
