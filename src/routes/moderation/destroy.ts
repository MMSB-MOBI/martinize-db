import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';

const DestroyStashedRouter = Router();

DestroyStashedRouter.delete('/:id', (req, res) => {
  (async () => {
    const id = req.params.id;

    if (!id || typeof id !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const user = req.full_user!;

    if (!user) {
      return Errors.throw(ErrorType.Forbidden);
    }

    if (user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    try {
      const mol = await Database.stashed.get(id);

      // Delete attached ZIP
      await MoleculeOrganizer.remove(mol.files);
      await Database.stashed.delete(mol);
    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    res.send();
  })().catch(errorCatcher(res));
});

DestroyStashedRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyStashedRouter;
