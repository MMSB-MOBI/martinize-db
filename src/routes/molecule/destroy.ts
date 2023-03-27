import { Router } from 'express';
import { errorCatcher, methodNotAllowed, deleteMolecule } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';

const DestroyMoleculeRouter = Router();

DestroyMoleculeRouter.delete('/:id', (req, res) => {
  (async () => {
    const id = req.params.id;

    if (!id || typeof id !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const user = req.full_user!;

    if (!user) {
      return Errors.throw(ErrorType.Forbidden);
    }

    const delResp = await deleteMolecule(id, user, false, true);

    res.send(delResp);
  })().catch(errorCatcher(res));
});

DestroyMoleculeRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyMoleculeRouter;
