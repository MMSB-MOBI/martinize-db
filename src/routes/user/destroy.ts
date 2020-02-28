import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, methodNotAllowed, deleteMolecule } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';

const DestroyUserRouter = Router();

DestroyUserRouter.get('/', (req, res) => {  
  (async () => {
    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const id = req.query.id;

    if (!id) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    // Check if user exists
    const to_delete_user = await Database.user.exists(id);
    if (!to_delete_user) {
      return Errors.throw(ErrorType.UserNotFound);
    }

    // Get and Delete all stashed molecules + published molecules
    const stashed_for_user = await Database.stashed.find({
      selector: {
        owner: id
      },
      limit: 99999
    });

    const molecule_for_user = await Database.molecule.find({
      selector: {
        owner: id
      },
      limit: 99999
    });

    await Promise.all(stashed_for_user.map(s => deleteMolecule(s.id, req.full_user!, true)));
    await Promise.all(molecule_for_user.map(s => deleteMolecule(s.id, req.full_user!, false, false)));

    res.json({
      deleted: true
    });
  })().catch(errorCatcher(res));
});

DestroyUserRouter.all('/', methodNotAllowed('DELETE'));

export default DestroyUserRouter;
