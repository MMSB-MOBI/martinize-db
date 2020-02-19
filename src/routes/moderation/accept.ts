import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import { StashedMolecule, Molecule } from '../../Entities/entities';

const AcceptModerationRouter = Router();

AcceptModerationRouter.post('/', (req, res) => {
  (async () => {
    const user = req.full_user!;

    if (!user || user.role !== 'admin') {
      return Errors.throw(ErrorType.Forbidden);
    }

    const id = req.body.id;
    if (!id) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    if (!await Database.stashed.exists(id)) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const original = await Database.stashed.get(id);
    const molecule = sanitize(original) as StashedMolecule;

    const full = molecule as Molecule;
    full.last_update = new Date().toISOString();
    full.approved_by = user.id;

    // Remove from stashed
    await Database.stashed.delete(original);
    // Add in full
    const response = await Database.molecule.save(full);

    if (!response.ok) {
      return Errors.throw(ErrorType.Server);
    }

    res.json(sanitize(full));
  })().catch(errorCatcher(res));
});

AcceptModerationRouter.all('/', methodNotAllowed('POST'));

export default AcceptModerationRouter;
