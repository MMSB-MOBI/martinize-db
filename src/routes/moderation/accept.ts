import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import { StashedMolecule, Molecule } from '../../Entities/entities';
import logger from '../../logger';
import Mailer from '../../Mailer/Mailer';
import { URLS } from '../../constants';

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
    const molecule = sanitize({ ...original }) as StashedMolecule;

    const full = molecule as Molecule;
    full.last_update = full.created_at;
    full.approved_by = user.id;

    // Add in full
    const response = await Database.molecule.save(full);

    if (!response.ok) {
      return Errors.throw(ErrorType.Server);
    }

    // Remove from stashed
    try {
      await Database.stashed.delete(original);
    } catch (e) {
      logger.error("Unable to delete stashed mol", e);
    }

    // Send mail to molecule sender
    const sender = await Database.user.get(full.owner);
    Mailer.send({
      to: sender.email,
      subject: "MArtini Database - Molecule approved",
    }, 'mail_molecule_accepted', {
      name: sender.name,
      molecule: full,
      molecule_url: URLS.SERVER + '/molecule/' + full.alias + '?version=' + full.id,
    }).catch(logger.error);

    res.json(sanitize(full));
  })().catch(errorCatcher(res));
});

AcceptModerationRouter.all('/', methodNotAllowed('POST'));

export default AcceptModerationRouter;
