import { Router } from 'express';
import { errorCatcher, methodNotAllowed, deleteMolecule } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import Mailer from '../../Mailer/Mailer';

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

    const molecule = await Database.stashed.get(id);
    const owner = await Database.user.get(molecule.owner);
    await deleteMolecule(id, user, true);

    // Inform user
    await Mailer.send(
      { 
        to: owner.email, 
        subject: "Molecule rejected - MArtini Database" 
      }, 
      'mail_molecule_rejected',
      {
        name: owner.name,
        molecule
      }  
    );

    res.send();
  })().catch(errorCatcher(res));
});

DestroyStashedRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyStashedRouter;
