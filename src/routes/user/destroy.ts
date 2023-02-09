import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, methodNotAllowed, deleteMolecule } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Mailer from '../../Mailer/Mailer';

const DestroyUserRouter = Router();

DestroyUserRouter.delete('/', (req, res) => {  
  (async () => {
    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const id = req.query.id as string;

    if (!id) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    // Check if user exists
    const to_delete_user = await Database.user.exists(id);
    if (!to_delete_user) {
      return Errors.throw(ErrorType.UserNotFound);
    }
    const user = await Database.user.get(id);

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

    await Database.user.delete(user);

    if(! user.approved) {
      await Mailer.send({ 
        to: user.email, 
        subject: "MArtini Database - " + user.name + ": Your account has been rejected" 
      }, 'mail_rejected', { 
        title: user.name + ": Your account has been rejected",
        new_user: {
          name: user.name,
        },
      });
    }

    res.json({
      deleted: true
    });
  })().catch(errorCatcher(res));
});

DestroyUserRouter.all('/', methodNotAllowed('DELETE'));

export default DestroyUserRouter;
