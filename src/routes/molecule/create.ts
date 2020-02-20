import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import { Molecule, StashedMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import checkCreateOrEditRequest from './checkCreateEditRequest';
import { DISABLE_MODERATION_PROCESS } from '../../constants';

const CreateMoleculeRouter = Router();

// Middleware that wipe uploaded files after request
CreateMoleculeRouter.use((req, res, next) => {
  function after() {
    // Response is sended
    cleanMulterFiles(req);
    res.removeListener('finish', after);
  }

  res.once('finish', after);
  next();
});

/**
 * Create a new molecule.
 * 
 * Expect a enctype `multipart/form-data`, because you should have two files in it.
 * 
 * Expect file `itp=itp_file`, `pdb=pdb_file`, `gro=gro_file`.
 * You must specify an ITP, and at least one PDB or one GRO file.
 * 
 * 
 */
CreateMoleculeRouter.post('/', Uploader.fields([
  { name: 'itp', maxCount: 99 }, 
  { name: 'gro', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req, res) => {
  (async () => {
    const molecule = await checkCreateOrEditRequest(req);

    const logged_user = req.full_user!;
    
    const user_role = DISABLE_MODERATION_PROCESS ? "admin" : logged_user.role;

    // Insert the molecule NOT STASHED //// TODO DEBUG REMOVE ////
    let response: nano.DocumentInsertResponse;
    if (user_role === "admin") {
      response = await Database.molecule.save(molecule as Molecule);
    }
    else {
      response = await Database.stashed.save(molecule as StashedMolecule);
    }

    if (response.ok) {
      res.json(molecule);
    }
    else {
      return Errors.throw(ErrorType.Server);
    }
  })().catch(errorCatcher(res));
});

CreateMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default CreateMoleculeRouter;
