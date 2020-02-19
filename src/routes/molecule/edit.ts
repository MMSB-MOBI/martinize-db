import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, isDebugMode } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import { Molecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import checkCreateOrEditRequest from './checkCreateEditRequest';

const EditMoleculeRouter = Router();

// Middleware that wipe uploaded files after request
EditMoleculeRouter.use((req, res, next) => {
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
EditMoleculeRouter.post('/', Uploader.fields([
  { name: 'itp', maxCount: 99 }, 
  { name: 'gro', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req, res) => {
  // Saving the file
  (async () => {
    const molecule = await checkCreateOrEditRequest(req, true);

    const logged_user = req.full_user!;
    const user_role = isDebugMode() ? "admin" : logged_user.role;

    // Insert the molecule NOT STASHED //// TODO DEBUG REMOVE ////
    let response = await Database.molecule.save(molecule as Molecule);

    if (response.ok) {
      res.json(molecule);
    }
    else {
      return Errors.throw(ErrorType.Server);
    }
  })().catch(errorCatcher(res));
});

EditMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default EditMoleculeRouter;
