import { Router, Request } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import { StashedMolecule, User } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import checkCreateOrEditRequest from '../molecule/checkCreateEditRequest';

const EditStashedRouter = Router();

// Middleware that wipe uploaded files after request
EditStashedRouter.use((req, res, next) => {
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
EditStashedRouter.post('/', Uploader.fields([
  { name: 'itp', maxCount: 99 }, 
  { name: 'gro', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req: Request, res) => {
  // Saving the file
  // @ts-ignore
  const logged_user = req.full_user! as User;
  (async () => {
    if (logged_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const molecule = await checkCreateOrEditRequest(req, true);

    let response = await Database.stashed.save(molecule as StashedMolecule);

    if (response.ok) {
      res.json(molecule);
    }
    else {
      return Errors.throw(ErrorType.Server);
    }
  })().catch(errorCatcher(res));
});

EditStashedRouter.all('/', methodNotAllowed(['POST']));

export default EditStashedRouter;
