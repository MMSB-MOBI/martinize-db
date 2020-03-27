import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, sanitize, informAdminFromNewMolecule } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import { Molecule, StashedMolecule, BaseMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import { DISABLE_MODERATION_PROCESS } from '../../constants';
import { MoleculeChecker } from './MoleculeChecker';
import logger from '../../logger';

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
  { name: 'top', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req, res) => {
  (async () => {
    const logged_user = req.full_user!;
    
    const user_role = DISABLE_MODERATION_PROCESS ? "admin" : logged_user.role;
    const force_moderation = req.query.force_moderation === "1";

    // Insert the molecule
    let response: nano.DocumentInsertResponse;
    let molecule: BaseMolecule;

    const checker = new MoleculeChecker(req);

    if (user_role === "admin" && !force_moderation) {
      // Inset directly in molecule db
      molecule = await checker.check();
      response = await Database.molecule.save(molecule as Molecule);
    }
    else {
      molecule = await checker.checkStashed();
      response = await Database.stashed.save(molecule as StashedMolecule);

      // Inform moderators
      informAdminFromNewMolecule(molecule as StashedMolecule, logged_user).catch(logger.error)
    }

    if (response.ok) {
      res.json(sanitize(molecule));
    }
    else {
      return Errors.throw(ErrorType.Server);
    }
  })().catch(errorCatcher(res));
});

CreateMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default CreateMoleculeRouter;
