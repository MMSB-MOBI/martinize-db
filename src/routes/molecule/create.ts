import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, sanitize, informAdminFromNewMolecule, dumpStdFromDir } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType, ApiError } from '../../Errors';
import { Molecule, StashedMolecule, BaseMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import { DISABLE_MODERATION_PROCESS } from '../../constants';
import { MoleculeChecker } from './MoleculeChecker';
import logger from '../../logger';
import { plainToInstance } from 'class-transformer'; 
import { validateOrReject } from 'class-validator'; 
import { FileDto } from './membrane_builder.dto'; 



//import * as test from "../../test";


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
  { name: 'map', maxCount: 99 },
]), (req, res) => {
  (async () => {  //console.log("HERE>>");console.log(req.files);

    //Validate file names and length
    const files = req.files as { [fieldname: string]: Express.Multer.File[] };
    if(files.pdb.length > 1) return Errors.throw(ErrorType.TooManyFiles)
    if(files.top.length > 1) return Errors.throw(ErrorType.TooManyFiles)
    const validatedPdb = plainToInstance(FileDto, files.pdb[0])
    const validatedTop = plainToInstance(FileDto, files.top[0])
    const validatedItps = files.itp.map(itp => plainToInstance(FileDto, itp))
    const validatedMaps = files.map ? files.map.map(map => plainToInstance(FileDto, map)) : undefined

    try {
      await validateOrReject(validatedPdb); 
      await validateOrReject(validatedTop); 
      await Promise.all(validatedItps.map(dto => validateOrReject(dto)))
      if(validatedMaps) await Promise.all(validatedMaps.map(dto => validateOrReject(dto)))
    } catch(e) {
      res.status(400).json({ error: true, statusCode: 400, errorCode: 'PARAMS_VALIDATION_ERROR', e })
      return; 
    }
    const logged_user = req.full_user!;
    
    const user_role = DISABLE_MODERATION_PROCESS ? "admin" : logged_user.role;
    const force_moderation = req.query.force_moderation === "1";

    // Insert the molecule
    let response: nano.DocumentInsertResponse;
    let molecule: BaseMolecule;

    req.body.category = typeof req.body.category === 'string' ? [req.body.category] : req.body.category

    try {
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
        informAdminFromNewMolecule(molecule as StashedMolecule, logged_user).catch(logger.error);
      }
  
      if (response.ok) {
        res.json(sanitize(molecule));
      }
      else {
        return Errors.throw(ErrorType.Server);
      }
    } catch (e) {
      if (
        e instanceof ApiError && 
        e.code === ErrorType.InvalidMoleculeFiles && 
        e.data && 
        e.data.dir
      ) {
        const { stdout, stderr } = await dumpStdFromDir(e.data.dir as string);

        return Errors.throw(ErrorType.InvalidMoleculeFiles, { stdout, stderr });
      }

      throw e;
    }
  })().catch(errorCatcher(res));
});

CreateMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default CreateMoleculeRouter;
