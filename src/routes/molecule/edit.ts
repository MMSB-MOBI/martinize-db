import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, sanitize, dumpStdFromDir } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType, ApiError } from '../../Errors';
import { Molecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import { MoleculeChecker } from './MoleculeChecker';
import { plainToInstance } from 'class-transformer'; 
import { validateOrReject } from 'class-validator'; 
import { FileDto } from './membrane_builder.dto'; 

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
  { name: 'top', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
  { name: 'map', maxCount: 99 },
]), (req, res) => {
  // Saving the file
  (async () => {

    //Validate file names and length
    //TO DO : factorize this, it's also used in create.ts
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

    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    try {
      const checker = new MoleculeChecker(req);
      const molecule = await checker.checkEdition();

      if (molecule.parent === null) {
        // Should refresh all the possible childs with name, alias, formula, category

        const children = (await Database.molecule.find({ limit: 99999, selector: { tree_id: molecule.tree_id } })).filter(m => m.id !== molecule.id);

        const saves: Promise<any>[] = [];

        for (const child of children) {
          child.name = molecule.name;
          child.alias = molecule.alias;
          child.smiles = molecule.smiles;
          child.category = molecule.category;

          saves.push(Database.molecule.save(child));
        }

        await Promise.all(saves);
      }

      // Save the edited molecule
      const response = await Database.molecule.save(molecule as Molecule);

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

EditMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default EditMoleculeRouter;
