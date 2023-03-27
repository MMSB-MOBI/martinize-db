import { Router, Request } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, sanitize, dumpStdFromDir } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType, ApiError } from '../../Errors';
import { StashedMolecule, User } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import { MoleculeChecker } from '../molecule/MoleculeChecker';
import { plainToInstance } from 'class-transformer'; 
import { validateOrReject } from 'class-validator'; 
import { FileDto } from '../molecule/membrane_builder.dto'; 

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
  { name: 'top', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
  { name: 'map', maxCount: 99 },
]), (req: Request, res) => {
  // Saving the file
  // @ts-ignore

  req.body.category = typeof req.body.category === 'string' ? [req.body.category] : req.body.category

  //Validate file names and length
    //TO DO : factorize this, it's also used in create.ts
    // const files = req.files as { [fieldname: string]: Express.Multer.File[] };
    // console.log(files); 
    // if(files.pdb.length > 1) return Errors.throw(ErrorType.TooManyFiles)
    // if(files.top.length > 1) return Errors.throw(ErrorType.TooManyFiles)
    // const validatedPdb = plainToInstance(FileDto, files.pdb[0])
    // const validatedTop = plainToInstance(FileDto, files.top[0])
    // const validatedItps = files.itp.map(itp => plainToInstance(FileDto, itp))
    // const validatedMaps = files.map ? files.map.map(map => plainToInstance(FileDto, map)) : undefined

    

  const logged_user = req.full_user! as User;
  (async () => {
    // try {
    //   await validateOrReject(validatedPdb); 
    //   await validateOrReject(validatedTop); 
    //   await Promise.all(validatedItps.map(dto => validateOrReject(dto)))
    //   if(validatedMaps) await Promise.all(validatedMaps.map(dto => validateOrReject(dto)))
    // } catch(e) {
    //   res.status(400).json({ error: true, statusCode: 400, errorCode: 'PARAMS_VALIDATION_ERROR', e })
    //   return; 
    // }

    if (logged_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const checker = new MoleculeChecker(req);

    try {
      const molecule = await checker.checkStashedEdition();
  
      let response = await Database.stashed.save(molecule as StashedMolecule);
  
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

EditStashedRouter.all('/', methodNotAllowed(['POST']));

export default EditStashedRouter;
