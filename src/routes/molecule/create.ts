import { Router } from 'express';
import { methodNotAllowed, errorCatcher, generateSnowflake, verifyAndCompleteMolecule } from '../../helpers';
import Uploader, { MAX_FILE_SIZE } from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import { Molecule, BaseMolecule, StashedMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');

const CreateMoleculeRouter = Router();

/**
 * Create a new molecule.
 * 
 * Expect a enctype `multipart/form-data`, because you should have two files in it.
 * 
 * Expect file `itp=itp_file`, `pdb=pdb_file`, `gro=gro_file`.
 * You must specify an ITP, and at least one PDB or one GRO file.
 * 
 */
CreateMoleculeRouter.post('/', Uploader.fields([
  { name: 'itp', maxCount: 1 }, 
  { name: 'gro', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req, res) => {
  if (!req.files || Array.isArray(req.files)) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  const { itp, gro, pdb } = req.files;

  if (gro.length && pdb.length) {
    return Errors.throw(ErrorType.TooManyFiles);
  }

  if (!itp.length) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  let pdb_file = pdb[0];

  if (gro.length) {
    pdb_file = gro[0];

    // TODO convert gro to pdb ! (?)
  }

  if (pdb_file.size > MAX_FILE_SIZE) {
    // > 1 GB
    return Errors.throw(ErrorType.FileTooLarge, { excepted_size: MAX_FILE_SIZE, actual_size: pdb_file.size });
  }

  if (!req.body) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  // Saving the file
  (async () => {
    const logged_user = req.full_user!;

    if (!logged_user) {
      return Errors.throw(ErrorType.Forbidden);
    }

    const save = await MoleculeOrganizer.save(itp[0].filename, pdb_file.filename);

    
    const b = req.body;
    // Create the molecule
    let molecule: BaseMolecule = {
      id: generateSnowflake(),
      name: b.name,
      alias: b.alias,
      formula: b.formula,
      version: b.version,
      category: b.category,
      parent: b.parent,
      tree_id: "",
      hash: save.infos.hash,
      owner: logged_user.id,
      files: save.id,
      comments: b.comments,
      created_at: new Date().toISOString(),
      command_line: b.command_line,
      martinize_version: b.martinize_version,
      force_field: b.force_field,
    };

    if (logged_user.role === "admin") {
      (molecule as Molecule).last_update = new Date().toISOString();
      (molecule as Molecule).approved_by = logged_user.id;
    }

    try {
      // this will define tree_id properly, and check every property
      molecule = await verifyAndCompleteMolecule(molecule);
    } catch (e) {
      return Errors.throw(ErrorType.Format, e);
    }

    let response: nano.DocumentInsertResponse;
    // Insert the molecule NOT STASHED //// TODO DEBUG REMOVE ////
    if (logged_user.role === "admin" || true) {
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
