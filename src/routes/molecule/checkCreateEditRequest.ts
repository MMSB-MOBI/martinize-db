import { Request } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { MAX_FILE_SIZE } from '../Uploader';
import { BaseMolecule, Molecule } from '../../Entities/entities';
import { generateSnowflake, isDebugMode, verifyAndCompleteMolecule } from '../../helpers';
import MoleculeOrganizer, { MoleculeSave } from '../../MoleculeOrganizer';

export default async function checkCreateOrEditRequest(req: Request, edit = false, stashed = false) {
  if (!req.files || Array.isArray(req.files)) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  let { itp, gro, pdb } = req.files;

  // If gro AND pdb
  if (gro && gro.length && pdb && pdb.length) {
    return Errors.throw(ErrorType.TooManyFiles);
  }

  if (!edit && !itp.length) {
    // If no PDB file
    return Errors.throw(ErrorType.MissingParameters);
  }  

  if (!pdb) {
    pdb = [];
  }

  // Take only one PDB file
  let pdb_file = pdb[0];

  // If there is one GRO file, take it
  if (gro && gro.length) {
    pdb_file = gro[0];

    // TODO convert gro to pdb ! (?)
  }

  if (!pdb_file && !edit) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  if (pdb_file) {
    if (pdb_file.size > MAX_FILE_SIZE) {
      // > 1 GB
      return Errors.throw(ErrorType.FileTooLarge, { excepted_size: MAX_FILE_SIZE, actual_size: pdb_file.size });
    }
  }

  if (!req.body) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  const logged_user = req.full_user!;

  if (!logged_user) {
    return Errors.throw(ErrorType.Forbidden);
  }

  const b = req.body;
  const id = edit ? b.id : generateSnowflake();
  // Create the molecule
  let molecule: BaseMolecule = {
    id,
    name: b.name,
    alias: b.alias,
    formula: b.formula,
    version: b.version,
    category: b.category,
    parent: b.parent,
    tree_id: b.tree_id ?? "",
    hash: b.hash ?? "",
    owner: logged_user.id,
    files: "",
    comments: b.comments,
    created_at: new Date().toISOString(),
    command_line: b.command_line,
    martinize_version: b.martinize_version,
    force_field: b.force_field,
  };

  if (b.parent === "null") {
    molecule.parent = null;
  }

  if (edit && b.files && typeof b.files === 'string') {
    molecule.files = b.files;
  }

  if (edit) {
    // Make some check todo more
    if (!molecule.id || !molecule.files || !molecule.hash) {
      return Errors.throw(ErrorType.MissingParameters);
    }
  }

  // Check the user role
  const user_role = isDebugMode() ? "admin" : logged_user.role;

  // If logged user is admin, set the required parameters
  if (user_role === "admin" && !stashed) {
    (molecule as Molecule).last_update = new Date().toISOString();
    (molecule as Molecule).approved_by = logged_user.id;
  }

  try {
    // this will define tree_id properly, and check every property
    molecule = await verifyAndCompleteMolecule(molecule, edit, stashed);
  } catch (e) {
    return Errors.throw(ErrorType.Format, e);
  }

  // If file changed
  if (pdb_file && itp.length) {
    // Save the molecule in ZIP format if files changed
    let save: MoleculeSave;
    try {
      save = await MoleculeOrganizer.save(itp, pdb_file);
    } catch (e) {
      return Errors.throw(ErrorType.InvalidMoleculeFiles, e);
    }
    molecule.hash = save.infos.hash;
    molecule.files = save.id;
  }

  return molecule;
}
