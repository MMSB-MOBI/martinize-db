import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles, sanitize } from '../../helpers';
import Uploader from '../Uploader';
import Errors, { ErrorType } from '../../Errors';
import { Molecule, StashedMolecule } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import { MoleculeChecker } from './MoleculeChecker';

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
    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const checker = new MoleculeChecker(req);
    const molecule = await checker.checkEdition();

    if (molecule.parent === null) {
      // Should refresh all the possible childs with name, alias, formula, category

      const children = (await Database.molecule.find({ limit: 99999, selector: { tree_id: molecule.tree_id } })).filter(m => m.id !== molecule.id);

      const saves: Promise<any>[] = [];

      for (const child of children) {
        child.name = molecule.name;
        child.alias = molecule.alias;
        child.formula = molecule.formula;
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
  })().catch(errorCatcher(res));
});

EditMoleculeRouter.all('/', methodNotAllowed(['POST']));

export default EditMoleculeRouter;
