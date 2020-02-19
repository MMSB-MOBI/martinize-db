import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import SearchWorker from '../../search_worker';

const DestroyMoleculeRouter = Router();

DestroyMoleculeRouter.delete('/:id', (req, res) => {
  (async () => {
    // TODO MAKE SEARCH
    // For now, it only returns every molecule
    const id = req.params.id;

    if (!id || typeof id !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const user = req.full_user!;

    if (!user) {
      return Errors.throw(ErrorType.Forbidden);
    }

    try {
      const mol = await Database.molecule.get(id);

      if (user.role !== "admin" && user.id !== mol.owner) {
        return Errors.throw(ErrorType.Forbidden);
      }

      // Recherche les versions attachées à cette molecule
      const versions_attached = await Database.molecule.find({
        limit: 99999,
        selector: {
          parent: mol.id,
          tree_id: mol.tree_id,
        },
      });

      // Met à jour les liens de parenté
      if (versions_attached.length) {
        if (mol.parent) {
          // Every molecule will be attached to the new parent
          for (const m of versions_attached) {
            m.parent = mol.parent;
          }
        }
        else {
          // The first molecule will be the new parent for everyone
          const first = versions_attached[0];
          for (const other of versions_attached.slice(1)) {
            other.parent = first.id;
          }
          first.parent = null;
        }

        // Save everyone
        await Promise.all(versions_attached.map(v => Database.molecule.save(v)));
        SearchWorker.clearCache();
      }

      // Delete attached ZIP
      await MoleculeOrganizer.remove(mol.files);
      await Database.molecule.delete(mol);
    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    res.send();
  })().catch(errorCatcher(res));
});

DestroyMoleculeRouter.all('/', methodNotAllowed(['DELETE']))

export default DestroyMoleculeRouter;
