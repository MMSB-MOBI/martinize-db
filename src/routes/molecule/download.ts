import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';

const DownloadMoleculeRouter = Router();

DownloadMoleculeRouter.get('/', (req, res) => {
  (async () => {
    const { id, filename } = req.query;
    
    if (!id) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    try {
      BigInt(id);
    } catch {
      return Errors.throw(ErrorType.Format);
    }

    if (!(await MoleculeOrganizer.exists(id))) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    // todo check filename

    res.download(MoleculeOrganizer.getFilenameFor(id), filename || undefined);
  })().catch(errorCatcher(res));
});

DownloadMoleculeRouter.all('/', methodNotAllowed('GET'))

export default DownloadMoleculeRouter;
