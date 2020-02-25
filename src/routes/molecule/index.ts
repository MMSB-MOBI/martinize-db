import { Router } from 'express';
import CreateMoleculeRouter from './create';
import ListMoleculeRouter from './list';
import DestroyMoleculeRouter from './destroy';
import GetMoleculeRouter from './get';
import DownloadMoleculeRouter from './download';
import EditMoleculeRouter from './edit';
import PdbGetterRouter from './pdb';

const MoleculeRouter = Router();

MoleculeRouter.use('/create', CreateMoleculeRouter);
MoleculeRouter.use('/list', ListMoleculeRouter);
MoleculeRouter.use('/destroy', DestroyMoleculeRouter);
MoleculeRouter.use('/edit', EditMoleculeRouter);
MoleculeRouter.use('/download', DownloadMoleculeRouter);
MoleculeRouter.use('/representation', PdbGetterRouter);
MoleculeRouter.use('/', GetMoleculeRouter);

export default MoleculeRouter;
