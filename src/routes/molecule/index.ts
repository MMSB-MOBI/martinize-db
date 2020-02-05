import { Router } from 'express';
import CreateMoleculeRouter from './create';
import ListMoleculeRouter from './list';
import DestroyMoleculeRouter from './destroy';

const MoleculeRouter = Router();

MoleculeRouter.use('/create', CreateMoleculeRouter);
MoleculeRouter.use('/list', ListMoleculeRouter);
MoleculeRouter.use('/destroy', DestroyMoleculeRouter);

export default MoleculeRouter;
