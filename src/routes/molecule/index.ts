import { Router } from 'express';
import CreateMoleculeRouter from './create';
import ListMoleculeRouter from './list';
import DestroyMoleculeRouter from './destroy';
import GetMoleculeRouter from './get';
import DownloadMoleculeRouter from './download';
import EditMoleculeRouter from './edit';
import PdbGetterRouter from './pdb';
import MartinizerRouter from './martinize';
import MembraneBuilderRouter from './membrane_builder';
import GetMoleculeAPI from './getMoleculeAPI';


const MoleculeRouter = Router();

MoleculeRouter.use('/create', CreateMoleculeRouter);
MoleculeRouter.use('/list', ListMoleculeRouter);
MoleculeRouter.use('/destroy', DestroyMoleculeRouter);
MoleculeRouter.use('/edit', EditMoleculeRouter);
MoleculeRouter.use('/download', DownloadMoleculeRouter);
MoleculeRouter.use('/representation', PdbGetterRouter);
MoleculeRouter.use('/martinize', MartinizerRouter);
MoleculeRouter.use('/membrane_builder', MembraneBuilderRouter);
MoleculeRouter.use('/', GetMoleculeRouter);
MoleculeRouter.use('/get', GetMoleculeAPI);

export default MoleculeRouter;
