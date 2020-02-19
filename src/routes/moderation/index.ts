import { Router } from 'express';
import GetStashedRouter from './get';
import EditStashedRouter from './edit';
import DestroyStashedRouter from './destroy';
import ListModerationRouter from './list';
import AcceptModerationRouter from './accept';

const ModerationRouter = Router();

ModerationRouter.use('/list', ListModerationRouter);
ModerationRouter.use('/destroy', DestroyStashedRouter);
ModerationRouter.use('/edit', EditStashedRouter);
ModerationRouter.use('/', GetStashedRouter);
ModerationRouter.use('/accept', AcceptModerationRouter);

export default ModerationRouter;
