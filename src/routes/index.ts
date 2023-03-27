import { Router } from 'express';
import cors from 'cors';
import bodyParser from 'body-parser';
import UserRouter from './user';
import cookieParser from 'cookie-parser';
import jwt from './jwt';
import Errors, { ErrorType } from '../Errors';
import MoleculeRouter from './molecule';
import SettingsRouter from './settings';
import ModerationRouter from './moderation';
import DownloadForceFieldRoute from './force_fields/download';
import HistoryRouter from './history'
import polymer from './polymergenerator';

const ApiRouter = Router();
ApiRouter.use(cors());
ApiRouter.use(bodyParser.urlencoded({ extended: true }));
ApiRouter.use(bodyParser.json());
ApiRouter.use(cookieParser());

ApiRouter.use(jwt);

// Subscribe to sub routers
ApiRouter.use('/user', UserRouter);
ApiRouter.use('/molecule', MoleculeRouter);
ApiRouter.use('/settings', SettingsRouter);
ApiRouter.use('/moderation', ModerationRouter);
ApiRouter.use('/force_fields', DownloadForceFieldRoute);
ApiRouter.use('/history', HistoryRouter)
ApiRouter.use('/polymergenerator', polymer)

// Catch all API invalid routes
ApiRouter.use(() => {
  Errors.throw(ErrorType.NotFound);
});

export default ApiRouter;
