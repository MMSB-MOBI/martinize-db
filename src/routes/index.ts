import { Router } from 'express';
import cors from 'cors';
import bodyParser from 'body-parser';
import UserRouter from './user';
import cookieParser from 'cookie-parser';
import jwt from './jwt';

const ApiRouter = ();
ApiRouter.use(cors());
ApiRouter.use(bodyParser.urlencoded());
ApiRouter.use(bodyParser.json());
ApiRouter.use(cookieParser());

ApiRouter.use(jwt);

ApiRouter.use('/user', UserRouter);

export default ApiRouter;
