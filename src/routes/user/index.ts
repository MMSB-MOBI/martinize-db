import { Router } from 'express';
import LoginUserRouter from './login';
import CreateUserRouter from './create';

const UserRouter = Router();

UserRouter.use('/login', LoginUserRouter);
UserRouter.use('/create', CreateUserRouter);

export default UserRouter;
