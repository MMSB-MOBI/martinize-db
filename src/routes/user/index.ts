import { Router } from 'express';
import LoginUserRouter from './login';
import CreateUserRouter from './create';
import RevokeTokenRouter from './revoke';
import ValidateUserRouter from './validate';

const UserRouter = Router();

UserRouter.use('/login', LoginUserRouter);
UserRouter.use('/create', CreateUserRouter);
UserRouter.use('/revoke', RevokeTokenRouter);
UserRouter.use('/validate', ValidateUserRouter);

export default UserRouter;
