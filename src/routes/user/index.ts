import { Router } from 'express';
import LoginUserRouter from './login';
import CreateUserRouter from './create';
import RevokeTokenRouter from './revoke';
import ValidateUserRouter from './validate';
import UpdateUserRouter from './update';
import ListUserRouter from './list';

const UserRouter = Router();

UserRouter.use('/login', LoginUserRouter);
UserRouter.use('/create', CreateUserRouter);
UserRouter.use('/revoke', RevokeTokenRouter);
UserRouter.use('/validate', ValidateUserRouter);
UserRouter.use('/update', UpdateUserRouter);
UserRouter.use('/list', ListUserRouter);

export default UserRouter;
