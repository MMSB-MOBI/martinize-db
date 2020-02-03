import { Router } from 'express';
import Errors, { ErrorType } from '../../../Errors';

const CreateUserRouter = Router();

CreateUserRouter.post('/', (req, res) => {
  if (!req.body) {
    Errors.send(ErrorType.MissingParameters, res);
    return;
  } 

  res.json({});
});

export default CreateUserRouter;
