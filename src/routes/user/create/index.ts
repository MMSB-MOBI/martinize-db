import { Router } from 'express';

const CreateUserRouter = Router();

CreateUserRouter.post('/', (req, res) => {
  res.json({});
});

export default CreateUserRouter;
