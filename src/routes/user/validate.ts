import { Router } from 'express';
import { errorCatcher, methodNotAllowed, sanitize } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';

const ValidateUserRouter = Router();

ValidateUserRouter.get('/', (req, res) => {
  (async () => {
    const user = await Database.user.fromToken(req.user!.jti);

    if (!user) {
      return Errors.throw(ErrorType.Forbidden);
    }
    if (!user.approved) {
      return Errors.throw(ErrorType.UserNotApproved);
    }

    res.json(sanitize({ ...user, password: null }));
  })().catch(errorCatcher(res));
});

ValidateUserRouter.all('/', methodNotAllowed('GET'));

export default ValidateUserRouter;
