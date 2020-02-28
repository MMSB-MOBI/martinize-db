import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';

const ListUserRouter = Router();

ListUserRouter.get('/', (req, res) => {  
  (async () => {
    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const all = await Database.user.all();

    res.json(all.map(u => sanitize(u)));
  })().catch(errorCatcher(res));
});

ListUserRouter.all('/', methodNotAllowed('GET'));

export default ListUserRouter;
