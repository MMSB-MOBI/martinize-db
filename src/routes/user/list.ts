import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import logger from '../../logger';

const ListUserRouter = Router();

ListUserRouter.get('/waiting', (req, res) => {  
  (async () => {
    if (!req.full_user || req.full_user.role !== "admin") {
      return Errors.throw(ErrorType.Forbidden);
    }

    const all = await Database.user.all();

    const waiting = all.filter(usr => usr.approved === false)

    res.json(waiting.map(u => sanitize(u)));
    
  })().catch(errorCatcher(res));
});

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
