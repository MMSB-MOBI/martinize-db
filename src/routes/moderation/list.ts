import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import Errors, { ErrorType } from '../../Errors';

const ListModerationRouter = Router();

ListModerationRouter.get('/', (req, res) => {
  (async () => {
    const user = req.full_user!;

    if (!user || user.role !== 'admin') {
      return Errors.throw(ErrorType.Forbidden);
    }

    const query: nano.MangoQuery = { selector: {}, limit: 25, skip: 0 };

    const { 
      skip,
      limit,
    } = req.query; 

    if (skip) {
      if (Number(skip) > 0) {
        query.skip = Number(skip);
      }
    }
    if (limit) {
      const l = Number(limit);

      if (l > 0 && l <= 200) {
        query.limit = l;
      }
    }

    // Search in stashed
    const molecules = sanitize(await Database.stashed.find(query));

    res.json({
      molecules,
      length: await Database.stashed.count()
    });
  })().catch(errorCatcher(res));
});

ListModerationRouter.all('/', methodNotAllowed('GET'));

export default ListModerationRouter;
