import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { Token } from '../../Entities/entities';

const RevokeTokenRouter = Router();

RevokeTokenRouter.delete('/', (req, res) => {
  (async () => {
    let token_to_revoke = req.user!.jti;

    if (req.query && req.query.token) {
      token_to_revoke = req.query.token as string;
    }

    if (!token_to_revoke) {
      Errors.throw(ErrorType.MissingParameters);
    }

    let token: Token;
    try  {
      token = await Database.token.get(token_to_revoke);
    } catch (e) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    if (token.user_id !== req.user!.user_id) {
      return Errors.throw(ErrorType.Forbidden);
    }

    // Token exists
    await Database.token.delete(token);

    res.send();
  })().catch(errorCatcher(res));
});

RevokeTokenRouter.all('/', methodNotAllowed('DELETE'));

export default RevokeTokenRouter;
