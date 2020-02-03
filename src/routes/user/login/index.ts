import { Router } from 'express';
import { errorCatcher, generateSnowflake, signToken } from '../../../helpers';
import Errors, { ErrorType } from '../../../Errors';
import { Database } from '../../../Entities/CouchHelper';
import { Token } from '../../../Entities/entities';

const LoginUserRouter = Router();

LoginUserRouter.post('/', (req, res) => {
  (async () => {
    if (!req.body) {
      Errors.throw(ErrorType.MissingParameters);
    } 

    const { username, password } = req.body as { username: string, password: string };

    if (!username || !password) {
      Errors.throw(ErrorType.MissingParameters);
    }

    const user = await Database.user.fromUsername(username);

    if (!user) {
      return Errors.throw(ErrorType.UserNotFound);
    }

    const is_connected = await Database.user.verifyPassword(user, password);

    if (!is_connected) {
      Errors.throw(ErrorType.InvalidPassword);
    }

    // Create a token for this user
    const token: Token = {
      id: generateSnowflake(),
      user_id: user.id,
      created_at: new Date().toISOString(),
    }; 
    
    await Database.token.save(token);

    // Create the encoded JSON Web Token
    const jwt = await signToken({ user_id: token.user_id, created_at: token.created_at }, token.id);

    res.json({
      token: jwt
    });
  })().catch(errorCatcher(res));
});

export default LoginUserRouter;
