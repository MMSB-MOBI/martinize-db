import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, generateSnowflake, signToken, sanitize, methodNotAllowed } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import { UserRole } from '../../types';
import { User, Token } from '../../Entities/entities';

const CreateUserRouter = Router();

CreateUserRouter.post('/', (req, res) => {
  if (!req.body) {
    Errors.throw(ErrorType.MissingParameters);
  } 

  // TODO CHECK VALIDITY OF ALL DATA
  let { username, password, email, role: choosen_role } = req.body as { username?: string, password?: string, email?: string, role?: string }; 
  
  (async () => {
    let role: UserRole = "admin"; 

    if (!username || !password || !email) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    if (req.user) {
      const user = await Database.user.fromToken(req.user.jti);

      // If user is admin, allow all roles to be defined.
      if (user && user.role === 'admin') {
        role = choosen_role as UserRole;
      }
    }

    // Check if username and email already exists
    const user1 = await Database.user.fromUsername(username!);
    if (user1) {
      Errors.throw(ErrorType.UsernameExists);
    }

    const user2 = await Database.user.fromEmail(email!);
    if (user2) {
      Errors.throw(ErrorType.EmailExists);
    }

    const user: User = {
      id: generateSnowflake(),
      email,
      name: username,
      role,
      created_at: new Date().toISOString(),
      password: "",
    };

    // Register the user
    await Database.user.setOrChangePassword(user, password);

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
      created: sanitize({ ...user, password: null }),
      token: jwt
    });
  })().catch(errorCatcher(res));
});

CreateUserRouter.all('/', methodNotAllowed('POST'));

export default CreateUserRouter;
