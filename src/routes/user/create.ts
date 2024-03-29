import { Router } from 'express';
import Errors, { ErrorType } from '../../Errors';
import { errorCatcher, generateSnowflake, signToken, sanitize, methodNotAllowed, informAdminFromAskCreation } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import { UserRole } from '../../types';
import { User, Token } from '../../Entities/entities';
import { USERNAME_REGEX, EMAIL_REGEX } from '../../constants';
import logger from '../../logger';
import util from 'util'; 

const CreateUserRouter = Router();

CreateUserRouter.post('/', (req, res) => {
  if (!req.body) {
    Errors.throw(ErrorType.MissingParameters);
  } 

  let { username, password, email, fullname, affiliation, role: choosen_role } = req.body as { username?: string, password?: string, email?: string, role?: string, fullname?: string, affiliation?:string }; 
  
  (async () => {
    let role: UserRole = "curator"; 
    let admin_logged = false;

    if (!username || !password || !email || !fullname || !affiliation) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    if (req.user) {
      const user = await Database.user.fromToken(req.user.jti);

      // If user is admin, allow all roles to be defined.
      if (user && user.role === 'admin') {
        admin_logged = true;
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

    if (!username.match(USERNAME_REGEX)) {
      return Errors.throw(ErrorType.InvalidUsername);
    }
    if (!email.match(EMAIL_REGEX)) {
      return Errors.throw(ErrorType.InvalidEmail);
    }

    const user: User = {
      id: generateSnowflake(),
      email,
      name: username,
      fullname,
      affiliation,
      role,
      created_at: new Date().toISOString(),
      password: "",
      approved: admin_logged,
    };

    // Register the user
    await Database.user.setOrChangePassword(user, password);

    if (user.approved) {
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
    }
    else {
      // Inform administrators (do not wait)
      informAdminFromAskCreation(user).catch(logger.error);

      res.json({
        created: sanitize({ ...user, password: null }),
        token: null
      });
    }
  })().catch(errorCatcher(res));
});

CreateUserRouter.all('/', methodNotAllowed('POST'));

export default CreateUserRouter;
