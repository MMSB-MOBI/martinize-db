import { Router } from 'express';
import { errorCatcher, methodNotAllowed, sanitize } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';

const UpdateUserRouter = Router();

UpdateUserRouter.post('/', (req, res) => {
  (async () => {
    const { username, old_password, password, email } = req.body;
    const usr_id = req.full_user!.id;

    let user = await Database.user.get(usr_id);
    let changed = false;

    if (username && typeof username === 'string') {
      // TODO check everything
      const existant = await Database.user.fromUsername(username);
      if (existant && existant.id !== usr_id) {
        return Errors.throw(ErrorType.UsernameExists);
      }

      user.name = username;
      changed = true;
    }
    if (email && typeof email === 'string') {
      const existant = await Database.user.fromEmail(email);
      if (existant && existant.id !== usr_id) {
        return Errors.throw(ErrorType.EmailExists);
      }

      user.email = email;
      changed = true;
    }
    if (old_password && password && typeof password === 'string' && typeof old_password === 'string') {
      const is_valid = await Database.user.verifyPassword(user, old_password);

      if (!is_valid) {
        return Errors.throw(ErrorType.InvalidPassword);
      }

      // @ts-ignore
      user = await Database.user.setOrChangePassword(user, password);
      changed = false;
    }

    if (changed) {
      await Database.user.save(user);
    }

    res.json(sanitize({ ...user, password: null }));
  })().catch(errorCatcher(res));
});

UpdateUserRouter.all('/', methodNotAllowed('POST'));

export default UpdateUserRouter;
