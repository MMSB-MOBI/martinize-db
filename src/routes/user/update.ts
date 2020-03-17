import { Router } from 'express';
import { errorCatcher, methodNotAllowed, sanitize } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { USERNAME_REGEX, EMAIL_REGEX } from '../../constants';
import Mailer from '../../Mailer/Mailer';

const UpdateUserRouter = Router();

UpdateUserRouter.post('/', (req, res) => {
  (async () => {
    const { username, old_password, password, email, approved } = req.body;
    let usr_id = req.full_user!.id;
    const usr_role = req.full_user!.role;

    if (usr_role === "admin") {
      if (req.body.id) {
        usr_id = req.body.id;

        const exists = await Database.user.exists(usr_id);

        if (!exists) {
          return Errors.throw(ErrorType.UserNotFound);
        }
      }
    }

    let user = await Database.user.get(usr_id);
    let changed = false;
    let now_approved = false;

    if (username && typeof username === 'string') {
      const existant = await Database.user.fromUsername(username);
      if (existant && existant.id !== usr_id) {
        return Errors.throw(ErrorType.UsernameExists);
      }

      if (!username.match(USERNAME_REGEX)) {
        return Errors.throw(ErrorType.InvalidUsername);
      }

      user.name = username;
      changed = true;
    }
    if (email && typeof email === 'string') {
      const existant = await Database.user.fromEmail(email);
      if (existant && existant.id !== usr_id) {
        return Errors.throw(ErrorType.EmailExists);
      }
      if (!email.match(EMAIL_REGEX)) {
        return Errors.throw(ErrorType.InvalidEmail);
      }

      user.email = email;
      changed = true;
    }
    if (approved && usr_role === "admin") {
      if (approved === 'true' && !user.approved) {
        now_approved = true;
        user.approved = approved === 'true';
        changed = true;
      }
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

    if (now_approved) {
      // User is now approved, please send him an email
      await Mailer.send({ 
        to: user.email, 
        subject: "MArtinize Database - " + user.name + ": Your account has been approved" 
      }, 'mail_created', { 
        title: user.name + ": Your account has been approved",
        site_url: "http://localhost:3000",
        new_user: {
          name: user.name,
        },
      });
    }

    res.json(sanitize({ ...user, password: null }));
  })().catch(errorCatcher(res));
});

UpdateUserRouter.all('/', methodNotAllowed('POST'));

export default UpdateUserRouter;
