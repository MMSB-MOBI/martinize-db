import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import Mailer from '../../Mailer/Mailer';

const ChangePasswordUserRouter = Router();

ChangePasswordUserRouter.post('/', (req, res) => {
  (async () => {
    const password = req.body.password;
    const token = req.body.token;

    if (!password || !token) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const users = await Database.user.find({
      selector: {
        lost_token: token
      }
    });

    if (!users.length) {
      return Errors.throw(ErrorType.UserNotFound);
    }

    const user = users[0];
    if (!user.approved) {
      return Errors.throw(ErrorType.UserNotApproved);
    }

    // Reset lost_token
    user.lost_token = undefined;

    // Save user
    await Database.user.setOrChangePassword(user, password);

    // Send an email
    Mailer.send({
      to: user.email,
      subject: "MArtini Database - Password changed",
    }, 'mail_changed_password', {
      user: {
        name: user.name
      },
    });

    res.send();
  })().catch(errorCatcher(res));
});

ChangePasswordUserRouter.all('/', methodNotAllowed('POST'));

export default ChangePasswordUserRouter;
