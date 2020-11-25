import { Router } from 'express';
import { errorCatcher, methodNotAllowed, generateSnowflake } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import Mailer from '../../Mailer/Mailer';

const LostPasswordUserRouter = Router();

LostPasswordUserRouter.post('/', (req, res) => {
  (async () => {
    const email = req.body.email;

    if (!email) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    const user = await Database.user.fromEmail(email);

    if (!user) {
      return Errors.throw(ErrorType.UserNotFound);
    }
    if (!user.approved) {
      return Errors.throw(ErrorType.UserNotApproved);
    }

    // Generate a token
    user.lost_token = generateSnowflake();

    // Save user
    await Database.user.save(user);

    // Send an email
    Mailer.send({
      to: user.email,
      subject: "MArtini Database - Restore your password",
    }, 'mail_lost_password', {
      user: {
        name: user.name
      },
      token: user.lost_token
    });

    res.send();
  })().catch(errorCatcher(res));
});

LostPasswordUserRouter.all('/', methodNotAllowed('POST'));

export default LostPasswordUserRouter;
