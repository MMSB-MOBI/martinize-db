import { Router } from 'express';
import { errorCatcher, methodNotAllowed, informAdminContact } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { EMAIL_REGEX } from '../../constants';

const ContactRouter = Router();

ContactRouter.post('/', (req, res) => {
  (async () => {
    const { content, email } = req.body;

    if (!content || typeof content !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }
    if (!email || typeof email !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    if (!email.match(EMAIL_REGEX)) {
      return Errors.throw(ErrorType.InvalidEmail);
    }

    await informAdminContact(content, email);

    res.send();
  })().catch(errorCatcher(res));
});

ContactRouter.all('/', methodNotAllowed('POST'));

export default ContactRouter;
