import jwt from 'express-jwt';
import { KEYS } from '../constants';

export default jwt({
  getToken: function fromHeaderOrQuerystring(req) {
    if (req.headers.authorization && req.headers.authorization.split(' ')[0] === 'Bearer') {
      return req.headers.authorization.split(' ')[1];
    } else if (req.cookies && req.cookies.login_token) {
      return req.cookies.login_token;
    }
    return null;
  },
  secret: KEYS.PRIVATE,
  credentialsRequired: true,
  isRevoked: (req, payload, done) => {
    isTokenInvalid(payload.jti, req, payload)
      .then(is_revoked => { done(null, is_revoked); })
      .catch(e => { logger.error("Unable to check token validity", e); done(e); });
  }
}).unless(
  { path: ["/api/user/request.json", "/api/users/access.json", "/api/callback_twitter", "/api", "/api/deleted_count.json"] }
);
