import jwt from 'express-jwt';
import { KEYS } from '../constants';
import { Database } from '../Entities/CouchHelper';
import { JSONWebToken } from '../types';

export default jwt({
  getToken: function fromHeaderOrQuerystring(req) {
    if (req.headers.authorization && req.headers.authorization.split(' ')[0] === 'Bearer') {
      return req.headers.authorization.split(' ')[1];
    } else if (req.cookies && req.cookies.login_token) {
      return req.cookies.login_token;
    }
    return null;
  },
  secret: KEYS.PUBLIC,
  credentialsRequired: true,
  isRevoked: (_, payload, done) => {
    // Get the token from string and call done(null, is_revoked)
    Database.token.get(payload.jti as string)
      .then(() => done(null, false))
      .catch(() => done(null, true));
  }
}).unless(
  { path: ["/api/molecule/list", "/api/user/login", "/api"] }
);

// Extends Express request
declare module 'express-serve-static-core' {
  interface Request {
    user?: JSONWebToken;
  }
}
