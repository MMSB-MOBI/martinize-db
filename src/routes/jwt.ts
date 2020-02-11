import jwt from 'express-jwt';
import { KEYS } from '../constants';
import { Database } from '../Entities/CouchHelper';
import { JSONWebToken } from '../types';
import nano from 'nano';
import { User } from '../Entities/entities';
import logger from '../logger';

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
  isRevoked: (req, payload, done) => {
    // Get the token from string and call done(null, is_revoked)
    Database.token.get(payload.jti as string)
      .then(() => Database.user.fromToken(payload.jti as string))
      .then(user => {
        req.full_user = user;
        
        if (!user) {
          // Token is orphan !
          logger.error("Orphan token ! This should not happen. ", payload);
          done(null, true);
        }
        else {
          done(null, false);
        }
      })
      .catch(() => done(null, true));
  }
}).unless(
  { path: ["/api/molecule/list", "/api/user/login", "/api/user/create", "/api", { url: "/api/settings", methods: ['GET'] }] }
);

// Extends Express request
declare module 'express-serve-static-core' {
  interface Request {
    user?: JSONWebToken;
    full_user?: (nano.DocumentGetResponse & User);
  }
}
