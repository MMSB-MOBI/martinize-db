import jwt from 'express-jwt';
import { KEYS } from '../constants';
import { JSONWebToken } from '../types';
import nano from 'nano';
import { User } from '../Entities/entities';
import logger from '../logger';
import Errors, { ErrorType } from '../Errors';
import { getUserFromToken } from '../helpers';

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
    getUserFromToken(payload.jti)
      .then(user => {
        req.full_user = user;
        
        if (!user) {
          // Token is orphan !
          logger.error("Orphan token ! This should not happen. ", payload);
          done(null, true);
        }
        else if (!user.approved) {
          // User is not approved yet
          done(Errors.make(ErrorType.UserNotApproved), true);
        }
        else {
          done(null, false);
        }
      })
      .catch(() => done(null, true));
  }
}).unless(
  { path: [
    "/api/molecule/list", 
    "/api/user/login", 
    "/api/user/create", 
    "/api/user/lost_password", 
    "/api/user/change_password", 
    "/api", 
    { url: "/api/settings", methods: ['GET'] },
    "/api/settings/lipids",
    "/api/molecule",
    "/api/molecule/download",
    "/api/molecule/get",
    "/api/molecule/martinize",
    "/api/molecule/membrane_builder",
    "/api/user/contact",
    /^\/api\/molecule\/representation\/.+$/,
  ] }
);

// Extends Express request
declare module 'express-serve-static-core' {
  interface Request {
    user?: JSONWebToken;
    full_user?: (nano.DocumentGetResponse & User);
  }
}
