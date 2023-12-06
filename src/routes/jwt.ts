//import jwt from 'express-jwt';
//import { expressjwt, ExpressJwtRequest } from "express-jwt";
import { expressjwt } from 'express-jwt'
import { KEYS } from '../constants';
import { JSONWebToken } from '../types';
import nano from 'nano';
import { User } from '../Entities/entities';
import logger from '../logger';
import Errors, { ErrorType } from '../Errors';
import { getUserFromToken } from '../helpers';
import { inspect } from 'util';
export default expressjwt({
  getToken: function fromHeaderOrQuerystring(req:any) {
    if (req.headers.authorization && req.headers.authorization.split(' ')[0] === 'Bearer') {
      return req.headers.authorization.split(' ')[1];
    } else if (req.cookies && req.cookies.login_token) {
      return req.cookies.login_token;
    }
    return null;
  },
  secret: KEYS.PUBLIC,
  algorithms: ["RS256"],
  credentialsRequired: true,
  isRevoked: async (req:any, packet:any) => {
    try {
    // Get the token from string and call done(null, is_revoked)
      logger.info(`[expressjwt::isRevoked] user request:\n${inspect(packet)}`);
      const user = await getUserFromToken(packet.payload.jti);
      logger.info(`[expressjwt::isRevoked] getUserFromToken value:\n${inspect(user)}`);
      req.full_user = user;
      if (!user) {
        // Token is orphan !
        logger.error("[expressjwt::isRevoked] Orphan token ! This should not happen. ", packet);
        return false;
      }
      if (!user.approved) {
    // User is not approved yet
        return false;
      //done(Errors.make(ErrorType.UserNotApproved), true);
      }
      logger.info("[expressjwt::isRevoked] iSRevoked returns false (access granted)");
      return false;
    } catch (e:any){
      logger.warn(`[expressjwt::isRevoked] error (access denied) : ${e}`);
      return true;
    }
  }
}).unless(
  {
    path: [
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
      "/api/molecule/martinize",
      "/api/molecule/membrane_builder",
      "/api/molecule/get",
      "/api/user/contact",
      "/api/force_fields/list",
      "/api/force_fields/download",
      "/api/polymergenerator/data",
      /^\/api\/molecule\/representation\/.+$/,
      /^\/api\/molecule\/get\/.+$/
    ]
  }
);

// Extends Express request
declare module 'express-serve-static-core' {
  interface Request {
    user?: JSONWebToken;
    full_user?: (nano.DocumentGetResponse & User);
  }
}
