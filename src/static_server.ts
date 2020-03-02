import path from 'path';
import express from 'express';

const StaticServer = express.Router();

function serveHtmlMiddleware(_: express.Request, res: express.Response, __: express.NextFunction) {
  // Aucune modification nécessaire
  res.sendFile(path.join(__dirname, "../static/www/index.html"));
}

const static_serve_middleware = express.static(path.join(__dirname, "../static/www"));

StaticServer.use('/', (req, res, next) => {
  if (req.path === "/" || req.path === "/index.html") {
      serveHtmlMiddleware(req, res, next);
  }
  else {
      static_serve_middleware(req, res, next);
  }
});

StaticServer.use('*', serveHtmlMiddleware);

export default StaticServer;
