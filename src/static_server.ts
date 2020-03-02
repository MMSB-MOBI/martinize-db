import path from 'path';
import express from 'express';

const StaticServer = express.Router();

const static_serve_middleware = express.static(path.join(__dirname, "../static/www"));

StaticServer.use('*', static_serve_middleware);

export default StaticServer;
