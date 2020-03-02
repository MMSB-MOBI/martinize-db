import express from 'express';
import commander from 'commander';
import { VERSION } from './constants';
import logger, { FORMAT_FILE } from './logger';
import Winston from 'winston';
import ApiRouter from './routes';
import Errors, { ErrorType, ApiError } from './Errors';
import { sendError } from './helpers';
import CouchHelper, { Database } from './Entities/CouchHelper';
import MOLECULE_CLI from './cli/molecule_cli';
import USER_CLI from './cli/user_cli';
import WORKER_CLI from './cli/worker_cli';
import { CLI } from './cli/cli';
import MAIL_CLI from './cli/mail_cli';
import StaticServer from './static_server';
import CliHelper from 'interactive-cli-helper';
import DATABASE_CLI from './cli/databases_cli';

commander
  .version(VERSION)
  .option('-c, --couchdb-url <url>', 'Couch DB URL', String, 'http://localhost:5984')
  .option('-p, --port <port>', 'Emit port', Number, 4123)
  .option('--wipe-init')
  .option('--init-db')
  .option('--quit-after-init')
  .option('-l, --log-level <logLevel>', 'Log level [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--file-log-level <logLevel>', 'Log level (written to file) [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--log-file <logFile>', 'File log level')
.parse(process.argv);

const app = express();

// Parse CLI args
if (commander.logLevel) {
  logger.level = commander.logLevel;
}

if (commander.logFile) {
  logger.add(new Winston.transports.File({ 
      filename: commander.logFile, 
      level: commander.fileLogLevel, 
      eol: "\n", 
      format: FORMAT_FILE 
  }));
}

if (commander.couchdbUrl) {
  Database.refresh(commander.couchdbUrl);
}

if (commander.wipeInit) {
  logger.info("Wiping databases and creating them again");
  Database.wipeAndCreate()
    .then(() => {
      if (commander.quitAfterInit) {
        logger.info("Exiting.");
        process.exit(0);
      }
    })
}

if (commander.initDb) {
  logger.info("Creating all databases");
  Database.createAll()
    .then(() => {
      if (commander.quitAfterInit) {
        logger.info("Exiting.");
        process.exit(0);
      }
    })
}

// Register API router
app.use('/api', ApiRouter);

// Catch API errors
app.use('/api', (err: any, req: express.Request, res: express.Response, next: Function) => {
  if (res.headersSent) {
    next(err);
    return;
  }

  if (err.name === 'UnauthorizedError') {
    logger.debug("Token identification error: " + err.name);
    Errors.send(ErrorType.TokenInvalid, res);
  }
  else if (err instanceof ApiError) {
    sendError(err, res);
  }
  // @ts-ignore Invalid field in request
  else if (req.field) {
    // @ts-ignore
    Errors.send(ErrorType.Format, res, { field: req.field });
  }
  else {
    next(err);
  }
});

// Serve the static folder
logger.debug("Serving static website");
// File should be in build/
app.use(StaticServer);


async function startCli() {
  // Cli starter
  CLI.addSubListener('exit', () => {
    CLI.onclose!();
    process.exit(0);
  });

  CLI.addSubListener('molecule', MOLECULE_CLI);
  CLI.addSubListener('user', USER_CLI);
  CLI.addSubListener('worker', WORKER_CLI);
  CLI.addSubListener('mail', MAIL_CLI);
  CLI.addSubListener('database', DATABASE_CLI);
  
  CLI.addSubListener(
    /^(\?|help)$/,  
    CliHelper.formatHelp("Martinize Database Server", {
      molecule: "Access and manage published / stashed molecules.",
      user: "Manage existing users, or create new ones.",
      worker: "View started search workers and kill existing instances.",
      mail: "Send test e-mails from defined templates.",
      exit: "Stop the server.",
    })
  );

  console.log("Welcome to Martinize server CLI. For help, type \"help\".");

  const db_exists = await Database.link.use(CouchHelper.USER_COLLECTION).info().catch(e => ({ not_found: true }));
  if ('not_found' in db_exists) {
    console.log("The database seems to be un-initialized. Please create all the databases by entering \"database create all\".");
  }

  CLI.listen();
}

async function main() {
  await Database.ping();

  app.listen(commander.port, () => {
    logger.info(`Martinize Database Server version ${VERSION} is listening on port ${commander.port}.`);
    startCli();
  });
}

main();
