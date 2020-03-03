import express from 'express';
import commander from 'commander';
import { VERSION, URLS } from './constants';
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
  .option('-c, --couchdb-url <url>', 'Couch DB URL', String, process.env.COUCHDB_HOST || URLS.COUCH)
  .option('--server-url <url>', 'Server URL', String, process.env.SERVER_URL || URLS.SERVER)
  .option('-p, --port <port>', 'Emit port', Number, 4123)
  .option('--wipe-init')
  .option('--init-db')
  .option('--quit-after-init')
  .option('-l, --log-level <logLevel>', 'Log level [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--file-log-level <logLevel>', 'Log level (written to file) [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--log-file <logFile>', 'File log level')
.parse(process.argv);

const app = express();


/* ------------------ */
/* - PARSE CLI ARGS - */
/* ------------------ */

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
  let url = commander.couchdbUrl;

  if (!url.startsWith('http://')) {
    if (process.env.COUCHDB_USER) {
      url = 'http://' + process.env.COUCHDB_USER + ':' + process.env.COUCHDB_PASSWORD + '@' + url;
    }
    else {
      url = 'http://' + url;
    }
  }

  Database.refresh(url);
  URLS.COUCH = url;
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


/* ------------------------------ */
/* - STARTING EXPRESS ENDPOINTS - */
/* ------------------------------ */

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


/* -------------------------- */
/* - COMMAND LINE INTERFACE - */
/* -------------------------- */

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
      database: "Create and wipe Couch databases.",
      exit: "Stop the server.",
    })
  );

  console.log("\nWelcome to Martinize server CLI. For help, type \"help\".");

  const db_exists = await Database.link.use(CouchHelper.USER_COLLECTION).info().catch(e => ({ not_found: true }));
  if ('not_found' in db_exists) {
    console.log("\nWARN: The database seems to be un-initialized. Please create all databases by entering \"database create all\".");
    console.log("WARN: Once database is created, you can create an administrator account with \"user create\".");
  }
  else {
    const user_db = await Database.user.find({ selector: { role: 'admin' } });
    if (!user_db.length) {
      console.log("\nWARN: Server doesn't seem to have an administrator account created. You can create an user with \"user create\".");
    }
  }

  console.log("");

  CLI.listen();
}


/* -------------------------------------- */
/* - HANDLE UNCATCHED REJECTED PROMISES - */
/* -------------------------------------- */

function propertiesValues(obj: any) {
  const data = Object.getOwnPropertyDescriptors(obj);

  for (const key in data) {
    data[key] = data[key].value;
  }

  return data;
}

process.on('unhandledRejection', reason => {
  const maximum_detail = typeof reason === 'object' && reason !== null ? JSON.stringify(propertiesValues(reason), null, 2)
    : (reason ? reason : "No rejection content.");

  logger.error("Unhandled rejected Promise handled: \n" + String(maximum_detail));
});


/* ------------------------------------------------- */
/* - STARTING THE SERVER AND LISTENING TO REQUESTS - */
/* ------------------------------------------------- */

async function main() {
  try {
    await Database.ping();
  } catch (e) {
    logger.error("CouchDB is not running or is unreachable. You must start Couch or specify a valid database URL.");
    console.log("Stack trace:", 'stack' in e ? e.stack : e);
    process.exit(2);
  }

  app.listen(commander.port, () => {
    logger.info(`Martinize Database Server version ${VERSION} is listening on port ${commander.port}.`);
    startCli();
  });
}

main();


/** MISC */
// Nano Type definition for MangoSelector is incorrect. If any TSC warning, remove type `MangoSelector` in `nano.d.ts`.
// This should be soonly corrected, see https://github.com/apache/couchdb-nano/issues/211
declare module 'nano' {
  type MangoSelector = {
    [K in MangoOperator | string]: MangoSelector | MangoValue | MangoValue[];
  }
}
