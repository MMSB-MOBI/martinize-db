import express from 'express';
import commander from 'commander';
import { VERSION, URLS, DEFAULT_TMP_BASE_DIR } from './constants';
import logger, { FORMAT_FILE } from './logger';
import Winston from 'winston';
import ApiRouter from './routes';
import Errors, { ErrorType, ApiError } from './Errors';
import { sendError } from './helpers';
import CouchHelper, { Database } from './Entities/CouchHelper';
import MOLECULE_CLI from './cli/molecule_cli';
import USER_CLI from './cli/user_cli';
import { CLI } from './cli/cli';
import MAIL_CLI from './cli/mail_cli';
import StaticServer from './static_server';
import CliHelper from 'interactive-cli-helper';
import DATABASE_CLI from './cli/databases_cli';
import TmpDirHelper from './TmpDirHelper';
import TEST_CLI from './cli/test.cli';
import { SocketIoMartinizer } from './routes/molecule/martinize';
import http from 'http';
import ShellManager from './Builders/ShellManager';

commander
  .version(VERSION)
  .option('-c, --couchdb-url <url>', 'Couch DB URL', String, process.env.COUCHDB_HOST || URLS.COUCH)
  .option('--server-url <url>', 'Server URL', String, process.env.SERVER_URL || URLS.SERVER)
  .option('-p, --port <port>', 'Emit port', Number, 4123)
  .option('--job-manager', 'Force using job manager as shell runner')
  .option('--os-tmp', 'Use automatic OSes temporary directory manager instead of ' + DEFAULT_TMP_BASE_DIR + ' base directory')
  .option('--child-process', 'Force using child process as shell runner')
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

// Log level
if (commander.logLevel) {
  logger.level = commander.logLevel;
}


// Shell handler mode
if (commander.jobManager && commander.childProcess) {
  logger.log("fatal", "You can't specify job manager AND child process as shell executor. Please select one of those.");
  process.exit(2);
}
else if (commander.jobManager) {
  ShellManager.mode = 'jm';
}
else if (commander.childProcess) {
  ShellManager.mode = 'child';
}

if (commander.osTmp) {
  TmpDirHelper.mode = 'os';
}
else {
  TmpDirHelper.mode = 'directory';
}

logger.silly(`Using ${ShellManager.mode === 'jm' ? 'job manager' : 'child processes'} as shell executor.`);
logger.silly(`Using ${TmpDirHelper.mode === 'os' ? 'os tmp dir manager' : 'custom tmp directory'} as base for creating temporary directories.`);


// Log files
if (commander.logFile) {
  logger.add(new Winston.transports.File({ 
      filename: commander.logFile, 
      level: commander.fileLogLevel, 
      eol: "\n", 
      format: FORMAT_FILE 
  }));
}


// CouchDB options
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


// Init options
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
logger.silly("Serving static website.");
// File should be in build/
app.use(StaticServer);

const HTTP_SERVER = http.createServer(app);

// Start socket.io
SocketIoMartinizer(HTTP_SERVER);


/* -------------------------- */
/* - COMMAND LINE INTERFACE - */
/* -------------------------- */

async function startCli() {
  const old_onclose = CLI.onclose.bind(CLI);

  CLI.onclose = async function() {
    await TmpDirHelper.clean();

    // this => attached to CLI; Should be fine
    old_onclose();
  };

  // Cli starter
  CLI.command('exit', async () => {
    await CLI.onclose();
    process.exit(0);
  });

  CLI.command('molecule', MOLECULE_CLI);
  CLI.command('user', USER_CLI);
  CLI.command('mail', MAIL_CLI);
  CLI.command('database', DATABASE_CLI);
  CLI.command('test', TEST_CLI);
  
  CLI.command(
    /^(\?|help)$/,  
    CliHelper.formatHelp("Martinize Database Server", {
      commands: {
        molecule: "Access and manage published / stashed molecules.",
        user: "Manage existing users, or create new ones.",
        worker: "View started search workers and kill existing instances.",
        mail: "Send test e-mails from defined templates.",
        database: "Create and wipe Couch databases.",
        exit: "Stop the server.",
      }
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

  HTTP_SERVER.listen(commander.port, () => {
    logger.info(`Martinize Database Server version ${VERSION} is listening on port ${commander.port}.`);
    startCli();
  });
}

main();

