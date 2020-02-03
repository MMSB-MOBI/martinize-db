import express from 'express';
import commander from 'commander';
import { VERSION } from './constants';
import logger, { FORMAT_FILE } from './logger';
import Winston from 'winston';
import ApiRouter from './routes';
import CliHelper from 'interactive-cli-helper';

commander
  .version(VERSION)
  .option('-c, --couchdb-url <url>', 'Couch DB URL', Number, 5984)
  .option('-p, --port <port>', 'Emit port', Number, 4123)
  .option('-l, --log-level <logLevel>', 'Log level [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--file-log-level <logLevel>', 'Log level (written to file) [debug|silly|verbose|info|warn|error]', /^(debug|silly|verbose|info|warn|error)$/, 'info')
  .option('--log-file <logFile>', 'File log level')
.parse(process.argv);

const app = express();

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

app.use('/api', ApiRouter);

app.listen(commander.port, () => {
  logger.info(`Martinize Database Server version ${VERSION} is listening on port ${commander.port}.`);
  startCli();
});

function startCli() {
  // Cli starter
  const CLI = new CliHelper("Command not found.");

  CLI.addSubListener('exit', () => {
    CLI.onclose!();
    process.exit(0);
  });

  CLI.listen();
}

