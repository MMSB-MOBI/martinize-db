import fs from 'fs';
import logger from './logger';

const json_package = JSON.parse(fs.readFileSync(__dirname + '/../package.json', 'utf-8'));
export const VERSION = json_package.version;

let KEYS: { PUBLIC: string, PRIVATE: string };
try {
  KEYS = {
    PUBLIC: fs.readFileSync(__dirname + "/../.keys/key_new", "utf-8"),
    PRIVATE: fs.readFileSync(__dirname + "/../.keys/key_new.pem", "utf-8"),
  };
} catch (e) {
  logger.log("fatal", "It seems you didn't properly set the RSA keys. Check the documentation.");
  process.exit(0);
}

export { KEYS };
