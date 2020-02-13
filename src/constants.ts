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

export const MOLECULE_ROOT_DIR = __dirname + "/../molecules/";
export const UPLOAD_ROOT_DIR = __dirname + "/../uploads/";
export const SETTINGS_FILE = __dirname + "/../settings.json";
export const COUCH_URL = "http://localhost:5984";
export const MINUTES_BEFORE_WORKER_KILL = 2;
export const MAX_POOL_SIZE = 10;
export const MAX_REQUEST_PER_WORKER_THRESHOLD = 3;
