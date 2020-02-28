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
export const TEMPLATE_DIR = __dirname + "/../templates/";
export const COUCH_URL = "http://localhost:5984";
export const MINUTES_BEFORE_WORKER_KILL = 2;
export const MAX_POOL_SIZE = 10;
export const MAX_REQUEST_PER_WORKER_THRESHOLD = 3;
export const DISABLE_MODERATION_PROCESS = false;
export const SERVER_URL = "http://localhost:4123";

// Regex user
export const USERNAME_REGEX = /^[a-z][a-z0-9_-]*[a-z0-9]$/i;
export const EMAIL_REGEX = /^(?:[a-z0-9!#$%&'*+/=?^_`{|}~-]+(?:\.[a-z0-9!#$%&'*+/=?^_`{|}~-]+)*|"(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21\x23-\x5b\x5d-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])*")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+[a-z0-9](?:[a-z0-9-]*[a-z0-9])?|\[(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?|[a-z0-9-]*[a-z0-9]:(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21-\x5a\x53-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])+)\])$/i;
