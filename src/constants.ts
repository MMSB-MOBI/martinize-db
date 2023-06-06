import fs from 'fs';
import logger from './logger';
import SMTPTransport from 'nodemailer/lib/smtp-transport';
import path from 'path';

const json_package = JSON.parse(fs.readFileSync(__dirname + '/../package.json', 'utf-8'));
export const VERSION = json_package.version as string;

/**
 * OpenSSL's SSH RSA keys used to make/decode JSON Web Tokens (JWT)
 */
let KEYS: { PUBLIC: string, PRIVATE: string };
try {
  KEYS = {
    PUBLIC: fs.readFileSync(__dirname + "/../.keys/key_new", "utf-8"),
    PRIVATE: fs.readFileSync(__dirname + "/../.keys/key_new.pem", "utf-8"),
  };
} catch (e) {
  logger.log("fatal", "It seems you didn't properly set the RSA keys, the server can't run. Check the documentation.");
  process.exit(0);
}

export { KEYS };

/* - DEFAULT DIRECTORIES - */
export const MOLECULE_ROOT_DIR    = process.env.MOLECULE_ROOT_DIR ?? path.resolve(__dirname, "../molecules/")
export const MARTINIZER_ROOT_DIR  = process.env.MARTINIZER_ROOT_DIR ?? path.resolve(__dirname, "../molecules/martinizer/")
export const UPLOAD_ROOT_DIR      = process.env.UPLOAD_ROOT_DIR ?? path.resolve(__dirname, "../uploads/")

export const LIPIDS_ROOT_DIR      = process.env.LIPIDS_ROOT_DIR ?? path.resolve(__dirname, "../lipids/")
export const SETTINGS_FILE        = path.resolve(__dirname, "../settings.json");
export const TEMPLATE_DIR         = path.resolve(__dirname, "../templates/") + "/";
export const FORCE_FIELD_DIR      = process.env.FORCE_FIELD_DIR ?? path.resolve(__dirname, "../force_fields/")
export const DEFAULT_TMP_BASE_DIR = process.env.DEFAULT_TMP_BASE_DIR ?? path.resolve(__dirname, "../tmp/")
export const HISTORY_ROOT_DIR = process.env.HISTORY_DIR ?? path.resolve(__dirname, "../history")

/* - Couch database - */
export const DB_PREFIX = process.env.DB_PREFIX ?? ""



/* - Job manager - */
export type JobMethod = 'jm' | 'child';
export const DEFAULT_JOB_METHOD: JobMethod = 'jm';

/* - SEARCH WORKERS SETTINGS - */
export const MINUTES_BEFORE_WORKER_KILL = 2;
export const MAX_POOL_SIZE = 10;
export const MAX_REQUEST_PER_WORKER_THRESHOLD = 3;

/** If `true`, every sended molecule to `/api/molecule/create` will be directly accepted without moderation. */
export const DISABLE_MODERATION_PROCESS = false;

/** Name of the e-mail sender. */
export const DEFAULT_MAILER_NAME = "MArtini Database";
/** E-mail address of the mail sender. */
export const DEFAULT_MAILER_ADDRESS = "martinize.db@ibcp.fr";
/** Debug purpose only; If `string`, all e-mails will be sent to the following address. */
export const MAILER_ENFORCE_RECIPIENT: false | string = false;
/** Default parameters for Mailer. See `nodemailer` package documentation. */
export const MAILER_TRANSPORT_SETTINGS: SMTPTransport.Options = {
  host: 'smtp.ibcp.fr',
  port: 587,
  secure: false, // true for 465, false for other ports
  tls: {
    // do not fail on invalid certs
    rejectUnauthorized: false
  }
};
/** Default parameters for Mailer. See `ms-jobmanager` package documentation. */
export const JOB_MANAGER_SETTINGS = {
  address: process.env.JM_ADRESS ?? "localhost",
  port: parseInt(process.env.JM_PORT || "") ?? "1234"
};

/* - Martinizer constants - */
/** DSSP Path. For now, due to a bug in Martinize2, it's optional. */
export const DSSP_PATH = process.env.DSSP_PATH; //if undefined => don't use dssp
/** Link to script used to create go virtual sites launcher. */
export const CREATE_GO_PATH = path.resolve(__dirname, "../utils/create_go_virt.sh");
export const CREATE_GO_PATH_JM = path.resolve(__dirname, "../utils/create_go_virt_coreScript.sh");
/** Link to script used to create go virtual sites. */
export const CREATE_GO_PY_SCRIPT_PATH = path.resolve(__dirname, "../utils/create_goVirt.py");
/** Link to script that can start ccmap */
export const CREATE_MAP_PATH = path.resolve(__dirname, "../utils/get_map.sh");
export const CREATE_MAP_PATH_JM = path.resolve(__dirname, "../utils/get_map_coreScript.sh");
export const CREATE_MAP_RCSU_PATH = path.resolve(__dirname, "../utils/get_map_rcsu.sh"); 
/** Associated python script to ccmap */
export const CREATE_MAP_PY_SCRIPT_PATH = path.resolve(__dirname, "../utils/get_map.py");
/** Link to script that can run GROMACS */
export const CONECT_PDB_PATH = path.resolve(__dirname, "../utils/create_conect_pdb.sh");
export const MINIMIZEPDB = path.resolve(__dirname, "../utils/minimize_and_create_pdb.sh");
export const CONECT_PDB_PATH_JM = path.resolve(__dirname, "../utils/create_conect_pdb_coreScript.sh");
export const MINIMIZEPDBBIS = path.resolve(__dirname, "../utils/minimize_and_create_pdb_coreScript.sh");
/** Link to MDP file needed for GROMACS's grompp */
export const CONECT_MDP_PATH = process.env.CONECT_MDP_PATH ?? path.resolve(__dirname, "../force_fields/run.mdp");
/** Link to Python 3 binary */
export const PYTHON_3_PATH = "python";
/** Path to script that starts martinize2 */
export const MARTINIZE_PATH = path.resolve(__dirname, "../utils/martinize.sh");
export const MARTINIZE_PATH_JM = process.env.MARTINIZE_PATH_JM ?? path.resolve(__dirname, "../utils/martinize_coreScript.sh"); //need to provide a script with chmod hack for old infra


export const RUN_POLYPLY_PATH = path.resolve(__dirname, "../utils/run_polyply.sh");
export const INIT_POLYPLY_PATH = path.resolve(__dirname, "../utils/init_polyply.sh");
export const POLYPLY_VENV = process.env.POLYPLY_VENV ?? path.resolve(__dirname, "/polyply_1.0/venv/bin/activate");
export const POLYPLYPATHDATA  = process.env.POLYPLYPATHDATA ?? path.resolve(__dirname, "../data"); 

/* - Membrane builder constants - */
/** Full path to insane start script */
export const INSANE_PATH = path.resolve(__dirname, "../utils/insane.sh");
export const INSANE_PATH_JM = path.resolve(__dirname, "../utils/insane_coreScript.sh");
export const INSANE_HACK_SCRIPT = {
  BEFORE: path.resolve(__dirname, "../utils/insane_hack.py"),  
  AFTER: path.resolve(__dirname, "../utils/insane_hack_reverse.py")
} //scripts to hack pdb coordinates when we have nan atoms



//Venv for local computation
export const MARTINIZE_VENV = process.env.MARTINIZE_VENV ?? path.resolve(__dirname, "../martinize2venv/bin/activate");
export const INSANE_VENV = process.env.INSANE_VENV ?? path.resolve(__dirname, "../insanevenv/bin/activate")

export const RCSU_PATH = process.env.RCSU_PATH


/**
 * Default URLs.
 * - `SERVER` is the public URL of the server. Don't forget to set it in order to have working URLs in e-mails !
 * - `COUCH` is default CouchDB URL. Usually, this URL is not used, the `--couchdb-url` parameter of server is used instead.
 */
export const URLS = {
  SERVER: process.env.SERVER_URL ?? "http://localhost:3003", //3003
  COUCH: process.env.COUCH_URL ?? "http://localhost:5984",
};

/* - Regular expressions used to check recieved parameters - */
export const USERNAME_REGEX = /^[a-z][a-z0-9_-]*[a-z0-9]$/i;
export const EMAIL_REGEX = /^(?:[a-z0-9!#$%&'*+/=?^_`{|}~-]+(?:\.[a-z0-9!#$%&'*+/=?^_`{|}~-]+)*|"(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21\x23-\x5b\x5d-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])*")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+[a-z0-9](?:[a-z0-9-]*[a-z0-9])?|\[(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?|[a-z0-9-]*[a-z0-9]:(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21-\x5a\x53-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])+)\])$/i;

// Slurm profile to use for job manager
export const SLURM_PROFILES = {
  JOB_PROFILE : process.env.JOB_PROFILE ?? "", 
  SYS_SETTINGS : process.env.JOB_SYS_SETTINGS ?? ""
}

export const MARTINIZE_VERSION = process.env.MARTINIZE_VERSION ?? "unknown"

//martinize-db/data/polyply-env$ polyply -V
// export const POLYPLY_VERSION = "1.5.0"

export const SEND_COMPLETION_MAIL = process.env.SEND_COMPLETION_MAIL === 'false' ? false : true; 

export const MAINTENANCE = {mode : false}; 
