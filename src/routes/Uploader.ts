import multer from 'multer';
import { UPLOAD_ROOT_DIR } from '../constants';

export const Uploader = multer({ dest: UPLOAD_ROOT_DIR });
export default Uploader;

// 100 MB
export const MAX_FILE_SIZE = 100 * 1024 * 1024;

export const MAX_ITP_FILE_SIZE = 10 * 1024 * 1024;
export const MAX_PDB_FILE_SIZE = 100 * 1024 * 1024;

export const NAME_REGEX = /^[a-z0-9_ :()\.\/\\-]+$/i;
export const ALIAS_REGEX = /^[a-z0-9_ \.-]+$/i;
export const VERSION_REGEX = /^[0-9\.a-z_-]+$/i;
