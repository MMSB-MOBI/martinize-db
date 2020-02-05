import multer from 'multer';
import { UPLOAD_ROOT_DIR } from '../constants';

export const Uploader = multer({ dest: UPLOAD_ROOT_DIR });
export default Uploader;

// 100 MB
export const MAX_FILE_SIZE = 100 * 1024 * 1024;
