import multer from 'multer';
import { UPLOAD_ROOT_DIR } from '../constants';

export const Uploader = multer({ dest: UPLOAD_ROOT_DIR });
export default Uploader;
