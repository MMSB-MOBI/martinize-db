import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import RadiusDatabase from '../../Entities/RadiusDatabase';
import { FORCE_FIELD_DIR } from '../../constants';
import JSZip from 'jszip';
import { promises as FsPromise } from 'fs';

const DownloadForceFieldRoute = Router();

DownloadForceFieldRoute.get('/list', (_, res) => {
  res.json(Object.keys(RadiusDatabase.FORCE_FIELD_TO_FILE_NAME).filter(e => e.startsWith('martini')));
});

DownloadForceFieldRoute.get('/download', (req, res) => {
  const name = req.query.name;

  (async () => {
    if (!name || typeof name !== 'string') {
      return Errors.throw(ErrorType.MissingParameters);
    }

    if (!(name in RadiusDatabase.FORCE_FIELD_TO_FILE_NAME)) {
      return Errors.throw(ErrorType.InvalidForceField);
    } 
    
    const filenames = RadiusDatabase.FORCE_FIELD_TO_FILE_NAME[name as string];
    const zip = new JSZip();

    for (const file of filenames) {
      zip.file(file, await FsPromise.readFile(FORCE_FIELD_DIR + file));
    }

    const zip_buffer = await zip.generateAsync({ compression: 'DEFLATE', compressionOptions: { level: 6 }, type: 'nodebuffer' });

    res
      .status(200)
      .type('zip')
      .header('Content-Disposition', 'attachement; filename=' + name + '.zip')
      .send(zip_buffer);
  })().catch(errorCatcher(res));
});

export default DownloadForceFieldRoute;
