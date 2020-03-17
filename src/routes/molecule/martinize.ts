import { Router } from 'express';
import { methodNotAllowed, cleanMulterFiles, errorCatcher } from '../../helpers';
import Uploader from '../Uploader';
import { Martinizer, MartinizeSettings } from '../../Martinizer/Martinizer';
import { SETTINGS_FILE } from '../../constants';
import { SettingsJson } from '../../types';
import { promises as FsPromise } from 'fs';
import Errors, { ErrorType } from '../../Errors';
import shellescape from 'shell-escape';
import logger from '../../logger';
import path from 'path';
import { Database } from '../../Entities/CouchHelper';

const MartinizerRouter = Router();

// Middleware that wipe uploaded files after request
MartinizerRouter.use((req, res, next) => {
  function after() {
    // Response is sended
    cleanMulterFiles(req);
    res.removeListener('finish', after);
  }

  res.once('finish', after);
  next();
});

MartinizerRouter.post('/', Uploader.single('pdb'), (req, res) => {
  function numberOrError(num: any) {
    const pos = Number(num);
    if (isNaN(pos) || pos < 0) {
      return Errors.throw(ErrorType.Format);
    }
    return pos;
  }

  const pdb_file = req.file;
  if (!pdb_file) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  const { ff, position, posref_fc, elastic, ef, el, eu, ea, ep, em, eb, use_go, sc_fix } = req.body;

  (async () => {
    const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));

    /// NON PRIS EN CHARGE (TODO) ///
    /**
     *  ignore?: string[];
     *  ignh?: boolean;
     *  collagen?: boolean;
     *  dihedral?: boolean;
     *  neutral_termini?: boolean;
     *  cystein_bridge?: string;
     */

    const runner = {
      input: pdb_file.path,
    } as Partial<MartinizeSettings>;

    if (ff) {
      if (settings.force_fields.includes(ff)) {
        runner.ff = ff;
      }
      else {
        return Errors.throw(ErrorType.InvalidForceField);
      }
    }
    if (position) {
      if (['none', 'all', 'backbone'].includes(position)) {
        runner.position = position;
      }
      else {
        return Errors.throw(ErrorType.Format);
      }
    }
    if (posref_fc) {
      runner.posref_fc = numberOrError(posref_fc);
    }
    if (elastic === "true") {
      runner.elastic = true;
    }
    if (ef) {
      runner.ef = numberOrError(ef);
    }
    if (el) {
      runner.el = numberOrError(el);
    }
    if (eu) {
      runner.eu = numberOrError(eu);
    }
    if (ea) {
      runner.ea = numberOrError(ea);
    }
    if (ep) {
      runner.ep = numberOrError(ep);
    }
    if (em) {
      runner.em = numberOrError(em);
    }
    if (eb && typeof eb === 'string') {
      const arg = shellescape([eb]);
      runner.eb = arg.split(',');
    }
    if (use_go === "true") {
      runner.use_go_virtual_sites = true;
    }
    if (sc_fix === "true") {
      runner.side_chain_fix = true;
    }

    try {
      const { pdb, itps, top } = await Martinizer.run(runner);
  
      // Formatting itps
      const res_itp: { content: string, name: string }[] = [];
      for (const itp of itps) {
        res_itp.push({
          content: await FsPromise.readFile(itp, 'utf-8'),
          name: path.basename(itp),
        });
      }

      // Get the radius for itps
      const radius = await Database.radius.getRadius(runner.ff || 'martini22', res_itp.map(e => e.content));

      res.json({
        pdb: { content: await FsPromise.readFile(pdb, 'utf-8'), name: path.basename(pdb) },
        itps: res_itp,
        top: { content: await FsPromise.readFile(top, 'utf-8'), name: path.basename(top) },
        radius
      });
    } catch (e) {
      if (e instanceof Error && e.stack?.startsWith('Error: Command failed')) {
        logger.error("Martinize Error.");
        logger.error(e.stack);

        return Errors.throw(ErrorType.MartinizeRunFailed);
      }
      throw e;
    }
  })().catch(errorCatcher(res));
});

MartinizerRouter.all('/', methodNotAllowed('POST'));

export default MartinizerRouter;
