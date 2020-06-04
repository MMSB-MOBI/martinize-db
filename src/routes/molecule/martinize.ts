import { Router } from 'express';
import { methodNotAllowed, cleanMulterFiles, errorCatcher, generateSnowflake, validateToken } from '../../helpers';
import Uploader from '../Uploader';
import { Martinizer, MartinizeSettings, ElasticOrGoBounds, GoMoleculeDetails } from '../../Builders/Martinizer';
import { SETTINGS_FILE } from '../../constants';
import { SettingsJson } from '../../types';
import { promises as FsPromise } from 'fs';
import Errors, { ErrorType, ApiError } from '../../Errors';
import shellescape from 'shell-escape';
import logger from '../../logger';
import path from 'path';
import { Database } from '../../Entities/CouchHelper';
import SocketIo from 'socket.io';
import TmpDirHelper from '../../TmpDirHelper';
import { Server } from 'http';

type MartinizeRunFailedPayload = { 
  error: string, 
  type: string, 
  stdout?: string, 
  stderr?: string, 
  dir: string, 
  code: ErrorType, 
  message: string, 
};

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

function numberOrError(num: any) {
  const pos = Number(num);
  if (isNaN(pos) || pos < 0) {
    return Errors.throw(ErrorType.Format);
  }
  return pos;
}

async function martinizeRun(parameters: any, pdb_path: string, onStep?: (step: string, ...data: any[]) => void) {
  const { ff, position, posref_fc, elastic, ef, el, eu, ea, ep, em, eb, use_go, sc_fix } = parameters;

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
    input: pdb_path,
  } as Partial<MartinizeSettings>;

  // Read all the settings
  if (ff) {
    if (settings.force_fields.includes(ff)) {
      runner.ff = ff;
    }
    else {
      return Errors.throw(ErrorType.InvalidForceField);
    }
  }
  if (position) {
    if (Martinizer.isMartinizePosition(position)) {
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
    const { pdb, itps, top } = await Martinizer.run(runner, onStep);

    // Create elastic if needed
    let elastic_bonds: ElasticOrGoBounds[] | undefined = undefined;

    if (runner.elastic) {
      elastic_bonds = await Martinizer.computeElasticNetworkBounds(top, itps);
    }

    return {
      pdb, 
      itps, 
      top, 
      elastic_bonds,
    };
  } catch (e) {
    if (e instanceof ApiError) {
      logger.error("Martinize Error.");
      logger.error(e.stack!);

      throw e;
    }
    if (e instanceof Error && e.stack?.startsWith('Error: Command failed')) {
      logger.error("Martinize Error.");
      logger.error(e.stack);

      return Errors.throw(ErrorType.MartinizeRunFailed);
    }

    // If any matches
    throw e;
  }
}

export function SocketIoMartinizer(app: Server) {
  const io = SocketIo(app);

  io.on('connection', socket => {
    socket.on('martinize', async (file: Buffer, run_id: string, settings: any) => {
      function sendFile(path: string, infos: { id?: string, name: string, type: string }) {
        return new Promise(async (resolve, reject) => {
          const timeout = setTimeout(reject, 1000 * 60 * 60);
          infos.id = run_id;

          socket.emit(
            'martinize download', 
            infos, 
            await FsPromise.readFile(path),
            () => {
              clearTimeout(timeout);
              resolve();
            }
          );
        });
      }

      if (!run_id || !file || !settings) {
        return;
      }
      if (run_id.length > 64) {
        return;
      }

      // // Verify token
      // try {
      //   await validateToken(token);
      // } catch {
      //   socket.emit('martinize error', { id: run_id, message: 'Invalid token.' });
      //   return;
      // }

      // Save to a temporary directory
      const tmp_dir = await TmpDirHelper.get();
      const INPUT = tmp_dir + '/input.pdb';

      try {
        await FsPromise.writeFile(INPUT, file);

        const { pdb, itps, top, elastic_bonds } = await martinizeRun(
          settings, 
          INPUT, 
          (step, ...data) => {
            socket.emit('martinize step', {
              id: run_id,
              step,
              data,
            });
          }
        );

        await sendFile(top, { 
          name: path.basename(top),
          type: 'chemical/x-topology' 
        });

        await sendFile(pdb, { 
          name: path.basename(pdb),
          type: 'chemical/x-pdb' 
        });

        for (const itp of itps) {
          await sendFile(itp, { 
            name: path.basename(itp),
            type: 'chemical/x-include-topology' 
          });
        }

        socket.emit('martinize before end', { id: run_id });

        const radius = await Database.radius.getRadius(
          settings.ff || 'martini22',
          itps,
        );

        socket.emit('martinize end', { id: run_id, elastic_bonds, radius });
      } catch (e) {
        // Error catch, test the error :D
        if (e instanceof ApiError && e.code === ErrorType.MartinizeRunFailed) {
          const { error, type, dir } = e.data as MartinizeRunFailedPayload;
          
          // Compress the directory
          const compressed_run = await Martinizer.zipDirectory(dir);

          socket.emit('martinize error', {
            id: run_id,
            error,
            type,
            stack: e.stack,
          }, compressed_run);
        }  
        else {
          socket.emit('martinize error', {
            id: run_id,
            error: e instanceof Error ? e.message : String(e),
            stack: e instanceof Error ? e.stack : "",
          });
        }
      } finally {
        await FsPromise.unlink(INPUT);
      }
    });
  });
}

MartinizerRouter.post('/', Uploader.single('pdb'), (req, res) => {
  const pdb_file = req.file;
  if (!pdb_file) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  (async () => {
    const { pdb, itps, top, elastic_bonds } = await martinizeRun(req.body, pdb_file.path);

    // Formatting itps
    const res_itp: { content: string, name: string, type: string, }[] = [];
    for (const itp of itps) {
      res_itp.push({
        content: await FsPromise.readFile(itp, 'utf-8'),
        name: path.basename(itp),
        type: 'chemical/x-include-topology',
      });
    }

    // todo: create the custom pdb? à voir

    // Get the radius for itps
    const radius = await Database.radius.getRadius(
      req.body.ff || 'martini22',
      itps,
    );

    res.json({
      pdb: { content: await FsPromise.readFile(pdb, 'utf-8'), name: path.basename(pdb), type: 'chemical/x-pdb' },
      itps: res_itp,
      top: { content: await FsPromise.readFile(top, 'utf-8'), name: path.basename(top), type: 'chemical/x-topology' },
      radius,
      elastic_bonds,
    });
  })().catch(errorCatcher(res));
});

MartinizerRouter.all('/', methodNotAllowed('POST'));

export default MartinizerRouter;
