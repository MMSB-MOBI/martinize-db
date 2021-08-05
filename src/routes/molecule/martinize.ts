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
import { History } from '../../Entities/entities'; 
import SocketIo from 'socket.io';
import TmpDirHelper from '../../TmpDirHelper';
import { Server } from 'http';
import ShellManager, { JobInputs } from '../../Builders/ShellManager';
import HistoryOrganizer from "../../HistoryOrganizer";

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

function createRunner(settings: any, parameters: any, pdb_path? : string) {

  let inp = '';
  if (pdb_path) {
    inp = pdb_path;
  } else {
    inp = 'input'
  }
  const runner = {
    input: inp,
  } as Partial<MartinizeSettings>;

    // Read all the settings
    if (parameters.ff) {
      if (settings.force_fields.includes(parameters.ff)) {
        runner.ff = parameters.ff;
      }
      else {
        return Errors.throw(ErrorType.InvalidForceField);
      }
    }
    if (parameters.position) {
      if (Martinizer.isMartinizePosition(parameters.position)) {
        runner.position = parameters.position;
      }
      else {
        return Errors.throw(ErrorType.Format);
      }
    }
    if (parameters.posref_fc) {
      runner.posref_fc = numberOrError(parameters.posref_fc);
    }
    if (parameters.elastic === "true") {
      runner.elastic = true;
    }
    if (parameters.ef) {
      runner.ef = numberOrError(parameters.ef);
    }
    if (parameters.el) {
      runner.el = numberOrError(parameters.el);
    }
    if (parameters.eu) {
      runner.eu = numberOrError(parameters.eu);
    }
    if (parameters.ea) {
      runner.ea = numberOrError(parameters.ea);
    }
    if (parameters.ep) {
      runner.ep = numberOrError(parameters.ep);
    }
    if (parameters.em) {
      runner.em = numberOrError(parameters.em);
    }
    if (parameters.eb && typeof parameters.eb === 'string') {
      const arg = shellescape([parameters.eb]);
      runner.eb = arg.split(',');
    }
    if (parameters.use_go === "true") {
      runner.use_go_virtual_sites = true;
    }
    if (parameters.sc_fix === "true") {
      runner.side_chain_fix = true;
    }
    if (parameters.cter !== '') {
      runner.cter = parameters.cter
    }
    if (parameters.nter !== '') {
      runner.nter = parameters.nter
    }
    if (parameters.neutral_termini === "true") {
      runner.neutral_termini = true
    }
    if (parameters.cystein_bridge) {
      runner.cystein_bridge = parameters.cystein_bridge
    }
    if (parameters.commandline !== undefined) {
      runner.commandline = parameters.commandline
    }

   return runner;
}

/*async function handleHistory(user_id:string, job_id: string){
  logger.debug("handleHistory")
  const exists = await Database.history.exists(user_id)
  logger.debug(`${exists}`)
  if(!exists){
    logger.debug("Create new user history")
    const doc: History = {
      id: user_id, 
      job_ids: [job_id]
    }
    await Database.history.save(doc)
  }
  else{

  }
  
}*/

async function martinizeRun(parameters: any, pdb_path: string, onStep?: (step: string, ...data: any[]) => void, path?: string) {
  //const { ff, position, posref_fc, elastic, ef, el, eu, ea, ep, em, eb, use_go, sc_fix, nter, cter, neutral_termini, commandline, cystein_bridge } = parameters;

  const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));

  /// NON PRIS EN CHARGE (TODO) ///
  /**
   *  ignore?: string[];
   *  ignh?: boolean;
   *  collagen?: boolean;
   *  dihedral?: boolean;
   *  cystein_bridge?: string;
   */

  const runner = createRunner(settings, parameters, pdb_path);

  try {
    const { pdb, itps, top, warns, dir } = await Martinizer.run(runner, onStep, path);

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
      warns,
      dir,
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

export async function SocketIoMartinizer(app: Server) {
  const io = SocketIo(app);

  let dir = await TmpDirHelper.get();
  //console.log(dir);

  const jobOpt: JobInputs = { 
    exportVar: {
      basedir: dir,
      martinizeArgs: "--version",
    },
    inputs: {},
  };

  await ShellManager.run(
    'martinize_version', 
    ShellManager.mode === "jm" ? jobOpt : "",  
    dir, 
    'martinize2'
  );
  let version = await FsPromise.readFile(dir+"/martinize_version.stdout", 'utf-8');


  io.on('connection', socket => {
    socket.on('previewMartinize', async (settings: any) => {

      const settings_file: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));
      const runner = createRunner(settings_file, settings);

      let {command_line} = Martinizer.settingsToCommandline(runner)
      
      socket.emit('martinizePreviewContent', command_line)
    })

    socket.emit('martinizeVersion', version);

    socket.on('martinize', async (file: Buffer, run_id: string, settings: any, userId:string) => {
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
      
      logger.info(`MARTINIZE USER ${userId}`); 
      try {
        await FsPromise.writeFile(INPUT, file);

        const { pdb, itps, top, elastic_bonds, warns, dir } = await martinizeRun(
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

        await sendFile(warns, { 
          name: path.basename(warns),
          type: 'martinize-warnings' 
        });

        socket.emit('martinize before end', { id: run_id });

        const radius = await Database.radius.getRadius(
          settings.ff || 'martini22',
          itps,
        );

        let stdout : string[] = [];
        await FsPromise.readFile(dir + '/martinize.stderr', 'utf-8')
          .then(function(result) {
            let tmp = result.split('\n');
            tmp.forEach(line => {
              if(line.match('WARNING')) {
                stdout.push(line);
              }
            });
          })
          .catch(function(error) {
            console.log("ERROR: " + error);
          });
        socket.emit('martinize stderr', stdout);

        socket.emit('martinize end', { id: run_id, elastic_bonds, radius });

        logger.info(`MARTINIZE END ${run_id}`)

        logger.info(`DIRECTORY ${dir}`)
        
        //ADD TO HISTORY
        logger.info("ADD TO HISTORY")
        
        const jobId = path.basename(dir); 

        HistoryOrganizer.save(jobId, [top, pdb, ...itps])
          .then(() => {
            Database.history.addToHistory(userId, jobId).then(ok => console.log("saved to db")).catch(e => console.log(e))
          })
          .catch(err => console.log(err))
        
        


        

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
