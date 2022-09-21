import { Router } from 'express';
import { methodNotAllowed, cleanMulterFiles, errorCatcher, generateSnowflake, validateToken, dateFormatter } from '../../helpers';
import Uploader from '../Uploader';
import { Martinizer, MartinizeSettings, ElasticOrGoBounds, GoMoleculeDetails } from '../../Builders/Martinizer';
import { SETTINGS_FILE, MARTINIZE_VERSION, URLS, SEND_COMPLETION_MAIL } from '../../constants';
import { SettingsJson } from '../../types';
import { promises as FsPromise } from 'fs';
import Errors, { ErrorType, ApiError } from '../../Errors';
import shellescape from 'shell-escape';
import logger from '../../logger';
import path from 'path';
import { Database } from '../../Entities/CouchHelper';
import { Job } from '../../Entities/entities';
import SocketIo from 'socket.io';
import TmpDirHelper from '../../TmpDirHelper';
import { Server } from 'http';
import ShellManager, { JobInputs } from '../../Builders/ShellManager';
import HistoryOrganizer from "../../HistoryOrganizer";
import Mailer from '../../Mailer/Mailer';
import { TopFile } from 'itp-parser-forked'
import { plainToInstance } from 'class-transformer';
import { ClientSettingsMartinize } from './molecule.dto';
import { validateOrReject } from 'class-validator';
import { MartinizeJobToSave } from './molecule.types'

type MartinizeRunFailedPayload = {
  error: string,
  type: string,
  stdout?: string,
  stderr?: string,
  dir: string,
  code: ErrorType,
  message: string,
};

interface ClientSettings {
  ff: string,
  position: string;
  cter: string
  nter: string
  sc_fix: string;
  cystein_bridge: string;
  elastic?: string;
  ef?: string;
  el?: string;
  eu?: string;
  ea?: string;
  ep?: string;
  em?: string;
  use_go?: string;
  builder_mode: string;
  send_mail: string;
  user_id?: string;
  pdb_name?: string;
}

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
  if (isNaN(pos) || pos < 0) {
    return Errors.throw(ErrorType.Format);
  }
  return pos;
}

function createRunner(settings: any, parameters: any, pdb_path?: string) {
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

    if (parameters.use_go === "true") {
      runner.use_go = true;
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
    if(parameters.builder_mode !== undefined){
      runner.builder_mode = parameters.builder_mode; 
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

async function martinizeRun(parameters: ClientSettingsMartinize, pdb_path: string, onStep?: (step: string, ...data: any[]) => void, path?: string) {
  //const { ff, position, posref_fc, elastic, ef, el, eu, ea, ep, em, eb, use_go, sc_fix, nter, cter, neutral_termini, commandline, cystein_bridge } = parameters;

  /// NON PRIS EN CHARGE (TODO) ///
  /**
   *  ignore?: string[];
   *  ignh?: boolean;
   *  collagen?: boolean;
   *  dihedral?: boolean;
   *  cystein_bridge?: string;
   */

  console.log("martinizeRun")

  const martinizeSettings: MartinizeSettings = Object.assign({}, {
    input: pdb_path,
    ff: 'martini22',
    position: 'none',
    commandline: ''
  }, parameters);

  try {
    const { pdb, itps, top, warns, dir, elastic_bonds } = await Martinizer.run(martinizeSettings, onStep, path);

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

async function sendMailMartinizeEnd(userId: string, jobId: string) {
  const user = await Database.user.get(userId);
  logger.debug(`Send an email to ${user.email} for job completion`)
  Mailer.send({
    to: user.email,
    subject: "MArtini Database - Job completed"
  },
    "mail_job_completed", {
    name: user.name,
    job_id: jobId,
    job_url: URLS.SERVER + '/builder/' + jobId
  }).catch(logger.error)
}

export async function SocketIoMartinizer(socket: SocketIo.Socket) {
  socket.on('martinize', async (file: Buffer, run_id: string, settings: ClientSettings) => {
    function sendFile(path: string, infos: { id?: string, name: string, type: string, mol_idx?: number }) {
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
      }) as Promise<void>;
    }




    if (!run_id || !file || !settings || !settings.user_id) {
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

    //SECURITY : CHECK PDB CONTENT BEFORE WRITE

    const tmp_dir = await TmpDirHelper.get();
    const INPUT = `${tmp_dir}/input.pdb`
    try {
      const validatedParams = plainToInstance(ClientSettingsMartinize, settings)
      try {
        await validateOrReject(validatedParams)
      } catch (e) {
        throw new Error(e as any);
      }

      logger.debug(`[MARTINIZE] save input to ${tmp_dir}`)


      await FsPromise.writeFile(INPUT, file);

      const { pdb, itps, top, elastic_bonds, warns, dir } = await martinizeRun(
        validatedParams,
        INPUT,
        (step, ...data) => {
          socket.emit('martinize step', {
            id: run_id,
            step,
            data,
          });
        },
        tmp_dir
      );

      await sendFile(top, {
        name: path.basename(top),
        type: 'chemical/x-topology'
      });

      await sendFile(pdb, {
        name: path.basename(pdb),
        type: 'chemical/x-pdb'
      });

      for (const [mol, itp_files] of itps.entries()) {
        for (const itp of itp_files) {
          await sendFile(itp, {
            name: path.basename(itp),
            type: 'chemical/x-include-topology',
            mol_idx: mol
          });
        }

      }

      await sendFile(warns, {
        name: path.basename(warns),
        type: 'martinize-warnings'
      });

      socket.emit('martinize before end', { id: run_id });
      const flatItps = itps.flat()
      const radius = await Database.radius.getRadius(
        validatedParams.ff || 'martini22',
        flatItps
      );


      const job: MartinizeJobToSave = {
        jobId: path.basename(dir),
        userId: validatedParams.user_id,
        type: "martinize",
        date: dateFormatter("Y-m-d H:i"),
        files: {
          all_atom: path.basename(INPUT),
          coarse_grained: path.basename(pdb),
          itp_files: itps.map(mol_itps => mol_itps.map(itp => path.basename(itp))),
          top_file: path.basename(top),
          warnings: path.basename(warns)
        },
        settings: validatedParams, //To avoid class into class that create conflicts for plainToInstance later (why ??)
        radius,
        name: validatedParams.pdb_name
      }

      let savedToHistory = false;
      try {
        await HistoryOrganizer.saveToHistory(job, [INPUT, top, pdb, ...itps.flat(), warns])
        savedToHistory = true;
      } catch (e) {
        logger.warn("error save to history", e)
      }

      finally {

        socket.emit('martinize end', { id: run_id, elastic_bonds, radius, savedToHistory, jobId: job.jobId });

        if (validatedParams.send_mail && job.userId) sendMailMartinizeEnd(job.userId, job.jobId);
      }


    } catch (e) {
      // Error catch, test the error :D
      console.error("ERROR", e)
      if (e instanceof ApiError) {
        if (e.code === ErrorType.MartinizeRunFailed) {
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
          const { error } = e.data
          socket.emit('martinize error', {
            id: run_id,
            error,
            stack: e.stack
          })
        }
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
}

MartinizerRouter.post('/', Uploader.single('pdb'), (req, res) => {
  const pdb_file = req.file;
  if (!pdb_file) {
    return Errors.throw(ErrorType.MissingParameters);
  }

  (async () => {
    const { pdb, itps, top, elastic_bonds } = await martinizeRun(req.body, pdb_file.path);

    // Formatting itps

    const res_itp: { content: string, name: string, type: string, }[][] = [];
    for (const mol_itps of itps) {

      const readed_itps: { content: string, name: string, type: string, }[] = []
      for (const itp of mol_itps) {
        readed_itps.push({
          content: await FsPromise.readFile(itp, 'utf-8'),
          name: path.basename(itp),
          type: 'chemical/x-include-topology',
        });
      }
      res_itp.push(readed_itps)
    }

    // todo: create the custom pdb? à voir

    // Get the radius for itps
    const radius = await Database.radius.getRadius(
      req.body.ff || 'martini22',
      itps.flat(),
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

MartinizerRouter.get('/version', (req, res) => {
  res.json({ version: MARTINIZE_VERSION })
})

MartinizerRouter.all('/version', methodNotAllowed('GET'))
MartinizerRouter.all('/', methodNotAllowed('POST'));

export default MartinizerRouter;
