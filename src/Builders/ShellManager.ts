import { JOB_MANAGER_SETTINGS, INSANE_PATH, INSANE_PATH_JM, POLYPLY_PATH, POLYPLY_PATH_JM, CONECT_PDB_PATH, CONECT_PDB_PATH_JM, CREATE_MAP_PATH, CREATE_MAP_PATH_JM, CREATE_GO_PATH, CREATE_GO_PATH_JM, MARTINIZE_PATH, MARTINIZE_PATH_JM, JobMethod, DEFAULT_JOB_METHOD, SLURM_PROFILES, MARTINIZE_VENV, INSANE_VENV, POLYPLY_VENV } from '../constants';
import { exec } from 'child_process';
import fs from 'fs';
import { ArrayValues } from '../helpers';
// @ts-ignore
import jmClient from 'ms-jobmanager'
import { inspect } from 'util';
import logger from '../logger';
import { Readable, Stream } from 'stream';
import { EventEmitter } from 'events';
import { dir } from 'console';
import { rejects } from 'assert';

const SupportedScripts = ['insane', 'conect', 'go_virt', 'ccmap', 'martinize', 'polyply'] as const;
export type SupportedScript = ArrayValues<typeof SupportedScripts>;

export interface JobInputs {
  exportVar?: { [key: string]: string },
  inputs: { [key: string]: any }
  modules?: string[]
};

export default new class ShellManager {
  /**
   * Which mode to use when tasks are started.
   */
  public mode: JobMethod = DEFAULT_JOB_METHOD;
  private _jm?: Promise<void>;
  private connect: boolean = false;
  /**
   * Link a script name `SupportedScript` to a .sh path.
   */
  protected readonly NAME_TO_PATH: { [scriptName in SupportedScript]: string } = {
    'conect': CONECT_PDB_PATH,
    'go_virt': CREATE_GO_PATH,
    'ccmap': CREATE_MAP_PATH,
    'insane': INSANE_PATH,
    'martinize': MARTINIZE_PATH,
    'polyply': POLYPLY_PATH
  };

  /**
   * Assign the following variables into child processes env.
   */
  protected readonly VARIABLES_TO_NAME: { [scriptName in SupportedScript]: object } = {
    'conect': {},
    'go_virt': {
      venv: MARTINIZE_VENV
    },
    'ccmap': {
      venv: MARTINIZE_VENV
    },
    'insane': {
      venv: INSANE_VENV
    },
    'martinize': {
      venv: MARTINIZE_VENV
    },
    'polyply': {
      venv: POLYPLY_VENV
    }
  };

  /**
   * Script name to jobOpt ? Specify here parameters to fill in job opt ?
   */
  private engine = {
    "engineSpecs": 'slurm',
    "binariesSpec": {
      "submitBin": "/usr/bin/sbatch",
      "cancelBin": "/usr/bin/scancel",
      "queueBin": "/usr/bin/squeue"
    }
  }

  protected readonly NAME_TO_ARGS: { [scriptName in SupportedScript]: any } = {
    'conect': {
      'script': CONECT_PDB_PATH_JM,
      'modules': ['gromacs'],
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    },
    'go_virt': {
      'script': CREATE_GO_PATH_JM,
      'modules': ['mad-utils'],
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    },
    'ccmap': {
      'script': CREATE_MAP_PATH_JM,
      'modules': ['mad-utils'],
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    },
    'insane': {
      'script': INSANE_PATH_JM,
      'modules': ['insane'],
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    },
    'martinize': {
      'script': MARTINIZE_PATH_JM,
      'modules': ['martinize2'],
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    },
    'polyply': {
      'script': POLYPLY_PATH_JM,
      'jobProfile': SLURM_PROFILES.JOB_PROFILE,
      'sysSettingsKey': SLURM_PROFILES.SYS_SETTINGS
    }
  };

  /**
   * Run a given script {script_name} with args {args} in {working_directory}, and save stdout/stderr to {save_std_name}.std<type>.
   */
  async run(script_name: SupportedScript, args: string | JobInputs, working_directory: string, save_std_name?: string | false, timeout?: number, mode: JobMethod = this.mode): Promise<void | any> {
    return new Promise(async (res, rej) => {
      try {
        if (mode === 'jm') {
          const myJob = await this.runWithJobManager(script_name, args as JobInputs, working_directory, save_std_name, timeout);
          res(myJob)
        } else {
          const vide = this.runWithChildProcess(script_name, args as string, working_directory, save_std_name, timeout);
          res(vide)
        }
      }
      catch (e) {
        rej(e)
      }
    }
    )
  }

  protected runWithChildProcess(script_name: SupportedScript, args: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const path = this.NAME_TO_PATH[script_name];
    const variables = this.VARIABLES_TO_NAME[script_name];
    if (!path) {
      throw new Error("Script is not supported.");
    }

    return new Promise((resolve, reject) => {
      let stdout: fs.WriteStream | undefined = undefined;
      let stderr: fs.WriteStream | undefined = undefined;

      if (save_std_name) {
        stdout = fs.createWriteStream(working_directory + '/' + save_std_name + '.stdout');
        stderr = fs.createWriteStream(working_directory + '/' + save_std_name + '.stderr');
      }
      logger.silly(`run command line : "${path}" ${args}`);
      const child = exec(`"${path}" ${args}`, {
        cwd: working_directory,
        maxBuffer: 1e9,
        timeout,
        env: Object.assign({}, process.env, variables)
      }, (err) => {
        stdout?.close();
        stderr?.close();
        child.stderr?.removeAllListeners();

        //logger.error(`${err}`)
        if (err) {
          reject({
            error: err,
            stdout: stdout ? (working_directory + '/' + save_std_name + '.stdout') : undefined,
            stderr: stderr ? (working_directory + '/' + save_std_name + '.stderr') : undefined,
          });
          return;
        }
        resolve();
      });

      if (stdout && stderr) {
        child.stdout?.pipe(stdout);
        child.stderr?.pipe(stderr);
      }
    }) as Promise<void>;
  }

  /**
   * Run a given script {script_name} with args {args} with the job manager inside {working_directory}.
   * 
   * {working_directory} may contain files needed by the job manager.
   * Save stdout/stderr (if {save_std_name} is string) to {working_directory}/{save_std_name}.<stdout/stderr>.
   * Autokill job (if {timeout} is defined) after {timeout} milliseconds.
   * 
   * If job create -new files- in its directory, they *must* be copied to {working_directory} after the job.
   */
  protected async runWithJobManager(script_name: SupportedScript, jobData: JobInputs, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const options = this.NAME_TO_ARGS[script_name];
    const jobOpt = { ...options, ...jobData };
    if (save_std_name) jobOpt.exportVar["OUTPUT_PREFIX"] = save_std_name

    logger.debug("JM PROCESS");

    if (!options) {
      throw new Error("Script is not known, unable to load job manager settings");
    }

    const dumpFile = (srcStream: any, targetPath: any) => new Promise((resolve, reject) => {
      const targetStream = fs.createWriteStream(targetPath);
      targetStream.on('finish', resolve);
      targetStream.on('error', reject);
      srcStream.pipe(targetStream);
    });

    logger.silly("Getting Job manager connection...");

    logger.silly(`Passing following job to ms-jobmanager: ${inspect(jobOpt)}`);

    try {
      jmClient.start(JOB_MANAGER_SETTINGS.address, JOB_MANAGER_SETTINGS.port)
    }
    catch (e) {
      throw new JMError(`Error with job manager : ${e}`)
    }

    // return new Promise(async (resolve, reject) => {
    try {
      const { stdout, jobFS } = await jmClient.pushFS(jobOpt)
      console.log(working_directory + '/' + save_std_name)


      // ##############################################
      // A retablir absolument pour les autres scripts le border de dumpfile pour ne pas causer de bug 

      // Chercher une regex dans jobFS list 
      //       '4dc2a08e-5600-4a4d-b767-39841693269e_coreScript.sh',
      // [1]   '4dc2a08e-5600-4a4d-b767-39841693269e.batch',
      // [1]   '4dc2a08e-5600-4a4d-b767-39841693269e.err',
      // [1]   '4dc2a08e-5600-4a4d-b767-39841693269e.out',

      //stdout, working_directory + '/' + save_std_name + '.stdout');
      const jobfilelist : string[] = await jobFS.list()

      let idJM = jobfilelist.map( (x) => {if(x.endsWith("_coreScript.sh")) return x} )[0]?.replace("_coreScript.sh","")
      console.log( "#JOBFS#",idJM)
      // Grosse regex pour retrouver 
      const stdout_fname = idJM+".out" // grace à la recherche dans list
      const stream_fname = idJM+".err"
      const stream_stdout: Readable = await jobFS.readToStream(stdout_fname);
      const stream_stderr: Readable = await jobFS.readToStream(stream_fname);
      // await dumpFile(stream_stdout, .....)
      await dumpFile(stream_stdout, working_directory + '/' + save_std_name + '.stdout');
      await dumpFile(stream_stderr, working_directory + '/' + save_std_name + '.stderr');
      return { stdout, jobFS }
    }
    catch (e) {
      console.log(e)
      throw (e)
    }

    // jobCreatePdbWithConect.on('completed', (stdout: Stream, stderr: Stream) => {
    //   logger.debug(`jobCreatePdbWithConect completed`);
    //   if (save_std_name) {
    //     (async () => {
    //       await dumpFile(stdout, working_directory + '/' + save_std_name + '.stdout');
    //       await dumpFile(stderr, working_directory + '/' + save_std_name + '.stderr');
    //     })()
    //       .then(resolve)
    //       .catch(reject);

    //     return;
    //   }
    //   resolve();
    // });

    // jobCreatePdbWithConect.on('error', reject);
    // jobCreatePdbWithConect.on("disconnect_error", () => {
    //   reject(new JMError(`Error with job manager : Job manager has been disconnected`))
    // })
    // jobCreatePdbWithConect.on("lostJob", () => {
    //   reject(new JMError(`Error with job manager : Job has been lost`))
    // })
  }


  // ) as Promise<void>;

}
// Job manager.start now return a promise, this function doesn't seem to be use anywhere 
// get job_manager() {
//   if (this._jm) {
//     return this._jm;
//   }
//   return this._jm = JobManager.start(JOB_MANAGER_SETTINGS.address,   JOB_MANAGER_SETTINGS.port)
// }
//}();

export class JMError extends Error {
}
