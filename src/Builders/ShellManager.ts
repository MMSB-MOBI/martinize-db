import { JOB_MANAGER_SETTINGS, INSANE_PATH, INSANE_PATH_JM, CONECT_PDB_PATH, CONECT_PDB_PATH_JM, CREATE_MAP_PATH, CREATE_GO_PATH, MARTINIZE_PATH, MARTINIZE_PATH_JM, JobMethod, DEFAULT_JOB_METHOD } from '../constants';
import { exec } from 'child_process';
import fs from 'fs';
import { ArrayValues } from '../helpers';
const jobManager = require('ms-jobmanager');
import { inspect } from 'util';
import logger from '../logger';
import { Stream } from 'stream';

const SupportedScripts = ['insane', 'conect', 'go_virt', 'ccmap', 'martinize'] as const;
export type SupportedScript = ArrayValues<typeof SupportedScripts>;

export interface jobInputs {
    exportVar?:  { [key: string]: string },
    inputs: { [key: string]: any }
};

export default new class ShellManager {
  /**
   * Which mode to use when tasks are started.
   */
  public mode: JobMethod = DEFAULT_JOB_METHOD;
  private _jm?: Promise<void>;
  /**
   * Link a script name `SupportedScript` to a .sh path.
   */
  protected readonly NAME_TO_PATH: { [scriptName in SupportedScript]: string } = {
    'conect': CONECT_PDB_PATH,
    'go_virt': CREATE_GO_PATH,
    'ccmap': CREATE_MAP_PATH,
    'insane': INSANE_PATH,
    'martinize': MARTINIZE_PATH,
  };

  /**
   * Script name to jobOpt ? Specify here parameters to fill in job opt ?
   */
  private engine = { 
    "engineSpecs" : 'slurm',
    "binariesSpec" : { 
      "submitBin" : "/data/www_dev/mad/bin/slurm/bin/sbatch",
      "cancelBin" : "/data/www_dev/mad/bin/slurm/bin/scancel",
      "queueBin"  : "/data/www_dev/mad/bin/slurm/bin/squeue"
    }
  }

  protected readonly NAME_TO_ARGS: { [scriptName in SupportedScript]: any } = {
    'conect': {
      'script' : CONECT_PDB_PATH_JM,
      'modules': ['gromacs'],
      'jobProfile' : "mad-dev",
      'engineOverride' : this.engine
    },
    'go_virt': {},
    'ccmap': {},
    'insane': {
      'script' : INSANE_PATH_JM,
      'modules': ['insane'],
      'jobProfile' : "mad-dev",
      'engineOverride' : this.engine
    },
    'martinize': {
      'script' : MARTINIZE_PATH_JM,
      'modules': ['martinize2'],
      'jobProfile' : "mad-dev",
      'engineOverride' : this.engine
    }
  };

  get job_manager() {
    if (this._jm) return this._jm;
    return this._jm = new Promise(resolve => jobManager.start({ 
      'port': JOB_MANAGER_SETTINGS.port,
      'TCPip': JOB_MANAGER_SETTINGS.address
    }).on('ready', resolve));
  }
  /**
   * Run a given script {script_name} with args {args} in {working_directory}, and save stdout/stderr to {save_std_name}.std<type>.
   */
  async run(script_name: SupportedScript, args: string|jobInputs, working_directory: string, save_std_name?: string | false, timeout?: number) {
    logger.debug(`HERE !!!${this.mode}`);
    if (this.mode === 'jm')
      return this.runWithJobManager(script_name, <jobInputs>args, working_directory, save_std_name, timeout);
    
    return this.runWithChildProcess(script_name, <string>args, working_directory, save_std_name, timeout);
  }

  protected runWithChildProcess(script_name: SupportedScript, args: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const path = this.NAME_TO_PATH[script_name];
    logger.debug("CHILD PROCESS");
    if (!path) {
      throw new Error("Script is not supported.");
    }

    return new Promise((resolve, reject) => {
      let stdout: fs.WriteStream | undefined = undefined;
      let stderr: fs.WriteStream | undefined = undefined;

      if (save_std_name) {
        stdout = fs.createWriteStream(working_directory + '/' + save_std_name + '.stdout');
        stderr = fs.createWriteStream(working_directory + '/' + save_std_name + '.stderr');
      }

      const child = exec(`"${path}" ${args}`, { cwd: working_directory, maxBuffer: 1e9, timeout }, (err) => {
        stdout?.close();
        stderr?.close();
        child.stderr?.removeAllListeners();
  
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
  protected async runWithJobManager(script_name: SupportedScript, jobData: jobInputs, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const options = this.NAME_TO_ARGS[script_name];
    const jobOpt = {...options, ...jobData};

    logger.debug("JM PROCESS");
    
    if (!options) {
      throw new Error("Script is not known, unable to load job manager settings");
    }

    const dumpFile = (srcStream:any, targetPath:any) => new Promise((resolve, reject) => {
      const targetStream = fs.createWriteStream(targetPath);
      targetStream.on("finish", resolve);
      srcStream.pipe(targetStream);
    });

    logger.debug("Getting Job manager connection...");
    await this.job_manager;
    logger.info(`About to run this: ${inspect(jobOpt)}`);
    
    return new Promise( (resolve, reject) => {
        let jobCreatePdbWithConect = jobManager.push(jobOpt);
        jobCreatePdbWithConect.on("completed",(stdout:Stream, stderr:Stream) => {
          logger.info(`jobCreatePdbWithConect completed`);
          if (save_std_name) {
           ( async () => {
              await dumpFile(stdout, working_directory + '/' + save_std_name + '.stdout' );
              await dumpFile(stderr, working_directory + '/' + save_std_name + '.stderr' );             
            } )().then(resolve);
            return;
          }
          resolve();
      });
  }) as Promise<void>;
}
}();
