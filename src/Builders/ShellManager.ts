import { INSANE_PATH, CONECT_PDB_PATH, CREATE_MAP_PATH, CREATE_GO_PATH, MARTINIZE_PATH, JobMethod, DEFAULT_JOB_METHOD } from '../constants';
import { exec } from 'child_process';
import fs from 'fs';
import { ArrayValues } from '../helpers';

const SupportedScripts = ['insane', 'conect', 'go_virt', 'ccmap', 'martinize'] as const;
export type SupportedScript = ArrayValues<typeof SupportedScripts>;


export default new class ShellManager {
  /**
   * Which mode to use when tasks are started.
   */
  public mode: JobMethod = DEFAULT_JOB_METHOD;

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

  protected readonly VARIABLES_TO_NAME: { [scriptName in SupportedScript]: {} } = {
    'conect': {},
    'go_virt': {},
    'ccmap': {},
    'insane': {},
    'martinize': {},
  };

  /**
   * Script name to jobOpt ? Specify here parameters to fill in job opt ?
   */
  protected readonly NAME_TO_ARGS: { [scriptName in SupportedScript]: any } = {
    'conect': {},
    'go_virt': {},
    'ccmap': {},
    'insane': {},
    'martinize': {},
  };

  /**
   * Run a given script {script_name} with args {args} in {working_directory}, and save stdout/stderr to {save_std_name}.std<type>.
   */
  async run(
    script_name: SupportedScript, 
    args: string, 
    working_directory: string,
    save_std_name?: string | false, 
    timeout?: number, 
    mode: JobMethod = this.mode
  ) {
    if (mode === 'jm') {
      return this.runWithJobManager(script_name, args, working_directory, save_std_name, timeout);
    }
    return this.runWithChildProcess(script_name, args, working_directory, save_std_name, timeout);
  }

  protected runWithChildProcess(script_name: SupportedScript, args: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const path = this.NAME_TO_PATH[script_name];
    const variables = this.VARIABLES_TO_NAME[script_name];

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

      const child = exec(`"${path}" ${args}`, { 
        cwd: working_directory, 
        maxBuffer: 1e9, 
        timeout, 
        env: Object.assign({}, process.env, variables) 
      }, (err) => {
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
  protected async runWithJobManager(script_name: SupportedScript, args: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    const options = this.NAME_TO_ARGS[script_name];

    if (!options) {
      throw new Error("Script is not known, unable to load job manager settings");
    }

    // todo
  }
}();
