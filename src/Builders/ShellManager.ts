import { ENABLE_JOB_MANAGER } from '../constants';
import { exec } from 'child_process';
import fs from 'fs';


export default new class ShellManager {
  /**
   * Run a given shell {command} in {working_directory}, and save stdout/stderr to {save_std_name}.std<type>.
   */
  async run(command: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    if (ENABLE_JOB_MANAGER) {
      return this.runWithJobManager(command, working_directory, save_std_name, timeout);
    }
    return this.runWithChildProcess(command, working_directory, save_std_name, timeout);
  }

  protected runWithChildProcess(command: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    return new Promise((resolve, reject) => {
      let stdout: fs.WriteStream | undefined = undefined;
      let stderr: fs.WriteStream | undefined = undefined;

      if (save_std_name) {
        stdout = fs.createWriteStream(working_directory + '/' + save_std_name + '.stdout');
        stderr = fs.createWriteStream(working_directory + '/' + save_std_name + '.stderr');
      }

      const child = exec(command, { cwd: working_directory, maxBuffer: 1e9, timeout }, (err) => {
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
   * Run a shell command {command} with the job manager inside {working_directory}.
   * 
   * {working_directory} may contain files needed by the job manager.
   * Save stdout/stderr (if {save_std_name} is string) to {working_directory}/{save_std_name}.<stdout/stderr>.
   * Autokill job (if {timeout} is defined) after {timeout} milliseconds.
   * 
   * If job create -new files- in its directory, they *must* be copied to {working_directory} after the job.
   */
  protected async runWithJobManager(command: string, working_directory: string, save_std_name?: string | false, timeout?: number) {
    // todo
  }
}();
