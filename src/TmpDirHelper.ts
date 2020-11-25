import fs, { promises as FsPromise } from 'fs';
import os from 'os';
import logger from './logger';
import { simpleflake } from 'simpleflakes';
import { DEFAULT_TMP_BASE_DIR } from './constants';

export type TmpDirMode = 'os' | 'directory';

/**
 * Singleton that helps managing temporary directories.
 */
export const TmpDirHelper = new class TmpDirHelper {
  protected cache: [string, number][] = [];

  public mode: TmpDirMode = 'os';

  constructor() {
    
  }

  async program_clean(){
    logger.debug("Program cache cleaning")
    setInterval(() => {

      // Remove every item aged more than now - 45 minutes
        this.clean(
          Date.now() - (1000 * 60 * 2)
        );
      }, 1000 * 60 * 2);
  }

  /**
   * Get a new temporary directory. Returns its path, without trailing /.
   */
  async get() {
    let dir: string;

    if (this.mode === 'os') {
      dir = await this.getRandomTmpDirFromOs();
    }
    else {
      dir = await this.getRandomTmpDirFromBaseDirectory();
    }

    this.cache.push([dir, Date.now()]);

    return dir;
  }

  protected async getRandomTmpDirFromOs() {
    const tmp_dir = os.tmpdir();
    const dir = await FsPromise.mkdtemp(tmp_dir + "/");

    return dir;
  }

  protected async getRandomTmpDirFromBaseDirectory() {
    const base = DEFAULT_TMP_BASE_DIR;
    const suffix = simpleflake().toString();

    const dir = base + suffix;

    await FsPromise.mkdir(dir, { recursive: true });
    // Mode does not work with mkdir, must do the chmod
    await FsPromise.chmod(dir, 0o777);

    return dir
  }

  /**
   * Revoke a specific directory.
   */
  revoke(dir: string) {
    const to_remove = this.cache.filter(e => e[0].startsWith(dir));
    this.cache = this.cache.filter(e => !e[0].startsWith(dir));

    return this.removeDirectories(to_remove.map(e => e[0]));
  }

  /**
   * Remove directories that has been created before {max_timestamp}.
   * If {max_timestamp} is not specified, it cleans everything.
   */
  clean(max_timestamp = Infinity) {
    const to_remove = this.cache.filter(e => e[1] < max_timestamp).map(e => e[0]);

    if (to_remove.length === 0) {
      return Promise.resolve();
    }

    this.cache = this.cache.filter(e => e[1] >= max_timestamp);

    logger.debug("Removing directories: " + to_remove.join(', '));

    return this.removeDirectories(to_remove);
  }

  protected removeDirectories(dir_entries: string[]) {
    for (const e of dir_entries) {
      // recursive does not work with non-sync method
      try {
        fs.rmdirSync(e, { recursive: true });
      } catch {
        logger.error('[TmpDirHelper] Unable to erase temporary directory. (' + e + ')');
      }
    }
  }
};

export default TmpDirHelper;
