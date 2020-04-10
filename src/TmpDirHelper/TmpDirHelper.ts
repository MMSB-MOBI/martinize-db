import { promises as FsPromise } from 'fs';
import os from 'os';
import logger from '../logger';

/**
 * Singleton that helps managing temporary directories.
 */
export const TmpDirHelper = new class TmpDirHelper {
  protected cache: [string, number][] = [];

  constructor() {
    // register the callback every 30 minutes
    setInterval(() => {

      // Remove every item aged more than now - 45 minutes
      this.clean(
        Date.now() - (1000 * 60 * 45)
      );
    }, 1000 * 60 * 30);
  }

  /**
   * Get a new temporary directory. Returns its path, without trailing /.
   */
  async get() {
    const tmp_dir = os.tmpdir();
    const dir = await FsPromise.mkdtemp(tmp_dir + "/");

    this.cache.push([dir, Date.now()]);

    return dir;
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
    return Promise.all(
      dir_entries.map(e => FsPromise.rmdir(e, { recursive: true }).catch(console.error))
    ).then(() => {});
  }
};

export default TmpDirHelper;
