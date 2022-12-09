import WorkerThreadManager, { WorkerPool } from 'worker-thread-manager';
import nano = require('nano');
import { Molecule } from '../Entities/entities';
import CouchHelper from '../Entities/CouchHelper';
import { MINUTES_BEFORE_WORKER_KILL, MAX_POOL_SIZE, MAX_REQUEST_PER_WORKER_THRESHOLD, URLS } from '../constants';

interface WorkerSendMessage {
  type: "new_search" | "clean";
}

interface WorkerStartTask extends WorkerSendMessage {
  query: nano.MangoQuery;
  as_all?: boolean;
  stashed?: boolean;
}

interface WorkerResult {
  molecules: Molecule[];
  length: number;
}

export const SearchWorker = new class SearchWorker {
  protected _pool?: WorkerPool<WorkerStartTask, WorkerResult>;

  protected get pool() {
    if (!this._pool) {
      this._pool = WorkerThreadManager.spawn<WorkerStartTask, WorkerResult>(
        __dirname + '/worker.js',
        {
          spawnerThreshold: MAX_REQUEST_PER_WORKER_THRESHOLD,
          poolLength: MAX_POOL_SIZE,
          // 2 minutes life max without task
          stopOnNoTask: MINUTES_BEFORE_WORKER_KILL * 60 * 1000,
          workerData: {
            molecule_collection: CouchHelper.MOLECULE_COLLECTION,
            molecule_stashed_collection : CouchHelper.STASHED_MOLECULE_COLLECTION,
            couch_url: URLS.COUCH,
          },
        }
      );
    }

    return this._pool;
  }

  async query(query: nano.MangoQuery, send_all = false, stashed = false) {
    const task: WorkerStartTask = {
      type: 'new_search',
      query,
      as_all: send_all,
      stashed
    };

    return this.pool.run(task);
  }

  clearCache() {
    return this.pool.send({ type: 'clean' } as WorkerStartTask);
  }
}

export default SearchWorker;
