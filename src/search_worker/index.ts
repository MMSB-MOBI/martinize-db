import nano = require('nano');
import { Molecule } from '../Entities/entities';
import CouchHelper from '../Entities/CouchHelper';
import { MINUTES_BEFORE_WORKER_KILL, MAX_POOL_SIZE, MAX_REQUEST_PER_WORKER_THRESHOLD, URLS } from '../constants';
import WorkerThreadManager, { WorkerPool } from 'worker-thread-manager';

interface WorkerSendMessage {
  type: "new_search" | "clean";
}

interface WorkerStartTask extends WorkerSendMessage {
  query: nano.MangoQuery;
  as_all?: boolean;
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
            couch_url: URLS.COUCH,
          },
        }
      );
    }

    return this._pool;
  }

  async query(query: nano.MangoQuery, send_all = false) {
    const task: WorkerStartTask = {
      type: 'new_search',
      query,
      as_all: send_all
    };

    return this.pool!.run(task);
  }

  clearCache() {
    return this.pool.send({ type: 'clean' } as WorkerStartTask);
  }
}

export default SearchWorker;

/* 
export const SearchWorker = new class SearchWorker {
  protected workers: Map<string, Worker> = new Map;
  protected workers_occupation: Map<string, [number, NodeJS.Timeout?]> = new Map;

  async query(query: nano.MangoQuery, send_all = false) {
    const [worker, worker_occ, worker_id] = this.getAvailableWorker();

    logger.debug(`New search task: worker ${worker_id} has been selected (occupation: ${worker_occ[0]}).`);

    const task_id = simpleflake().toString();

    const mols = await new Promise((resolve, reject) => {
      this.sendMessageToWorker({
        type: "new_search",
        uuid: task_id,
        query,
        as_all: send_all
      } as WorkerStartTask, worker);

      // Ajoute une tâche au worker
      worker_occ[0]++;
  
      const worker_fn = (msg: WorkerIncomingMessage) => {
        if (msg.uuid !== task_id) 
          return;
        
        if (msg.type === "search_end" && msg.molecules) {
          resolve({
            molecules: msg.molecules,
            length: msg.length!,
          });

          // Remove the worker
          worker.off("message", worker_fn);
          // On a fini, lui enlève une tâche
          worker_occ[0]--;

          if (worker_occ[0] <= 0) {
            logger.silly(`Child worker ${worker_id} has no search task left: it will be killed in ${MINUTES_BEFORE_WORKER_KILL} minutes.`);
            // Timeout to kill worker
            this.assignKillTimeout(worker_id);
          }
        }
        if (msg.type === "error") {
          reject(msg.error);
        }
      };
  
      worker.on('message', worker_fn);
    }) as { molecules: Molecule[], length: number };

    return mols;
  } 

  protected getAvailableWorker() : [Worker, [number, NodeJS.Timeout?], string] {
    let selected_worker_id: string;
    if (!this.workers.size) {
      // Spawn a new worker
      selected_worker_id = this.spawn();
    }
    else {
      const smaller_worker_occupied = [...this.workers_occupation.entries()].reduce((prev, val) => prev[1][0] < val[1][0] ? prev : val);
      
      // Si il y a moins de MAX_POOL_SIZE workers et si le worker le moins chargé a plus de MAX_REQUEST_PER_WORKER_THRESHOLD requêtes en cours
      if (
        this.workers_occupation.size < MAX_POOL_SIZE && 
        smaller_worker_occupied[1][0] > MAX_REQUEST_PER_WORKER_THRESHOLD
      ) {
        logger.debug("Spawing a new worker due to load balancing.")
        selected_worker_id = this.spawn();
      } 
      else {
        selected_worker_id = smaller_worker_occupied[0];
      }
    }

    const selected_worker_occ = this.workers_occupation.get(selected_worker_id)!;
    const selected_worker = this.workers.get(selected_worker_id)!;

    // Remove the killer timeout
    if (selected_worker_occ[1]) {
      logger.silly(`Cancelling kill for worker ${selected_worker_id}.`);
      clearTimeout(selected_worker_occ[1]);
      selected_worker_occ[1] = undefined;
    }

    return [selected_worker, selected_worker_occ, selected_worker_id];
  }

  spawn(with_log = true) {
    const worker_id = simpleflake().toString();
    const worker = new Worker(__dirname + '/worker.js', {
      workerData: { 
        id: worker_id.toString(), 
        molecule_collection: CouchHelper.MOLECULE_COLLECTION,
        couch_url: URLS.COUCH,
      }
    });

    if (with_log)
      logger.debug(`Worker ${worker_id} has been spawned.`);
    
    this.workers.set(worker_id, worker);
    this.workers_occupation.set(worker_id, [0, undefined]);
    return worker_id;
  }

  assignKillTimeout(worker_id: string) {
    const selected_occ = this.workers_occupation.get(worker_id);

    if (selected_occ) {
      // Remove the previous timeout
      if (selected_occ[1]) {
        clearTimeout(selected_occ[1]);
        selected_occ[1] = undefined;
      }
  
      selected_occ[1] = setTimeout(() => {
        this.kill(worker_id);
        logger.verbose("Worker killed : " + worker_id);
      }, MINUTES_BEFORE_WORKER_KILL * 60 * 1000);
    }
  }

  protected sendMessageToWorkers(message: WorkerSendMessage) {
    for (const worker of this.workers.values()) {
      this.sendMessageToWorker(message, worker);
    }
  }
  
  protected sendMessageToWorker(message: WorkerSendMessage, worker: Worker) {
    worker.postMessage(message);
  }

  clearCache() {
    this.sendMessageToWorkers({ type: "clean" });
  }

  kill(id: string) {
    const worker = this.workers.get(id)!;

    if (worker) {
      const occ = this.workers_occupation.get(id);
      if (occ && occ[1]) {
        clearTimeout(occ[1]);
      }

      this.workers.delete(id);
      this.workers_occupation.delete(id);
  
      worker.terminate();
    }
  }

  get available() {
    return Object.fromEntries(this.workers_occupation.entries());
  }

  exists(id: string) {
    return this.workers.has(id);
  }
}();
 */
