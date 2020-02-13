import { parentPort, workerData } from 'worker_threads';
import { Molecule } from '../Entities/entities';
import nano = require('nano');
import md5 from 'md5';

interface WorkerIncomingMessage {
  type: "search_end" | "error";
  uuid: string;
  molecules?: Molecule[];
  length?: number;
}

interface WorkerSendMessage {
  type: "new_search" | "clean";
}

interface WorkerStartTask extends WorkerSendMessage {
  uuid: string;
  query: nano.MangoQuery;
  as_all?: boolean;
}

interface MoleculeResultTree {
  [tree_id: string]: Molecule[];
}

// Read from dat invoked by parent
const WORKER_ID: string = workerData.id;
const MOLECULE_COLLECTION: string = workerData.molecule_collection;
const COUCH_URL: string = workerData.couch_url;

// Create couch connection
const CONNECTION = nano({ url: COUCH_URL, requestDefaults: { proxy: null } });
const MOLECULE_DATABASE = CONNECTION.use<Molecule>(MOLECULE_COLLECTION);
const QUERY_LIMIT = 9999999;

// Create cache
let Cache: { [hash: string]: MoleculeResultTree } = {};

setInterval(() => {
  // Keep the 10 most recent cached
  const keys = Object.keys(Cache);
  const end = keys.length > 10 ? keys.length - 10 : 0;

  // Keys are always ordered in insert order !
  for (const key of keys.slice(0, end)) {
    delete Cache[key];
  }
}, 2 * 60 * 1000);

// -----------------------
// - LISTEN FOR MESSAGES -
// -----------------------
parentPort!.on('message', (msg: WorkerSendMessage) => {
  if (msg.type === "new_search") {
    startNewSearch(msg as WorkerStartTask);
  }
  if (msg.type === "clean") {
    cleanCache();
  }
});


// Worker functions
function cleanCache() {
  Cache = {};
}

async function startNewSearch(msg: WorkerStartTask) {
  const task_id = msg.uuid;
  const query_hash = md5(JSON.stringify(msg.query));

  // Do the search
  let molecule_tree: MoleculeResultTree;
  if (query_hash in Cache) {
    molecule_tree = Cache[query_hash];
  }
  else {
    const query = { ...msg.query };
    query.limit = QUERY_LIMIT;
    // This will be applied later
    query.skip = 0;

    const molecules = await MOLECULE_DATABASE.find(query)
      .catch(() => {
        sendErrorMessage(task_id);
      });

    if (!molecules) {
      return;
    }

    // Create the trees
    molecule_tree = {};
    for (const mol of molecules.docs) {
      if (mol.tree_id in molecule_tree) {
        molecule_tree[mol.tree_id].push(mol);
      }
      else {
        molecule_tree[mol.tree_id] = [mol];
      }
    }

    // Sort every molecule in created at date, the lastest the first
    for (const tree_id in molecule_tree) {
      molecule_tree[tree_id] = molecule_tree[tree_id].sort((a, b) => +new Date(b.created_at) - +new Date(a.created_at));
    }
  }

  // Molecule tree is okay
  const skip = msg.query.skip || 0;
  const limit = msg.query.limit || 25;

  Cache[query_hash] = molecule_tree;

  if (msg.as_all) {
    // Fusion all tree an returns the .slice(skip, skip+limit)
    // This requires Node 11 !
    const all_mols = Object.values(molecule_tree).flat();
  
    sendEndMessage(task_id, all_mols.slice(skip, skip + limit), all_mols.length);
    return;
  }

  // Return tree skip, skip+limit
  const trees = Object.values(molecule_tree);

  const concerned_molecules = trees
    .slice(skip, skip + limit)
    .map(e => e[0]); // Keep the most recent
  
  sendEndMessage(task_id, concerned_molecules, trees.length);
}


// Utilities
function sendEndMessage(task_id: string, molecules: Molecule[], full_length: number) {
  parentPort!.postMessage({
    type: "search_end",
    uuid: task_id,
    molecules,
    length: full_length,
  } as WorkerIncomingMessage);
}

function sendErrorMessage(task_id: string) {
  parentPort!.postMessage({
    type: "error",
    uuid: task_id,
  } as WorkerIncomingMessage);
}
