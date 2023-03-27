import { WorkerChild } from 'worker-thread-manager';
import nano = require('nano');
import md5 from 'md5';
import { Molecule } from '../Entities/entities';

/* TYPES */

interface WorkerInitData {
  molecule_collection: string;
  molecule_stashed_collection: string; 
  couch_url: string;
}

interface WorkerResult {
  molecules: Molecule[];
  length: number;
}

interface WorkerStartTask {
  query: nano.MangoQuery;
  as_all?: boolean;
  stashed? : boolean; 
}

interface MoleculeResultTree {
  [tree_id: string]: Molecule[];
}


/* WORKER DATA */

// Init variables
const CONSTANTS = {
  molecule_db: '',
  molecule_stashed_db : '',
  couch: '',
  query_limit: 999999999,
};

function getDb(stashed = false) {
  // Create couch connection
  const db_to_connect = stashed ? CONSTANTS.molecule_stashed_db : CONSTANTS.molecule_db
  return nano({ 
    url: CONSTANTS.couch, 
    requestDefaults: { proxy: null } 
  }).use<Molecule>(db_to_connect);
}

// Create cache
let Cache: { [hash: string]: MoleculeResultTree } = {};

// Keep only the most recent cached results
setInterval(() => {
  // Keep the 10 most recent cached
  const keys = Object.keys(Cache);
  const end = keys.length > 10 ? keys.length - 10 : 0;

  // Keys are always ordered in insert order !
  for (const key of keys.slice(0, end)) {
    delete Cache[key];
  }
}, 2 * 60 * 1000);

/* WORKER LIFECYCLE */

async function onTask(data: WorkerStartTask) {
  const query_hash = md5(JSON.stringify(data.query));

  const db = getDb();
  const db_stashed = data.stashed ? getDb(true) : undefined; 

  // Do the search
  let molecule_tree: MoleculeResultTree;
  if (query_hash in Cache) {
    molecule_tree = Cache[query_hash];
  }
  else {
    const query = { ...data.query };
    query.limit = CONSTANTS.query_limit;
    // This will be applied later
    query.skip = 0;

    const molecules = await db.find(query);
    const molecules_stashed = db_stashed ? await db_stashed.find(query) : undefined

    if (!molecules && !molecules_stashed) {
      return {
        molecules: [],
        length: 0,
      };
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
    if(molecules_stashed){
      for (const mol of molecules_stashed.docs){
        const new_mol = {...mol, from_stashed : true}
        if (new_mol.tree_id in molecule_tree) {
          molecule_tree[new_mol.tree_id].push(new_mol);
        }
        else {
          molecule_tree[new_mol.tree_id] = [new_mol];
        }
      }
    }
   

    // Sort every molecule in created at date, the lastest the first
    for (const tree_id in molecule_tree) {
      molecule_tree[tree_id] = molecule_tree[tree_id].sort((a, b) => +new Date(b.created_at) - +new Date(a.created_at));
    }
  }

  // Molecule tree is okay
  const skip = data.query.skip || 0;
  const limit = data.query.limit || 25;

  Cache[query_hash] = molecule_tree;

  if (data.as_all) {
    // Fusion all tree an returns the .slice(skip, skip+limit)
    // This requires Node 11 !
    const all_mols = Object.values(molecule_tree).flat();
  
    const molecules = all_mols.slice(skip, skip + limit);

    return {
      molecules,
      length: all_mols.length,
    };
  }

  // Return tree skip, skip+limit
  const trees = Object.values(molecule_tree);

  const molecules = trees
    .slice(skip, skip + limit)
    .map(e => e[0]); // Keep the most recent
  
  return {
    molecules,
    length: trees.length,
  };
}

function onStartup(data: WorkerInitData) {
  CONSTANTS.couch = data.couch_url;
  CONSTANTS.molecule_db = data.molecule_collection;
  CONSTANTS.molecule_stashed_db = data.molecule_stashed_collection
}

function onMessage(data: { type: string }) {
  if (data.type !== 'clean')
    return;

  Cache = {};
}

new WorkerChild<WorkerStartTask, WorkerResult, WorkerInitData>({
  onTask,
  onStartup,
  onMessage
}).listen();
