import CliHelper, { CliListener } from "interactive-cli-helper";
import SearchWorker from "../search_worker";
import { MINUTES_BEFORE_WORKER_KILL } from "../constants";

const WORKER_CLI = new CliListener(
  CliHelper.formatHelp("worker", {
    spawn: 'Spawn a new worker',
    'assign timeout <id>': 'Set a kill timeout for worker <id>',
    list: 'List started workers',
    'get <id>': 'Get details about worker <id>',
    'kill <id>/all': 'Kill worker <id> / Kill all workers',
  })
);

WORKER_CLI.addSubListener('spawn', () => {
  const id = SearchWorker.spawn(false);

  return `A new worker has been spawned with ID ${id}.`;
});

WORKER_CLI.addSubListener('assign timeout', rest => {
  const id = rest.trim();

  if (!id) {
    return `You must specify a worker ID.`;
  }

  if (!SearchWorker.exists(id)) {
    return `Worker ${id} does not exists.`;
  }

  SearchWorker.assignKillTimeout(id);

  return `Worker ${id} will be killed in ${MINUTES_BEFORE_WORKER_KILL} minutes.`;
});

WORKER_CLI.addSubListener('list', () => {
  const workers = SearchWorker.available;

  if (Object.keys(workers).length === 0) {
    return "There aren't any started worker.";
  }

  return `Available workers are \n- ${
    Object.entries(workers)
      .map(w => `Worker ${w[0]} (occupation ${w[1][0]})`)
      .join('\n- ')
  }`;
});

WORKER_CLI.addSubListener('get', rest => {
  const workers = SearchWorker.available;
  const id = rest.trim();

  if (!id) {
    return `Please specify a worker id.`;
  }

  if (!(id in workers)) {
    return "This worker does not exists.";
  }

  const w = workers[id];
  return `Worker ${id}: Occupation ${w[0]}`;
});

WORKER_CLI.addSubListener('kill', rest => {
  const id = rest.trim();
  
  if (!id) {
    return `Please specify a worker id or "all".`;
  }

  if (rest === "all") {
    const available = SearchWorker.available;

    for (const w_id of Object.keys(available)) {
      SearchWorker.kill(w_id);
    }

    return `All child workers has been killed.`;
  }

  if (SearchWorker.exists(id)) {
    SearchWorker.kill(id);
    return `Child worker ${id} has been killed.`;
  }
  return `Unable to find specified worker.`
});

export default WORKER_CLI;

