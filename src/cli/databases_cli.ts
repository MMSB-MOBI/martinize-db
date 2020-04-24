import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import AbstractDatabase from "../Entities/AbstractDatabase";

const DATABASE_CLI = new CliListener(
  CliHelper.formatHelp("database", {
    'create <name>/all': 'Create a single or all databases. Available names are: ' + CouchHelper.DBS.join(', ') + '.',
    'wipe <name>/all': 'Delete a single or all databases.',
    info: 'Check existence of each database and show their document count.',
    'get <database> <docId>': 'Get a document by id in a selected database.',
  }, "Command is incorrect. Type \"database\" for help.")
);

DATABASE_CLI.addSubListener('create', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a database name or "all".`;
  }

  if (rest === "all") {
    await Database.createAll();
    return "All databases has been created.";
  }
  if (!CouchHelper.DBS.includes(rest)) {
    return `This database name is not authorized. Available names are ${CouchHelper.DBS.join(', ')}.`;
  }
  await Database.create(rest);
  return `Database ${rest} has been created.`;
});

DATABASE_CLI.addSubListener('wipe', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a database name or "all".`;
  }

  if (rest === "all") {
    await Database.deleteAll();
    return "All databases has been wiped.";
  }
  if (!CouchHelper.DBS.includes(rest)) {
    return `This database name is not authorized. Available names are ${CouchHelper.DBS.join(', ')}.`;
  }
  await Database.delete(rest);
  return `Database ${rest} has been wiped.`;
});

DATABASE_CLI.addSubListener('get', async rest => {
  const [database, id] = rest.split(/ +/);
  if (!database || !id) {
    return `You must specify database name and id. Available names are ${CouchHelper.DBS.join(', ')}.`;
  }

  return Database.link.use(database).get(id);
});

DATABASE_CLI.addSubListener('info', async () => {
  // Show: Database info (document count in each)
  async function infoAbout(database: AbstractDatabase<any>) {
    return {
      count: await database.count().catch(e => 0),
      created: await database.isCreated(),
    };
  }
  
  function formatInfo(name: string, infos: { count: number, created: boolean }) {
    return `${name}\n\tcreated:\t${infos.created}\n\tdoc_count:\t${infos.count}`;
  }

  return '\n' + (
    await Promise.all(
      [
        [CouchHelper.USER_COLLECTION, Database.user],
        [CouchHelper.TOKEN_COLLECTION, Database.token],
        [CouchHelper.MOLECULE_COLLECTION, Database.molecule],
        [CouchHelper.STASHED_MOLECULE_COLLECTION, Database.stashed],
        [CouchHelper.RADIUS_COLLECTION, Database.radius],
        [CouchHelper.LIPID_COLLECTION, Database.lipid],
      ]
      .map(async e => formatInfo(e[0] as string, await infoAbout(e[1] as AbstractDatabase<any>)))
    )
  ).join('\n');
});

export default DATABASE_CLI;
