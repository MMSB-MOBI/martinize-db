import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";

const DATABASE_CLI = new CliListener(
  CliHelper.formatHelp("database", {
    'create <name>/all': 'Create a single or all databases. Available names are: ' + CouchHelper.DBS.join(', ') + '.',
    'wipe <name>/all': 'Delete a single or all databases.',
  })
);

DATABASE_CLI.addSubListener('create', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a database name or "all".`;
  }

  if (rest === "all") {
    return Database.createAll();
  }
  if (!CouchHelper.DBS.includes(rest)) {
    return `This database name is not authorized. Available names are ${CouchHelper.DBS.join(', ')}.`;
  }
  return Database.create(rest);
});

DATABASE_CLI.addSubListener('wipe', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a database name or "all".`;
  }

  if (rest === "all") {
    return Database.deleteAll();
  }
  if (!CouchHelper.DBS.includes(rest)) {
    return `This database name is not authorized. Available names are ${CouchHelper.DBS.join(', ')}.`;
  }
  return Database.delete(rest);
});

export default DATABASE_CLI;
