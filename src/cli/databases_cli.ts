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

export default DATABASE_CLI;
