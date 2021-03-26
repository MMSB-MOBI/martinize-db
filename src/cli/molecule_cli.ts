import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import MoleculeOrganizer from "../MoleculeOrganizer";
import { Molecule, StashedMolecule } from "../Entities/entities";
import { parser_files } from "../routes/molecule/parser/parser_files";
import { CreateMoleculeFromJson, InfosJson } from "../routes/molecule/CreateMoleculeJson";

const MOLECULE_CLI = new CliListener(
  CliHelper.formatHelp("molecule", {
    commands: {
      list: 'List registred molecules',
      'get <id>': 'Get details about molecule <id>',
      'wipe <id>/all': 'Delete registred molecule <id> / all molecules',
    },
    onNoMatch: "Command is incorrect. Type \"molecule\" for help.",
  })
);

MOLECULE_CLI.command('list', async () => {
  const mols = await Database.molecule.all();
  const stashed = await Database.stashed.all();

  if (!mols.length && !stashed.length) {
    return `This server does not contain any molecule.`;
  }

  const normal = `Available molecules are \n- ${mols.map(m => m._id).join('\n- ')}`;
  const stash = `Available stashed molecules are \n- ${stashed.map(m => m._id).join('\n- ')}`;

  return (mols.length ? normal : "") + (mols.length && stashed.length ? "\n" : "") + (stashed.length ? stash : "");
});

MOLECULE_CLI.command('get', async rest => {
  rest = rest.trim();

  if (!rest) {
    return `Please specify a molecule id.`;
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  try {
    return Database.molecule.get(rest);
  } catch (e) {
    return Database.stashed.get(rest);
  }
});

MOLECULE_CLI.command('wipe', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a molecule id or "all".`;
  }

  if (rest === "all") {
    await Database.delete(CouchHelper.MOLECULE_COLLECTION);
    await Database.delete(CouchHelper.STASHED_MOLECULE_COLLECTION);
    await MoleculeOrganizer.removeAll();
    await Database.create(CouchHelper.MOLECULE_COLLECTION);
    await Database.create(CouchHelper.STASHED_MOLECULE_COLLECTION);
    return `Molecule database is wiped`;
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  let mol: Molecule | undefined = undefined;
  let stash: StashedMolecule | undefined = undefined;
  try {
    mol = await Database.molecule.get(rest);
  } catch (e) {
    try {
      stash = await Database.stashed.get(rest);
    } catch {
      return `Unable to get molecule (${rest})`;
    }
  }

  if (mol) {
    await MoleculeOrganizer.remove(mol.files);
    return Database.molecule.delete(mol);
  }
  else if (stash) {
    await MoleculeOrganizer.remove(stash.files);
    return Database.stashed.delete(stash);
  }
  return `Unable to find molecule.`
});

export let BATCH_MOLECULES: InfosJson[];

MOLECULE_CLI.command('batch', async rest=> {
  rest = rest.trim();

  if (!rest){
    return 'please specify a molecule files path';
  } else {
    BATCH_MOLECULES = parser_files(rest);
  }
});

MOLECULE_CLI.command('push', async => {
  CreateMoleculeFromJson(BATCH_MOLECULES);
  console.log('done');
})

export default MOLECULE_CLI;
