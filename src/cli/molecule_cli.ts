import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import MoleculeOrganizer from "../MoleculeOrganizer";
import { Molecule, StashedMolecule } from "../Entities/entities";
import { parser_files } from "../routes/molecule/parser/parser_files";
import { CreateMoleculeFromJson, InfosJson } from "../routes/molecule/CreateMoleculeJson";
import logger from "../logger";
import { create_top_in_dir } from "../routes/molecule/parser/create_topFile";

import cliFileSuggestor from '@interactive-cli-helper/file-suggestor';
import { ErrorType } from "../Errors";
import { GoTerms } from "../types";
import { CONNECTED_USER_CLI } from "./user_cli";
const fs = require('fs');

const MOLECULE_CLI = new CliListener(
  CliHelper.formatHelp("molecule", {
    commands: {
      list: 'List registred molecules',
      'get <id>': 'Get details about molecule <id>',
      'wipe <id>/all': 'Delete registred molecule <id> / all molecules',
      'load <path>': 'Load in memory all the molecules in the directory to insert them in the database',
      'push': 'Insert the molecules in memory in the database',
      'top <path>': 'Create top files for all the molecules in the directory',
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

MOLECULE_CLI.command('load', rest => {
  if (!CONNECTED_USER_CLI) {
    return 'Please connect before using this command by using user connect';
  }
  else {
    if (CONNECTED_USER_CLI.role != 'admin') {
      return 'You must be admin to add molecules to the database'
    }
    else {
      // TODO ajouter autres types de goterms
      

      if (!rest){
        return 'please specify a molecule files path';
      
      } else {

        let params = rest.split(' ');

        if (params.length < 2) {
          return 'please specify a type of molecules inserted (lipids, sugars)';
        
        } else {
          let path = params[0].trim();
          let type = params[1].trim();
          let go : keyof typeof GoTerms;

          if (!fs.existsSync(path)) {
            return 'Path does not exist';
          }
          else {
            if (type == 'lipids') {
              go = 'lipids';
            }
            else {
              return 'Invalid type of moecule';
            }
  
            try {
              BATCH_MOLECULES = parser_files(path, go);
              logger.info('load done');
            } catch (e) {
              logger.warn(e.data);   
            }
          }
        }
      }
    }
  }
}, {
  onSuggest: cliFileSuggestor,
});



MOLECULE_CLI.command('push', async() => {
  if (!CONNECTED_USER_CLI) {
    return 'Please connect before using this command by using user connect';
  }
  else {
    if (BATCH_MOLECULES) {
      try {
        await CreateMoleculeFromJson(BATCH_MOLECULES);
        logger.info('push done');
      } catch (e) {
        logger.warn(e.data);
      }
    } else {
      logger.warn('Missing molecules in memory. Please insert them by using molecule load.')
    }
  }
});



MOLECULE_CLI.command('top', rest => {
  rest = rest.trim();

  if (!fs.existsSync(rest)) {
    return 'Path does not exist';
  }
  else {
    try {
      create_top_in_dir(rest);
      logger.info('top files created');
    } catch (e) {
      logger.warn(e.data);
    }
  }
}, {
  onSuggest: cliFileSuggestor,
})

export default MOLECULE_CLI;
