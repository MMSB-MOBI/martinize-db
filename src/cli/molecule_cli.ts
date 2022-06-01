import CliHelper, { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import MoleculeOrganizer from "../MoleculeOrganizer";
import { Molecule, StashedMolecule } from "../Entities/entities";
import { parser_files, completeItpFiles, decodeCategory } from "../routes/molecule/parser/parser_files";
import { CreateMoleculeFromJson, InfosJson } from "../routes/molecule/CreateMoleculeJson";
import logger from "../logger";
import { create_top_in_dir } from "../routes/molecule/parser/create_topFile";

import cliFileSuggestor from '@interactive-cli-helper/file-suggestor';
import { ErrorType } from "../Errors";
import { GoTerms } from "../types";
import { CONNECTED_USER_CLI } from "./user_cli";
import { correctVersions} from '../routes/molecule/tmp_version'
const fs = require('fs');

const MOLECULE_CLI = new CliListener(
  CliHelper.formatHelp("molecule", {
    commands: {
      list: 'List registred molecules',
      'longlist': 'list registered molecules with informations',
      'stats': 'basic stats on molecules that are in the database',
      'get <id>': 'Get details about molecule <id>',
      'wipe <id>/all': 'Delete registred molecule <id> / all molecules',
      'load <path>': 'Load in memory all the molecules in the directory to insert them in the database',
      'push <log_path>': 'Insert the molecules in memory in the database, write recap in log_path',
      'top <path>': 'Create top files for all the molecules in the directory',
      'itp <path>' : 'Modify itp to include more informations. Only work with a specific directory organization.'
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


MOLECULE_CLI.command('longlist', async () => {
  const mols = await Database.molecule.all();
  const stashed = await Database.stashed.all();

  if (!mols.length && !stashed.length) {
    return `This server does not contain any molecule.`;
  }

  const normal = `Available molecules are \n- ${mols.map(m => 'id: ' + m._id + ', name: ' + m.name + ', version: '+ m.version +', type: ' + m.category).join('\n- ')}`;
  const stash = `Available stashed molecules are \n- ${stashed.map(m => 'id: ' + m._id + ', name: ' + m.name+ ', version: '+ m.version +', type: ' + m.category).join('\n- ')}`;

  return (mols.length ? normal : "") + (mols.length && stashed.length ? "\n" : "") + (stashed.length ? stash : "");
})

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

MOLECULE_CLI.command('stats', async () => {
  const {byCategories, byForceField, all} = await Database.molecule.stats()
  const molCount = `# Number of molecules : ${all.length}`
  const categories = `# ${Object.keys(byCategories).length} categories :`
  let categoriesPrint = ''
  for (const cat in byCategories){
    categoriesPrint += `${decodeCategory(cat)} : ${byCategories[cat].length}\n`
  }
  const ffHeader = `# ${Object.keys(byForceField).length} force fields :`
  let ffPrint = ''
  for (const ff in byForceField){
    ffPrint += `${ff} : ${byForceField[ff].length}\n`
  }
  return(`${molCount}\n\n${categories}\n${categoriesPrint}\n${ffHeader}\n${ffPrint}`)
})

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



export let BATCH_MOLECULES: InfosJson[] | undefined;

export class Excel {
  public static text: string = 'git repository = https://forge.ibcp.fr/glaunay/martini-molecule-repository,,,\nMolecule,No GRO,ITPs incompatible with GRO,Specific force_field not in server,Inserted successfuly\n'
}

MOLECULE_CLI.command('load', rest => {
  if (!CONNECTED_USER_CLI) {
    return 'Please connect before using this command by using user connect';
  }
  else {
    if (CONNECTED_USER_CLI.role != 'admin') {
      return 'You must be admin to add molecules to the database'
    }
    else {
      

      if (!rest){
        return 'please specify a molecule type and a files path';
      
      } else {

      let params = rest.split(' ');
      let path = params[0].trim();

      if (!fs.existsSync(path)) {
        return 'Path does not exist';
      }

        try {
          BATCH_MOLECULES = parser_files(path);
          logger.info(`batched ${BATCH_MOLECULES.length} molecules`);
          

        } catch (e) {
          logger.error("Error during load");
          console.log(e)
        }
      }
    }
  }
}, {
  onSuggest: cliFileSuggestor,
});

MOLECULE_CLI.command('tmpversion', async () => {
  await correctVersions()

})


MOLECULE_CLI.command('push', async rest => {
  rest = rest.trim();
  let logged = ''
  if(! rest) {
    logger.warn("No log file to write insertion recap")
  }
  else {
    logged = "# MAD molecules batch insertion"
  }

  if (!CONNECTED_USER_CLI) {
    return 'Please connect before using this command by using user connect';
  }
  else {
    if (BATCH_MOLECULES) {
      try {
        const recapInsertion = await CreateMoleculeFromJson(BATCH_MOLECULES);
        logger.info(`${recapInsertion.inserted.length} molecules inserted`);
        if (logged !== '') logged += "\n## Inserted \n" + recapInsertion.inserted.join("\n")
        logger.warn(`Molecules not inserted :`)
        if (logged !== '') logged += "\n## Not inserted\n"
        for(const reason in recapInsertion.not_inserted){
          console.log('##', reason)
          if (logged !== '') logged += `\n### ${reason}\n`
          if(logged === '') console.log(recapInsertion.not_inserted[reason].join("\n"))
          else {
            console.log(recapInsertion.not_inserted[reason].length)
            logged += recapInsertion.not_inserted[reason].join("\n")
          }
        }

        if(logged !== '') fs.writeFileSync(rest, logged)
        
        //logger.debug(Excel.text);
        //fs.writeFileSync('/home/achopin/Documents/molecules.csv', Excel.text);
      } catch (e) {
        logger.warn(e.data !== undefined ? e.data.message : e);
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

MOLECULE_CLI.command('itp', path =>  {
  path = path.trim()
  if (!fs.existsSync(path)){
    return 'Path does not exist'
  }
  try {
    completeItpFiles(path, true)
  } catch(e) {
    logger.error('Error while complete itp')
  }
}, {
  onSuggest: cliFileSuggestor,
})

export default MOLECULE_CLI;
