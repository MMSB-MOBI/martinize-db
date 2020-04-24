import CliHelper, { CliListener } from "interactive-cli-helper";
import { Martinizer } from "../Builders/Martinizer";
import path from "path";
import { promises as FsPromise, fstat } from 'fs';
import { LIPIDS_ROOT_DIR } from "../constants";
import { Lipid } from "../Entities/entities";
import { Database } from "../Entities/CouchHelper";
import RadiusDatabase from "../Entities/RadiusDatabase";

const TEST_CLI = new CliListener(CliHelper.formatHelp('test', {
  map: 'Create a map with the web server with the specified file',
  ccmap: 'Create a map with the Python CCMap package with the specified file',
  'elastic-bonds': 'Generate the elastic bonds definitions with the given folder. It must contain a TOP file, PDB file and ITP(s)',
  'go-sites': 'Generate the elastic bonds definitions with the given folder. It must contain a ITPs files describing a Go model',
  'add-lipid': 'Add a lipid, followed by a lipid filename (search in {workdir}/lipids/2_2).',
  'lipid-prefixes': 'Get lipid directory prefix name by force field.',
  'auto-import-lipid <forceField>': 'Import all the itp files inside a lipid directory found by its related force field.'
}));

TEST_CLI.addSubListener('map', rest => Martinizer.getMap(rest));

TEST_CLI.addSubListener('ccmap', rest => Martinizer.getCcMap(rest));

TEST_CLI.addSubListener('elastic-bonds', async rest => {
  const itps: string[] = [];
  let top_file: string = "", pdb_file: string = "";

  rest = path.resolve(rest) + "/";
  console.log("Given path:", rest);

  for (const element of await FsPromise.readdir(rest)) {
    const basename = path.basename(element);
    if (basename.endsWith('.pdb')) {
      pdb_file = rest + basename;
    }
    else if (basename.endsWith('.top')) {
      top_file = rest + basename;
    }
    else if (basename.endsWith('.itp')) {
      itps.push(rest + basename);
    }
  }

  if (!top_file || !pdb_file || !itps.length) {
    return "Given folder must have a ITP file, TOP file and PDB file.";
  }

  const relations = await Martinizer.computeElasticNetworkBounds(top_file, itps);

  await FsPromise.writeFile(rest + 'relations.json', JSON.stringify(relations, null, 2));

  return "Relations has been written to '" + (rest + 'relations.json') + "'.";
});

TEST_CLI.addSubListener('go-sites', async rest => {
  const itps: string[] = [];
  let top_file: string = "";

  rest = path.resolve(rest) + "/";
  console.log("Given path:", rest);

  for (const element of await FsPromise.readdir(rest)) {
    const basename = path.basename(element);
    if (basename.endsWith('.top')) {
      top_file = rest + basename;
    }
    else if (basename.endsWith('.itp')) {
      itps.push(rest + basename);
    }
  }

  if (!itps.length) {
    return "Given folder must have ITP files.";
  }

  const relations = await Martinizer.computeGoModelBounds(top_file, itps);

  await FsPromise.writeFile(rest + 'relations.json', JSON.stringify(relations, null, 2));

  return "Relations has been written to '" + (rest + 'relations.json') + "'.";
});

TEST_CLI.addSubListener('add-lipid', async rest => {
  const name = LIPIDS_ROOT_DIR + "2_2/" + rest + '.itp';

  const lipid: Lipid = {
    id: '',
    name: rest,
    itp: await FsPromise.readFile(name, 'utf-8'),
  };

  await Database.lipid.add(lipid, 'martini22');

  return lipid;
});

TEST_CLI.addSubListener('lipid-prefixes', () => {
  return CliHelper.formatHelp("    force field  prefix", RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION); 
});

TEST_CLI.addSubListener('auto-import-lipid', async force_field => {
  const prefix = RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field];

  if (!prefix) {
    return "Unknown force field.";
  }

  const dir = LIPIDS_ROOT_DIR + prefix + "/";

  for (const itp of await FsPromise.readdir(dir)) {
    if (!itp.endsWith('.itp')) {
      continue;
    }

    const lipid: Lipid = {
      id: '',
      name: itp.split('.itp')[0].toLocaleUpperCase(),
      itp: await FsPromise.readFile(dir + itp, 'utf-8'),
    };

    await Database.lipid.add(lipid, force_field);
  }
});

export default TEST_CLI;
