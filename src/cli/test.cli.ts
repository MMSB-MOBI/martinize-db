import CliHelper, { CliListener } from "interactive-cli-helper";
import { Martinizer } from "../Builders/Martinizer";
import path from "path";
import { promises as FsPromise, fstat } from 'fs';

const TEST_CLI = new CliListener(CliHelper.formatHelp('test', {
  map: 'Create a map with the web server with the specified file',
  ccmap: 'Create a map with the Python CCMap package with the specified file',
  'elastic-bonds': 'Generate the elastic bonds definitions with the given folder. It must contain a TOP file, PDB file and ITP(s)',
  'go-sites': 'Generate the elastic bonds definitions with the given folder. It must contain a ITPs files describing a Go model',
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

export default TEST_CLI;
