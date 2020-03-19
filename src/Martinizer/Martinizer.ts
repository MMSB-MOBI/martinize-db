import { exec } from 'child_process';
import fs, { promises as FsPromise } from 'fs';
import os from 'os';
import Errors, { ErrorType } from '../Errors';
import axios, { AxiosResponse } from 'axios';
import FormData from 'form-data';
import path from 'path';
import TarStream from 'tar-stream';
import zlib from 'zlib';

const DSSP_PATH = "/Users/alki/opt/anaconda3/bin/mkdssp";
const CREATE_GO_PATH = "/Users/alki/IBCP/create_goVirt.py";
const CREATE_MAP_PATH = path.resolve(__dirname, "../../utils/get_map.py");

interface CCMapResChain {
  resID: string;
  chainID: string;
}

interface CCMapWithDistance extends CCMapResChain {
  distance: number;
}

interface ContactMapCCMap {
  type: 'contactList';
  data: {
    root: CCMapResChain;
    partners: CCMapWithDistance[];
  }[];
}

export interface MartinizeSettings {
  /** PDB file path */
  input: string;
  
  /** Ignore residues */
  ignore?: string[];

  /** Ignore hydrogens */
  ignh?: boolean;

  /** Force field */
  ff: string;

  /** Position restrains */
  position: 'none' | 'all' | 'backbone';

  /** Position restrain force const */
  posref_fc?: number;

  /** Use collagen */
  collagen?: boolean;

  /** Use dihedral */
  dihedral?: boolean;

  /** Elastic bounds */
  elastic?: boolean;
  /** Elastic force const */
  ef?: number;
  /** Elastic lower bound */
  el?: number;
  /** Elastic upper bound */
  eu?: number;
  /** Elastic decay alpha */
  ea?: number;
  /** Elastic decay power */
  ep?: number;
  /** Elastic remover minimum force */
  em?: number;
  /** List of bead names for elastic bound (comma separated in martinize) */
  eb?: string[];

  /** Use govs */
  use_go_virtual_sites?: boolean;

  /** Set neutral termini */
  neutral_termini?: boolean;
  /** Apply side chains corrections */
  side_chain_fix?: boolean;
  /** Cystein bounds */
  cystein_bridge?: string;
}

export const Martinizer = new class Martinizer {
  async run(settings: Partial<MartinizeSettings>) {
    const full: MartinizeSettings = Object.assign({}, {
      input: '',
      ff: 'martini22',
      position: 'none'
    }, settings);

    if (!full.input.trim()) {
      throw new Error("Invalid input file");
    }

    const basename = full.input;
    const with_ext = basename + '.pdb';

    // Check dssp ps
    let command_line = "martinize2 -f " + with_ext + " -x output.pdb -o system.top -dssp " + DSSP_PATH + " -ff " + full.ff + " -p " + full.position + " ";

    if (full.ignore) {
      command_line += " " + full.ignore.join(',');
    }
    if (full.ignh) {
      command_line += " -ignh";
    }

    if (full.posref_fc) {
      command_line += " -pf " + full.posref_fc.toString();
    }
    if (full.collagen) {
      command_line += " -collagen ";
    }
    if (full.dihedral) {
      command_line += " -ed ";
    }
    if (full.elastic) {
      command_line += " -elastic ";
    }
    if (full.ef) {
      command_line += " -ef " + full.ef.toString();
    }
    if (full.el) {
      command_line += " -el " + full.el.toString();
    }
    if (full.eu) {
      command_line += " -eu " + full.eu.toString();
    }
    if (full.ea) {
      command_line += " -ea " + full.ea.toString();
    }
    if (full.ep) {
      command_line += " -ep " + full.ep.toString();
    }
    if (full.em) {
      command_line += " -em " + full.em.toString();
    }
    if (full.eb) {
      command_line += " -eb " + full.eb.toString();
    }
    if (full.use_go_virtual_sites) {
      command_line += " -govs-include ";
    }
    if (full.neutral_termini) {
      command_line += " -nt ";
    }
    if (full.side_chain_fix) {
      command_line += " -scfix ";
    }
    if (full.cystein_bridge) {
      command_line += " -cys " + full.cystein_bridge;
    }

    await FsPromise.rename(basename, with_ext);

    try {
      // Run the command line
      const tmp_dir = os.tmpdir();
      const dir = await FsPromise.mkdtemp(tmp_dir + "/");
  
      const exists = await FsPromise.access(dir, fs.constants.F_OK).then(() => true).catch(() => false);
      if (!exists) {
        await FsPromise.mkdir(dir, { recursive: true });
      }
  
      await new Promise((resolve, reject) => {
        exec(command_line, { cwd: dir }, (err, stdout, stderr) => {
          if (err) {
            reject(err);
            return;
          }
          resolve([stdout, stderr]);
        });
      }) as [string, string];
  
      // Scan for files in result dir
      let pdb_file: string = dir + "/output.pdb";
      const itp_files: string[] = [];
      console.log("Tmp directory for Martinize job:", dir);
  
      const exists_pdb = await FsPromise.access(pdb_file, fs.constants.F_OK).then(() => true).catch(() => false);
      if (!exists_pdb) {
        return Errors.throw(ErrorType.Server, { error: "Unable to find created pdb." });
      }

      // If go mode, we should compute map + run a python script to refresh ITPs files.
      if (settings.use_go_virtual_sites) {
        // Must create the go sites
        const [out, ] = await new Promise((resolve, reject) => {
          exec('cut -f1 -d \' \' output.pdb | grep -c ATOM', { cwd: dir }, (err, stdout, stderr) => {
            if (err) {
              reject(err);
              return;
            }
            resolve([stdout, stderr]);
          });
        }) as [string, string];
        
        // GET THE MAP FILE FROM A CUSTOM WAY.
        // Use the original pdb file !!!
        const map_filename = await this.getMap(with_ext, dir);

        // Ensure to have the script at the correct place...
        await new Promise((resolve, reject) => {
          exec(`python ${CREATE_GO_PATH} -s output.pdb -f ${map_filename} --Natoms ${out.trim()} --moltype molecule_0`, { cwd: dir }, (err, stdout, stderr) => {
            if (err) {
              // Sometimes, the program crashes, I don't know why. To investigate
              reject(err);
              return;
            }
            resolve([stdout, stderr]);
          });
        }) as [string, string];

        // Ok, ITP file refreshed.
      }

      for (const file of await FsPromise.readdir(dir)) {
        if (file.endsWith('.itp')) {
          itp_files.push(dir + "/" + file);
        }
      }

      // TODO: modify system.top to include the right ITPs !

      return {
        pdb: pdb_file,
        itps: itp_files,
        top: dir + '/system.top',
      };
    } finally {
      await FsPromise.rename(with_ext, basename);
    }
  }

  async getCcMap(pdb_filename: string, use_tmp_dir?: string) {
    const [strmap, ] = await new Promise((resolve, reject) => {
      exec(`python ${CREATE_MAP_PATH} -f ${pdb_filename}`, (err, stdout, stderr) => {
        if (err) {
          reject(err);
          return;
        }
        resolve([stdout, stderr]);
      });
    }) as [string, string];

    const map: ContactMapCCMap = JSON.parse(strmap);

    // Save the map file inside a temporary directory
    if (!use_tmp_dir) {
      const tmp_dir = os.tmpdir();
      use_tmp_dir = await FsPromise.mkdtemp(tmp_dir + "/");
    }

    // Prepare the write stream for filesave
    const map_filename = path.resolve(use_tmp_dir + '/output.map');
    const map_stream = fs.createWriteStream(map_filename);

    map_stream.write(`            I1  AA  C I(PDB)    I2  AA  C I(PDB)    DISTANCE       CMs    rCSU    aSurf    rSurf    nSurf\n`);
    map_stream.write(`==========================================================================================================\n`);

    for (const atom of map.data) {
      for (const partner of atom.partners) {
        map_stream.write(`R      1     1  XXX ${atom.root.chainID}    ${atom.root.resID}        2  XXX ${partner.chainID}    ${partner.resID}       ${partner.distance}     1 1 1 1    16   2.6585   0.0000  60.5690\n`);
      }
    }

    map_stream.close();

    return map_filename;
  }

  protected findUrlInRedirect(data: string) {
    // Find the redirect in page
    const rest = data.split('Content="0; URL=')[1];
    if (!rest) {
      throw new Error('Unable to find map URL.');
    }
    
    return 'http://info.ifpan.edu.pl/~rcsu/rcsu/' + rest.split('">')[0];
  }

  async getMap(pdb_filename: string, use_tmp_dir?: string) {
    // Create the initial form data job
    const form = new FormData;

    form.append('filename', fs.createReadStream(pdb_filename));
    form.append('radii', 'tsai');
    form.append('fib', '14');
    form.append('allchains', '1');
    form.append('PDB_ID', '');

    // Start the job
    let res: AxiosResponse<string> = await axios.post('http://info.ifpan.edu.pl/~rcsu/rcsu/prepare.php', form, {
      headers: {
        ...form.getHeaders()
      },
      responseType: 'text'
    });

    // Follow the first redirect
    res = await axios.get(this.findUrlInRedirect(res.data), { responseType: 'text' });

    // Follow the second redirect, to the TAR file.
    // To get the tar, you should only replace .html to .tgz
    const url = this.findUrlInRedirect(res.data).replace('.html', '.tgz');

    // Save the map file inside a temporary directory
    if (!use_tmp_dir) {
      const tmp_dir = os.tmpdir();
      use_tmp_dir = await FsPromise.mkdtemp(tmp_dir + "/");
    }

    // Prepare the write stream for filesave
    const map_filename = path.resolve(use_tmp_dir + '/output.map');
    const map_stream = fs.createWriteStream(map_filename);

    // Download the tgz via a stream
    const map_response = await axios.get(url, { responseType: 'stream' });
    
    // Prepare the targz extractor, then pipe it to response stream
    const extractor = TarStream.extract();

    map_response.data
      .pipe(zlib.createGunzip())
      .pipe(extractor);

    // Assign stream download to extract only the .map file,
    // pipe the file content to the {map_stream}
    await new Promise((resolve, reject) => {
      extractor.on('entry', (header, stream, next) => {
        // {stream} is the file body, 
        // call {next} entry is read

        if (header.name.endsWith('.map')) {
          stream.on('data', chunk => {
            map_stream.write(chunk);
          });

          stream.on('end', () => {
            // ready for next entry
            next(); 
          });

          // Start the read stream
          stream.resume();
        }
        else {
          next();
        }
      });
      
      extractor.on('finish', resolve);
      extractor.on('error', reject);
    });

    map_stream.close();

    return map_filename;
  }
}();
