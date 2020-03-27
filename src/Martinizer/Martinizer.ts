import { exec } from 'child_process';
import fs, { promises as FsPromise } from 'fs';
import os from 'os';
import Errors, { ErrorType } from '../Errors';
import axios, { AxiosResponse } from 'axios';
import FormData from 'form-data';
import path from 'path';
import TarStream from 'tar-stream';
import zlib from 'zlib';
import { ArrayValues, fileExists } from '../helpers';
import logger from '../logger';
import RadiusDatabase from '../Entities/RadiusDatabase';
import readline from 'readline';
import { FORCE_FIELD_DIR } from '../constants';

const DSSP_PATH = "/Users/alki/opt/anaconda3/bin/mkdssp";
const CREATE_GO_PATH = "/Users/alki/IBCP/create_goVirt.py";
const CREATE_MAP_PATH = path.resolve(__dirname, "../../utils/get_map.py");
const CONECT_PDB_PATH = path.resolve(__dirname, "../../utils/create_conect_pdb.sh");
const CONECT_MDP_PATH = path.resolve(__dirname, "../../utils/run.mdp");

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

const MARTINIZE_POSITIONS = ['none', 'all', 'backbone'] as const;
export type MartinizePosition = ArrayValues<typeof MARTINIZE_POSITIONS>;

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
  position: MartinizePosition;

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
  isMartinizePosition(value: any) : value is MartinizePosition {
    return MARTINIZE_POSITIONS.includes(value);
  }
  
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
  
      const exists = await fileExists(dir);
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
  
      const exists_pdb = await fileExists(pdb_file);
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

      // Modify system.top to include the right ITPs !
      const { top: full_top, itps: full_itps } = await this.createTopFile(dir, dir + '/system.top', itp_files, full.ff);

      // Generate the right PDB file (with conect entries)
      const pdb_with_conect = await this.createPdbWithConect(pdb_file, full_top, dir, false);

      return {
        pdb: pdb_with_conect,
        itps: full_itps,
        top: full_top,
      };
    } finally {
      await FsPromise.rename(with_ext, basename);
    }
  }

  /**
   * Modify the original top file in order to include the right ITPs, 
   * and create a link of the force field martini file into current directory.
   * 
   * ITPs in {itps_path} MUST be in {current_directory} !
   * 
   * Returns new TOP filename and all the used ITPs to generate top.
   */
  async createTopFile(current_directory: string, original_top_path: string, itps_path: string[], force_field: string) {
    let itps_ff = RadiusDatabase.FORCE_FIELD_TO_FILE_NAME[force_field];

    if (!itps_ff) {
      throw new ReferenceError("Your force field is invalid: can't find related ITPs.");
    }

    itps_ff = typeof itps_ff === 'string' ? [itps_ff] : itps_ff;

    const itps = [...itps_path, ...itps_ff.map(e => FORCE_FIELD_DIR + e)];
    const base_ff_itps = [] as string[];

    // Create everysym link
    for (const itp of itps_ff) {
      const itp_path = FORCE_FIELD_DIR + itp;
      const name = itp;
      const dest = path.resolve(current_directory + "/" + name);
      base_ff_itps.push(name);

      await FsPromise.symlink(itp_path, dest);
    }

    const real_itps = [...base_ff_itps, ...itps_path.map(e => path.basename(e))];
    const top = current_directory + "/full.top";

    const top_write_stream = fs.createWriteStream(top);
    
    // Add every #include first
    for (const itp of real_itps) {
      top_write_stream.write(`#include "${itp}"\n`);
    }

    const top_read_stream = readline.createInterface({
      input: fs.createReadStream(original_top_path),
      crlfDelay: Infinity,
    });

    // Remove every #include line
    for await (const line of top_read_stream) {
      if (line.startsWith('#include')) {
        continue;
      }
      top_write_stream.write(line + '\n');
    }

    top_write_stream.close();

    return {
      top,
      itps,
    };
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

  /**
   * Create the conect entries of the desired PDB/GRO.
   * 
   * This function should be done ONCE: After a Martinize Run / A INSANE Run / A molecule insert in database
   * 
   * Don't do it at each call!
   * 
   * Need the TOP topology file. 
   * ITP includes should be able to be resolved, use the {base_directory} parameter
   * in order to set the used current directory path.
   */
  async createPdbWithConect(pdb_or_gro_filename: string, top_filename: string, base_directory: string, remove_water = false) {
    const command = `${CONECT_PDB_PATH} "${pdb_or_gro_filename}" "${top_filename}" "${CONECT_MDP_PATH}" ${remove_water ? "--remove-water" : ""}`;
    const pdb_out = base_directory + "/output-conect.pdb";

    const [, err] = await new Promise((resolve, reject) => {
      exec(command, { cwd: base_directory }, (err, stdout, stderr) => {
        if (err) {
          reject(err);
          return;
        }
        resolve([stdout, stderr]);
      });
    }) as [string, string];

    if (err) {
      logger.warn("An error occured during run of pdb conect: " + err);
    }

    const exists = await FsPromise.access(pdb_out, fs.constants.F_OK).then(() => true).catch(() => false);

    if (!exists) {
      throw new Error("PDB could not be created for an unknown reason. Check the logs.");
    }

    return pdb_out;
  }
}();
