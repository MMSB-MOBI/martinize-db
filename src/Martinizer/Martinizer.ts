import axios, { AxiosResponse } from 'axios';
import { exec, ExecException } from 'child_process';
import FormData from 'form-data';
import fs, { promises as FsPromise } from 'fs';
import path from 'path';
import readline from 'readline';
import TarStream from 'tar-stream';
import zlib from 'zlib';
import { FORCE_FIELD_DIR } from '../constants';
import RadiusDatabase from '../Entities/RadiusDatabase';
import Errors, { ErrorType } from '../Errors';
import { ArrayValues, fileExists } from '../helpers';
import logger from '../logger';
import TmpDirHelper from '../TmpDirHelper/TmpDirHelper';
import { TopFile, ItpFile } from '../ItpParser';

const DSSP_PATH = "/Users/alki/opt/anaconda3/bin/mkdssp";
const CREATE_GO_PATH = "/Users/alki/IBCP/create_goVirt.py";
const CREATE_MAP_PATH = path.resolve(__dirname, "../../utils/get_map.py");
const CONECT_PDB_PATH = path.resolve(__dirname, "../../utils/create_conect_pdb.sh");
const CONECT_MDP_PATH = path.resolve(__dirname, "../../utils/run.mdp");

/**
 * Tuple of two integers: [{from} atom index, {to} atom index]
 */
export type ElasticOrGoBounds = [number, number];

interface ContactMapCCMap {
  type: 'atomic';
  data: [ 
    /* Atom name, residue name, residue number, chain */
    [string, string, string, string],
    /* Atom name, residue name, residue number, chain */
    [string, string, string, string],
    /* Distance */
    number,
  ][];
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
  protected MAX_JOB_EXECUTION_TIME = 5 * 60 * 1000;

  isMartinizePosition(value: any) : value is MartinizePosition {
    return MARTINIZE_POSITIONS.includes(value);
  }
  
  /**
   * Create a martinize run.
   * Returns created path to created PDB, TOP and ITP files.
   */
  async run(settings: Partial<MartinizeSettings>) {
    const full: MartinizeSettings = Object.assign({}, {
      input: '',
      ff: 'martini22',
      position: 'none'
    }, settings);
    
    if (!full.input.trim()) {
      throw new Error("Invalid input file");
    }

    logger.debug(`Starting a Martinize run for ${path.basename(full.input)}.`);

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
      const dir = await TmpDirHelper.get();

      logger.debug("[MARTINIZER-RUN] Tmp directory for Martinize job: " + dir);
  
      const exists = await fileExists(dir);
      if (!exists) {
        await FsPromise.mkdir(dir, { recursive: true });
      }
      
      try {
        await new Promise((resolve, reject) => {
          const stdout_martinize = fs.createWriteStream(dir + '/martinize.stdout');
          const stderr_martinize = fs.createWriteStream(dir + '/martinize.stderr');
          const MTZ_STEP_INDICATOR = 'INFO - step - ';

          const child = exec(command_line, { cwd: dir }, err => {
            stdout_martinize.close();
            stderr_martinize.close();
            child.stderr?.removeAllListeners();
  
            if (err) {
              reject(err);
              return;
            }
            resolve();
          });

          child.stdout?.pipe(stdout_martinize);
          child.stderr?.on('data', (chunk: Buffer) => {
            const line = chunk.toString().trim();

            if (line.startsWith(MTZ_STEP_INDICATOR)) {
              // Get the step name (every step is ended with a dot, so line.length-1 remove it)
              const step_info = line.slice(MTZ_STEP_INDICATOR.length, line.length - 1);

              // Todo send to socket.io ?
              logger.debug("[MARTINIZER-RUN] [Martinize Step] " + step_info);
            }

            stderr_martinize.write(chunk);
          });
        });
      } catch (e) {
        const err: ExecException = e;

        return Errors.throw(ErrorType.MartinizeRunFailed, { error: "Martinize as failed with an exit code.", misc: err, cwd: dir });
      }
  
      // Scan for files in result dir
      let pdb_file: string = dir + "/output.pdb";
      const itp_files: string[] = [];
  
      const exists_pdb = await fileExists(pdb_file);
      if (!exists_pdb) {
        return Errors.throw(ErrorType.MartinizeRunFailed, { error: "Unable to find created pdb." });
      }

      logger.debug(`[MARTINIZER-RUN] Generated PDB is found, run should be fine.`);

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
        // todo change (ccmap create way too much distances, so the shell script takes forever)
        const map_filename = await this.getCcMap(with_ext, dir);

        logger.debug("[MARTINIZER-RUN] Creating Go virtual bonds");

        // Ensure to have the script at the correct place...
        await new Promise((resolve, reject) => {
          const stdout = fs.createWriteStream(dir + '/contact-map.stdout');
          const stderr = fs.createWriteStream(dir + '/contact-map.stderr');

          const child = exec(`python ${CREATE_GO_PATH} -s output.pdb -f ${map_filename} --Natoms ${out.trim()} --moltype molecule_0`, { cwd: dir, maxBuffer: 1e9 }, err => {
            stdout.close();
            stderr.close();

            if (err) {
              // Sometimes, the program crashes, I don't know why. To investigate
              reject(err);
              return;
            }
            resolve();
          });

          child.stdout?.pipe(stdout);
          child.stderr?.pipe(stderr);
        });

        // Ok, ITP file refreshed.
      }

      for (const file of await FsPromise.readdir(dir)) {
        if (file.endsWith('.itp')) {
          itp_files.push(dir + "/" + file);
        }
      }

      // Modify system.top to include the right ITPs !
      logger.debug("[MARTINIZER-RUN] Creating full TOP file for Martinize built molecule.");
      // Full itps contain the force field itps.
      const { top: full_top, itps: full_itps } = await this.createTopFile(dir, dir + '/system.top', itp_files, full.ff);

      // Generate the right PDB file (with conect entries)
      logger.debug("[MARTINIZER-RUN] Creating PDB with CONECT entries for Martinize built molecule.");
      const pdb_with_conect = await this.createPdbWithConect(pdb_file, full_top, dir, false);

      logger.debug("[MARTINIZER-RUN] Run is complete, everything seems to be fine :)");

      return {
        pdb: pdb_with_conect,
        itps: itp_files,
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
    
    const includes: string[] = [];

    // Define the includes
    for (const itp of real_itps) {
      // Exclude the GO ITPs, they're already included in martini_304.itp
      if (itp.endsWith('VirtGoSites.itp') || itp.endsWith('go4view_harm.itp')) {
        continue;
      }

      includes.push(`#include "${itp}"`);
    }

    const top_write_stream = fs.createWriteStream(top);
    const top_read_stream = readline.createInterface({
      input: fs.createReadStream(original_top_path),
      crlfDelay: Infinity,
    });

    let includes_included = false;

    // Remove every #include line
    for await (const line of top_read_stream) {
      if (line.startsWith('#include')) {
        if (!includes_included) {
          // Include the hand-crafted includes
          top_write_stream.write(includes.join('\n') + '\n');
          includes_included = true;
        }

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

  /**
   * Get a contact map using ccmap python package.
   */
  async getCcMap(pdb_filename: string, use_tmp_dir?: string) {
    // Save the map file inside a temporary directory
    if (!use_tmp_dir) {
      use_tmp_dir = await TmpDirHelper.get();
      logger.debug(`Created tmp directory for ccmap: ${use_tmp_dir}.`);
    }

    const distances_file = use_tmp_dir + '/distances.json';

    logger.debug("[CCMAP] Calculating distances for backbone atoms.");

    // Extract only the CA atoms from original PDB, as they're the only to count for contact map.
    await new Promise((resolve, reject) => {
      exec(`grep 'CA' "${path.resolve(pdb_filename)}" > "${use_tmp_dir + "/_backbones.pdb"}"`, { maxBuffer: 1e9 }, err => {
        if (err) {
          reject(err);
          return;
        }

        resolve();
      });
    });

    // Compute contacts with the CA pdb
    await new Promise((resolve, reject) => {
      const stdout = fs.createWriteStream(use_tmp_dir + '/distances.stdout');
      const stderr = fs.createWriteStream(use_tmp_dir + '/distances.stderr');

      const child = exec(`python ${CREATE_MAP_PATH} -f "${use_tmp_dir}/_backbones.pdb" -o "${distances_file}"`, (err) => {
        stdout.close();
        stderr.close();

        if (err) {
          reject(err);
          return;
        }
        resolve();
      });

      child.stdout?.pipe(stdout);
      child.stderr?.pipe(stderr);
    }) as [string, string];

    const distances_exists = await fileExists(distances_file);

    if (!distances_exists) {
      throw new Error("Distance file not found. Did the ccmap run has been done fine? Check distances.stdout and distances.stderr in tmp dir.");
    }

    const map: ContactMapCCMap = JSON.parse(await FsPromise.readFile(distances_file, 'utf-8'));

    // Remove the distances.json

    logger.debug("[CCMAP] Creating fake rCSU contact map.");

    // Prepare the write stream for filesave
    const map_filename = path.resolve(use_tmp_dir + '/output.map');
    const map_stream = fs.createWriteStream(map_filename);

    map_stream.write(`            I1  AA  C I(PDB)    I2  AA  C I(PDB)    DISTANCE       CMs    rCSU    aSurf    rSurf    nSurf\n`);
    map_stream.write(`==========================================================================================================\n`);

    for (const [atom1, atom2, distance] of map.data) {
      const chain_1 = atom1[3].trim();
      const chain_2 = atom2[3].trim();

      const res_1 = atom1[2].trim();
      const res_2 = atom2[2].trim();

      map_stream.write(`R      1     1  XXX ${chain_1}    ${res_1}        2  XXX ${chain_2}    ${res_2}       ${distance}     1 1 1 1    16   2.6585   0.0000  60.5690\n`);
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

  /**
   * Get a contact map using the contact map web-server.
   */
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
      use_tmp_dir = await TmpDirHelper.get();
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

    await new Promise((resolve, reject) => {
      const stdout = fs.createWriteStream(base_directory + "/gromacs.stdout");
      const stderr = fs.createWriteStream(base_directory + "/gromacs.stderr");

      const child = exec(command, { cwd: base_directory, timeout: this.MAX_JOB_EXECUTION_TIME }, err => {
        stdout.close();
        stderr.close();

        if (err) {
          reject(err);
          return;
        }
        resolve();
      });

      child.stdout?.pipe(stdout);
      child.stderr?.pipe(stderr);
    });

    const exists = await FsPromise.access(pdb_out, fs.constants.F_OK).then(() => true).catch(() => false);

    if (!exists) {
      throw new Error("PDB could not be created for an unknown reason. Check the files gromacs.stdout and gromacs.stderr in directory " + base_directory + ".");
    }

    return pdb_out;
  }
  
  /**
   * Generate an object than contain "additionnal" bounds between atoms, that are created with the elastic network (option -elastic).
   * 
   * Specificy only ITP files that are generated by Martinize, not the force field itself!!
   * 
   * TODO: worker thread
   */
  async computeElasticNetworkBounds(top_file: string, itp_files: string[]) {
    logger.verbose("[ELASTIC-BUILD] Constructing elastic network bonds.");

    const bounds: ElasticOrGoBounds[] = [];

    logger.debug("[ELASTIC-BUILD] Reading TOP+ITP files.");
    const top = new TopFile(top_file, itp_files);
    await top.read();

    logger.debug(`[ELASTIC-BUILD] Available molecules: ${
      top.molecule_list
        .map(([name, itp]) => `${name} (${itp.length} time${itp.length > 1 ? 's' : ''})`)
        .join(', ')
    }`);

    // Incrementer for designating PDB line
    let i = 0;

    for (const [name, itps] of top.molecule_list) {
      logger.debug("[ELASTIC-BUILD] Reading molecule " + name + ".");

      // Get the number of atoms in a single chain of this molecule
      const atom_count = itps[0].atoms.filter(line => line && !line.startsWith(';')).length;
      let chain_n = 0;

      // There is one ITP per chain
      for (const itp of itps) {
        let should_read = false;
        chain_n++;
        
        logger.verbose(`[ELASTIC-BUILD] Chain ${chain_n} of molecule "${name}": ${itp.bonds.length} bond lines.`);

        for (const band of itp.bonds) {
          if (band.startsWith(';')) {
            // Find the elastic related comments
            if (
              band.startsWith('; Long elastic bonds for extended regions') ||
              band.startsWith('; Rubber band') || 
              band.startsWith('; Short elastic bonds for extended regions')
            ) {
              should_read = true;
            }
            // This is ALWAYS after rubber band/elastic bonds comments, we can stop here
            else if (
              band.startsWith('; Side chain bonds')
            ) {
              break;
            }
            else {
              should_read = false;
            }

            // Its a comment, skip it
            continue;
          }

          if (!should_read || !band) {
            continue;
          }
  
          // Here is a line with a band bound.
          // The two first numbers are the two concern atom by the bound.
          // Here: Atom 1 and 11 of the PDB file, part {itp_index}.
          // MARTINI22: 1 11 6 0.598 500.0
          // MARTINI30: 1 11 6 0.59901 500.0
  
          const [atom_from, atom_to, ] = band.split(/\s+/g);
  
          bounds.push([
            Number(atom_from) + i,
            Number(atom_to) + i,
          ]);
        }

        // Add atom count of this molecule to i
        i += atom_count;
      }
    }

    return bounds;
  }

  /**
   * Generate an object than contain "additionnal" bounds between atoms, that are created with the Go model (option -govs-include).
   * 
   * Specificy only ITP files that are generated by Martinize, not the force field itself!!
   * 
   * TODO: worker thread
   */
  async computeGoModelBounds(top_file: string, itp_files: string[], remove_duplicates = true) {
    /*
     *  Même si il y a deux chaînes, en mode go elles seront dans une seule molécule (normalement).
     *
     *  1: Lire les noms des atoms Go dans [atoms]: {molecule_type}_{i}
     *  2: Lire les associations go_atom => real atom index dans [virtual_sitesn] après comment "; Virtual go site"
     *  3: Lire le fichier {molecule_type}_go-table-VirtGoSites.itp qui définit les liaisons go_atom <=> go_atom
     *  4: Convertir ces liaisons go <=> go en "additionnal" bonds entre real atoms
     */

    logger.verbose("[GO-VIRT-SITES] Reading system topology.");
    const top = new TopFile(top_file, []);
    const molecule_types: string[] = [];

    await top.read();

    for (const molecule of top.getField('molecules')) {
      // molecule_0    1 (count is always 1 for Go model.)
      const [name, ] = molecule.split(ItpFile.BLANK_REGEX);
      molecule_types.push(name);
    }

    // Here, we store bounds
    const bounds: ElasticOrGoBounds[] = [];
    // Increment counter for bounds add.
    let i = 0;

    for (const molecule_type of molecule_types) {
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Finding files used to describe the Go model.`);

      // Find the [molecule_type].itp and [molecule_type]gotable ITP in ITP files
      const molecule_itp_index = itp_files.findIndex(e => path.basename(e) === molecule_type + '.itp');
      const go_table_index = itp_files.findIndex(e => path.basename(e) === molecule_type + '_go-table_VirtGoSites.itp');
  
      if (molecule_itp_index === -1 || go_table_index === -1) {
        logger.error(`[GO-VIRT-SITES] [${molecule_type}] Molecule ITP file or Go Virt Table not found.`);
        continue;
      }
  
      // Instanciate ITP
      const molecule_itp = new ItpFile(itp_files[molecule_itp_index]);
      const go_table = new ItpFile(itp_files[go_table_index]);

      // Read the molecule file
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Reading ITP files.`);
      await molecule_itp.read();
  
      const prefix = molecule_type + '_';
      /** Link go atom index to real atom index. */
      const index_to_real: { [index: number]: number } = {};
      /** Link go atom name to go atom index. */
      const name_to_index: { [name: string]: number } = {};

      // Step 1: Find atoms that name start by "{molecule_type}_" in category "atoms"
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Looking for virtual atoms.`);
      // Count all atoms, used to increment atom counter at the end of loop
      const all_atom_count = molecule_itp.atoms.filter(line => line && !line.startsWith(';')).length;

      for (const atom_line of molecule_itp.atoms) {
        // Typical line is : 
        // 2575 molecule_0_9       9 LYS CA  2575    0 
  
        const [index, name, ] = atom_line.split(ItpFile.BLANK_REGEX);
  
        if (name.startsWith(prefix)) {
          name_to_index[name] = Number(index);
        }
      }
  
      // Step 2: Associate go atom index => real atom index
      let seen_virt_comment = false;

      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Looking for virtual sites description.`);
      for (const virt_line of molecule_itp.virtual_sites) {
        if (virt_line.startsWith('; Virtual go site')) {
          seen_virt_comment = true;
          continue;
        }
  
        if (!seen_virt_comment) {
          continue;
        }
  
        // Typical line is: 
        // 2575 1    1
        const [go_index, , real_index] = virt_line.split(ItpFile.BLANK_REGEX);
        index_to_real[Number(go_index)] = Number(real_index);
      }
  
      const n_atoms = Object.keys(name_to_index).length;
      const n_sites = Object.keys(index_to_real).length;

      if (n_sites !== n_atoms) {
        // Print number of atoms only if useful.
        logger.verbose(`[GO-VIRT-SITES] [${molecule_type}] ${n_atoms} virtual atoms found.`);
        logger.verbose(`[GO-VIRT-SITES] [${molecule_type}] ${n_sites} virtual sites found.`);
        logger.warn(`[GO-VIRT-SITES] [${molecule_type}] Number of sites does not match number of atoms. Some bonds may not be linked correclty.`);
      }
  
      // Clean the ITP (we don't need it anymore)
      molecule_itp.dispose();

      // Read the go table
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Reading virtual Go sites table.`);
      await go_table.read();
  
      logger.verbose(`[GO-VIRT-SITES] [${molecule_type}] Atom bonds described: ${go_table.headlines.length - 2}.`);
  
      // Step 3+4: Read bonds between go atoms and associate them

      // To remove duplicates (that, unfortunately, exists...), we use a map of set
      const local_bonds: { [index: number]: Set<number> } = {}; 

      for (const line of go_table.headlines) {
        if (line.startsWith(';')) {
          continue;
        }
  
        // Typical line is (may begin by spaces.)
        // molecule_0_9  molecule_0_14    1  0.7369739126  9.4140000000  ;  24  36  0.827 
  
        // filter trim blank spaces created by regex
        const [name1, name2, ] = line.split(ItpFile.BLANK_REGEX).filter(e => e);
  
        const go_index_1 = name_to_index[name1], go_index_2 = name_to_index[name2];
  
        if (go_index_1 === undefined || go_index_2 === undefined) {
          logger.warn(`[GO-VIRT-SITES] [${molecule_type}] Undefined go indexes for names ${name1}-${name2}. This should not happen...`);
          continue;
        }
  
        const real_index_1 = index_to_real[go_index_1], real_index_2 = index_to_real[go_index_2];
  
        if (real_index_1 === undefined || real_index_2 === undefined) {
          logger.warn(`[GO-VIRT-SITES] [${molecule_type}] Undefined real indexes for names ${name1}(${go_index_1})-${name2}(${go_index_2}). This should not happen...`);
          continue;
        }
  
        // We add the bonds in the set
        if (remove_duplicates) {
          const [computed_1, computed_2] = real_index_1 < real_index_2 ? [
            real_index_1 + i, real_index_2 + i
          ] : [
            real_index_2 + i, real_index_1 + i
          ];
  
          if (computed_1 in local_bonds) {
            local_bonds[computed_1].add(computed_2);
          }
          else {
            local_bonds[computed_1] = new Set([computed_2]);
          }
        }
        else {
          bounds.push([
            real_index_1 + i, real_index_2 + i
          ]);
        }
      }

      // ...Then, add the bonds
      if (remove_duplicates) {
        for (const atom in local_bonds) {
          for (const linked of local_bonds[atom]) {
            bounds.push([
              Number(atom), linked
            ]);
          }
        }
      }

      // Increment i by number of atoms
      i += all_atom_count;
    }

    return bounds;
  } 
}();
