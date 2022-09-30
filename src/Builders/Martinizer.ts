import axios, { AxiosResponse } from 'axios';
import { ExecException } from 'child_process';
import FormData from 'form-data';
import fs, { exists, promises as FsPromise } from 'fs';
import path from 'path';
import readline from 'readline';
import TarStream from 'tar-stream';
import zlib from 'zlib';
import { FORCE_FIELD_DIR, CONECT_MDP_PATH, CREATE_MAP_PY_SCRIPT_PATH, CREATE_GO_PY_SCRIPT_PATH, DSSP_PATH } from '../constants';
import RadiusDatabase from '../Entities/RadiusDatabase';
import Errors, { ErrorType } from '../Errors';
import { ArrayValues, fileExists } from '../helpers';
import logger from '../logger';
import TmpDirHelper from '../TmpDirHelper';
import { TopFile, ItpFile } from 'itp-parser-forked';
import JSZip from 'jszip';
import ShellManager, { JobInputs, JMError } from './ShellManager';
import JMSurcouche from './JMSurcouche';
import { pathsToInputs, str_to_stream } from './JMSurcouche';
import { cpuUsage } from 'process';
import { Readable } from 'stream';

/**
 * Tuple of two integers: [{from} atom index, {to} atom index]
 */
export type ElasticOrGoBounds = [number, number];

export interface GoBoundsDetails {
  index_to_real: { [index: number]: number };
  name_to_index: { [name: string]: number };
  index_to_name: { [index: number]: string };
  real_to_index: { [index: number]: number };
  /** Atom count */
  count: number;
}

export type GoMoleculeDetails = { [moleculeType: string]: GoBoundsDetails };

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


const MARTINIZE_POSITIONS = ['none', 'all', 'backbone'] as const;
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
  use_go?: boolean;

  /** Set neutral termini */
  neutral_termini?: boolean;
  /** Apply side chains corrections */
  side_chain_fix?: boolean;
  /** Cystein bounds */
  cystein_bridge?: string;

  cter?: string;
  nter?: string;
  commandline: string;
  advanced?: boolean;
  builder_mode?: "elastic" | "go" | "classic"
}

export const Martinizer = new class Martinizer {
  STEP_MARTINIZE_INIT = 'init';
  STEP_MARTINIZE_RUNNING = 'internal';
  STEP_MARTINIZE_ENDED_FINE = 'martinize-end';
  STEP_MARTINIZE_GET_CONTACTS = 'contacts';
  STEP_MARTINIZE_GO_SITES = 'go-sites';
  STEP_MARTINIZE_GROMACS = 'gromacs';


  protected MAX_JOB_EXECUTION_TIME = 5 * 60 * 1000;

  isMartinizePosition(value: any): value is MartinizePosition {
    return MARTINIZE_POSITIONS.includes(value);
  }

  settingsToCommandline(settings: Partial<MartinizeSettings>) {
    const full: MartinizeSettings = Object.assign({}, {
      input: '',
      ff: 'martini22',
      position: 'none',
      commandline: ''
    }, settings);

    /*
    if (!full.input.trim()) {
      throw new Error("Invalid input file");
    }
    */

    const basename = full.input;
    const with_ext = basename + '.pdb';

    // Check dssp ps
    // TODO: DSSP gives bad results... this should not append
    let command_line = " -f " + with_ext + " -x output.pdb -o system.top -ff " + full.ff + " -p " + full.position
    if (DSSP_PATH) command_line += " -dssp " + DSSP_PATH
    // let command_line = "martinize2 -f " + with_ext + " -x output.pdb -o system.top -dssp " + DSSP_PATH + " -ff " + full.ff + " -p " + full.position + " ";

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
    if (full.use_go) {
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
    if (full.cter) {
      command_line += " -cter " + full.cter
    }
    if (full.nter) {
      command_line += " -nter " + full.nter
    }

    return { command_line: command_line, basename: basename, with_ext: with_ext, full: full }
  }

  /**
   * Create a martinize run.
   * Returns created path to created PDB, TOP and ITP files.
   */
  async run(settings: Partial<MartinizeSettings>, onStep?: (step: string, ...data: any[]) => void, path?: string) {

    const MARTINIZE_WARN = "martinize_warnings.log"
    const OUTPUT_PDB = "output.pdb"
    const OUTPUT_TOP = "system.top"

    // try-finally to ensure rename even if error happend

    // Step: Martinize Init
    onStep?.(this.STEP_MARTINIZE_INIT);
    let jobOpt: JobInputs = {
      exportVar: {
        MARTINIZE_WARN,
        OUTPUT_PDB,
        OUTPUT_TOP,
        DSSP_PATH: DSSP_PATH ?? "",
        FORCE_FIELD: settings.ff ?? "martini22",
        POSITION: settings.position ?? "none",
        CYSTEIN_BRIDGE: settings.cystein_bridge ?? "none",
        C_TER: settings.cter ?? "[]",
        N_TER: settings.nter ?? "[]",
        SIDE_CHAIN_FIX: settings.side_chain_fix ? "true" : ""
      },
      inputs: {
        "input.pdb": settings.input
      },
    };

    console.log("run", jobOpt)
    try {
      const { stdout, jobFS } = await JMSurcouche.run(
        'martinize', jobOpt
      );

      const pdbStream = await jobFS.readToStream(OUTPUT_PDB)
      console.log("output pdb", OUTPUT_PDB)
      console.log(pdbStream.readableLength);
      console.log(pdbStream.readable) 

      if (!pdbStream.readable) {
        return Errors.throw(ErrorType.MartinizeNoOutput, {error : "Pdb file is empty"})
      }

      const itp_files = await jobFS.list("*.itp")
      const itpContents : {[name : string] : Readable} = {}
      for (const itpName of itp_files){
        itpContents[itpName] = await jobFS.readToStream(itpName)
      }

      logger.debug(`[MARTINIZER-RUN] Generated PDB is found, run should be fine.`);
      onStep?.(this.STEP_MARTINIZE_ENDED_FINE);


      logger.debug("[MARTINIZER-RUN] Creating full TOP file for Martinize built molecule.");
      const topStream = await jobFS.readToStream(OUTPUT_TOP)
      if (!topStream.readable) {
        return Errors.throw(ErrorType.MartinizeNoOutput, { error: "Top file is empty" })
      }

      const top = await this.createTopFileFromStream(topStream, settings.ff!, itp_files )

      console.log("resulting top", top)

      logger.debug("[MARTINIZER-RUN] Creating PDB with CONECT entries for Martinize built molecule.");
      onStep?.(this.STEP_MARTINIZE_GROMACS);
      //TAKE ELASTIC TOP IF ELASTIC IS DONE
      const { pdb: pdb_with_conect, top: full_top } = await this.createPdbWithConectFromStream(pdbStream, "pdb", top, false, settings.ff, itpContents);

      logger.debug("[MARTINIZER-RUN] Run is complete, everything seems to be fine :)");
      console.log(top)

      const warnString = await jobFS.readToString(MARTINIZE_WARN)

      const sortedItps = settings.builder_mode === "elastic" || settings.ff?.includes("elnedyn") ? await this.sortItpsFromTop(full_top, itp_files) : [itp_files]

      return {
        pdb: pdb_with_conect,
        itps: sortedItps,
        top: full_top,
        warns: warnString,
        jobId: jobFS.job.id,
        //dir: dir,
        elastic_bonds: undefined,
      };




    } catch (e) {
      console.error(e)
      if (e instanceof JMError) return Errors.throw(ErrorType.JMError, { error: e.message })
      return Errors.throw(ErrorType.MartinizeRunFailed, {
        error: "Martinize has failed with an non-zero exit code.",
        type: "non-zero",
        //dir,
      });
    }

    // Scan for files in result dir












    // If go mode, we should compute map + run a python script to refresh ITPs files.
    // if (settings.use_go) {
    //   let map_filename: string;
    //   // GET THE MAP FILE FROM A CUSTOM WAY.
    //   // Use the original pdb file !!!
    //   // todo change (ccmap create way too much distances, so the shell script takes forever)
    //   onStep?.(this.STEP_MARTINIZE_GET_CONTACTS);
    //   try {
    //     map_filename = await this.getCcMapRCSU(with_ext, dir);
    //   } catch {
    //     return Errors.throw(ErrorType.MartinizeRunFailed, {
    //       error: "Unable to create contact map.",
    //       type: "contact-map",
    //       stdout: dir + '/rcsu.stdout',
    //       stderr: dir + '/rcsu.stderr',
    //       dir,
    //     });
    //   }

    //   const itp = await ItpFile.read(dir + "/molecule_0.itp")

    //   const firstResidueNumber: number = parseInt(itp.atoms[0].split(" ").filter(splitElmt => splitElmt !== '')[2])

    //   const nbAtomsWithoutGO = itp.atoms.filter(atomLine => {
    //     const splittedLine = atomLine.split(" ").filter(splitElmt => splitElmt !== "")
    //     if (splittedLine[1].includes("molecule")) return false
    //     else return true
    //   }).length

    //   logger.debug("[MARTINIZER-RUN] Creating Go virtual bonds");
    //   onStep?.(this.STEP_MARTINIZE_GO_SITES);

    //   // Ensure to have the script at the correct place...
    //   try {
    //     // Must create the go sites
    //     const moltype = "molecule_0"
    //     const goArgs = `-s output.pdb -f ${map_filename} --moltype ${moltype} --Natoms ${nbAtomsWithoutGO} --missres ${firstResidueNumber - 1}`

    //     const jobOptGo: JobInputs = {
    //       exportVar: {
    //         WORKDIR: dir,
    //         GO_ARGS: goArgs
    //       },
    //       inputs: {},
    //     };

    //     const command_line_go = `"${CREATE_GO_PY_SCRIPT_PATH} ${goArgs}"`;

    //     await ShellManager.run(
    //       'go_virt',
    //       ShellManager.mode === "jm" ? jobOptGo : command_line_go,
    //       dir,
    //       'go-virt-sites',
    //     );
    //   } catch (e) {
    //     if (e instanceof JMError) return Errors.throw(ErrorType.JMError, { error: e.message })
    //     return Errors.throw(ErrorType.MartinizeRunFailed, {
    //       error: "Unable to create go virtual sites.",
    //       type: "create-go-virt",
    //       stdout: dir + '/go-virt-sites.stdout',
    //       stderr: dir + '/go-virt-sites.stderr',
    //       dir,
    //     });
    //   }

    //   // Ok, ITP file refreshed.
    // }



    // let elasticBounds: ElasticOrGoBounds[] | undefined = undefined;
    // let elasticTop: string | undefined = undefined;
    // if (settings.builder_mode === "elastic" || settings.ff?.includes("elnedyn")) {
    //   console.log("elastic bounds", dir)
    //   const { elastic_bounds, elastic_itps, itp_without_elastic } = await this.computeElasticNetworkBounds(full_top, itp_files, dir);
    //   for (const itp of elastic_itps) {
    //     console.log("elastic itp", itp)
    //     itp_files.push(itp);
    //   }
    //   elasticBounds = elastic_bounds;
    //   //Create top file without elastic bounds to create pdb without elastic in CONECT fields
    //   const { top } = await this.createTopFile(dir, dir + '/system.top', itp_without_elastic, full.ff, "full_without_elastic.top");
    //   elasticTop = top;
    // }

    // // Generate the right PDB file (with conect entries)


    // logger.debug("[MARTINIZER-RUN] Sort itps")

    // //const sortedItps = settings.builder_mode === "elastic" || settings.ff?.includes("elnedyn") ? await this.sortItpsFromTop(full_top, itp_files) : [itp_files]

    // logger.debug("[MARTINIZER-RUN] Run is complete, everything seems to be fine :)");

    // return {
    //   pdb: pdb_with_conect,
    //   itps: sortedItps,
    //   top: full_top,
    //   warns: dir + "/" + martinizeWarnOut,
    //   dir: dir,
    //   elastic_bonds: elasticBounds,
    // };

  }

  async sortItpsFromTop(topFile: string, itp_files: string[]) {
    const sortedItps: string[][] = []
    const system = await TopFile.read(topFile);
    for (const mol of system.molecules) {
      const name = mol.type;
      const itps = itp_files.filter(itp_filename => itp_filename.includes(name))
      sortedItps.push(itps)
    }
    return sortedItps
  }

  /**
   * Modify the original top file in order to include the right ITPs, 
   * and create a link of the force field martini file into current directory.
   * 
   * ITPs in {itps_path} MUST be in {current_directory} !
   * 
   * Returns new TOP filename and all the used ITPs to generate top.
   */
  async createTopFile(current_directory: string, original_top_path: string | undefined, itps_path: string[] | undefined, force_field: string, top_name: string = "full.top") {
    let itps_ff = RadiusDatabase.FORCE_FIELD_TO_FILE_NAME[force_field];

    if (!itps_ff) {
      throw new ReferenceError("Your force field is invalid: can't find related ITPs.");
    }

    itps_ff = typeof itps_ff === 'string' ? [itps_ff] : itps_ff;

    let itps = undefined;
    if (itps_path !== undefined) {
      itps = [...itps_path, ...itps_ff.map(e => FORCE_FIELD_DIR + e)];
    }
    else {
      itps = [...itps_ff.map(e => FORCE_FIELD_DIR + e)];
    }
    const base_ff_itps = [] as string[];

    // Create everysym link
    for (const itp of itps_ff) {
      const itp_path = FORCE_FIELD_DIR + itp;
      const dest = path.resolve(current_directory + "/" + itp);
      base_ff_itps.push(itp);
      try {
        await FsPromise.symlink(itp_path, dest);
      } catch (e) {
        const err = e as any
        if (err.code === "EEXIST") logger.verbose(`${itp_path} symlink already exists`)
        else throw new Error(err);
      }

    }

    let real_itps = undefined;
    if (itps_path !== undefined) {
      real_itps = [...base_ff_itps, ...itps_path.map(e => path.basename(e))];
    }
    else {
      real_itps = [...base_ff_itps];
    }
    const top = current_directory + "/" + top_name;

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
    if (original_top_path !== "") {
      const top_read_stream = readline.createInterface({
        //@ts-ignore
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
    }
    else {
      top_write_stream.write(includes.join('\n') + '\n');
    }

    top_write_stream.close();


    return {
      top,
      itps,
    };

  }

  async createTopFileFromStream(originalTopStream: Readable, force_field: string, itps_path?: string[]) {
    console.log("create top from stream")
    const itps_ff = RadiusDatabase.getFilesForForceField(force_field)

    if (!itps_ff) {
      throw new ReferenceError("Your force field is invalid: can't find related ITPs.");
    }

    let real_itps = undefined;
    if (itps_path !== undefined) {
      real_itps = [...itps_ff, ...itps_path.map(e => path.basename(e))];
    }
    else {
      real_itps = [...itps_ff];
    }

    const includes: string[] = [];

    // Define the includes
    for (const itp of real_itps) {
      // Exclude the GO ITPs, they're already included in martini_304.itp
      if (itp.endsWith('VirtGoSites.itp') || itp.endsWith('go4view_harm.itp')) {
        continue;
      }

      includes.push(`#include "${itp}"`);
    }
    
    console.log("includes", includes)

    let newTop = ""

    if (originalTopStream.readable) {
      const topString = originalTopStream.read().toString()
      //console.log("TEST READ STREAM", test.toString())
      let includes_included = false;

      // Remove every #include line
      for await (const line of topString.split("\n")) {
        if (line.startsWith('#include')) {
          if (!includes_included) {
            // Include the hand-crafted includes
            newTop = newTop + includes.join('\n') + '\n'
            includes_included = true;
          }

          continue;
        }
        newTop = newTop + line + '\n'

      }
    }
    else {
      newTop = newTop + includes.join('\n') + '\n'
    }

    console.log("newTop", newTop)

    return newTop;

  }

  async createTopFileToString(original_top_path: string, itps_path: string[] | undefined, force_field: string) {
    console.log("create top file from string")
    const itps_ff = RadiusDatabase.getFilesForForceField(force_field)

    if (!itps_ff) {
      throw new ReferenceError("Your force field is invalid: can't find related ITPs.");
    }

    const base_ff_itps = [] as string[];

    let real_itps = undefined;
    if (itps_path !== undefined) {
      real_itps = [...itps_ff, ...itps_path.map(e => path.basename(e))];
    }
    else {
      real_itps = [...itps_ff];
    }

    //const top = current_directory + "/" + top_name;

    const includes: string[] = [];

    // Define the includes
    for (const itp of real_itps) {
      // Exclude the GO ITPs, they're already included in martini_304.itp
      if (itp.endsWith('VirtGoSites.itp') || itp.endsWith('go4view_harm.itp')) {
        continue;
      }

      includes.push(`#include "${itp}"`);
    }

    let newTop = ""

    //const top_write_stream = fs.createWriteStream(top);

    if (original_top_path !== "") {
      const top_read_stream = readline.createInterface({
        //@ts-ignore
        input: fs.createReadStream(original_top_path),
        crlfDelay: Infinity,
      });

      let includes_included = false;

      // Remove every #include line
      for await (const line of top_read_stream) {
        console.log("tl", line)
        if (line.startsWith('#include')) {
          if (!includes_included) {
            // Include the hand-crafted includes
            newTop = newTop + includes.join('\n') + '\n'
            includes_included = true;
          }

          continue;
        }
        newTop = newTop + line + '\n'

      }
    }
    else {
      newTop = newTop + includes.join('\n') + '\n'
    }

    console.log("newTop", newTop)

    return newTop;

  }

  async getCcMapRCSU(pdb_filename: string, use_tmp_dir?: string) {
    logger.debug("GET MAP RCSU")
    console.log(pdb_filename, use_tmp_dir)
    if (!use_tmp_dir) {
      use_tmp_dir = await TmpDirHelper.get();
      logger.debug(`Created tmp directory for rcsu: ${use_tmp_dir}.`);
    }
    const outputFile = "output.map"
    const jobOpt: JobInputs = {
      exportVar: {
        WORKDIR: use_tmp_dir,
        PDB: path.resolve(pdb_filename),
        OUTPUT: outputFile,
      },
      inputs: {},
    };

    try {
      await ShellManager.run(
        'map_rcsu',
        jobOpt,
        use_tmp_dir,
        'rcsu',
      );
    }
    catch (e) {
      if (e instanceof JMError) return Errors.throw(ErrorType.JMError, { error: e.message })
      throw new Error(e)
    }

    return outputFile

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

    const jobOpt: JobInputs = {
      exportVar: {
        WORKDIR: use_tmp_dir,
        INPUT_PDB: path.resolve(pdb_filename),
        DISTANCES: distances_file,
      },
      inputs: {},
    };

    const command_line = `${CREATE_MAP_PY_SCRIPT_PATH} "${path.resolve(pdb_filename)}" "${distances_file}"`

    // Compute contacts with the CA pdb
    try {
      await ShellManager.run(
        'ccmap',
        ShellManager.mode === "jm" ? jobOpt : command_line,
        use_tmp_dir,
        'distances',
      );
    }
    catch (e) {
      if (e instanceof JMError) return Errors.throw(ErrorType.JMError, { error: e.message })
    }


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
  async createPdbWithConect(pdb_or_gro_path: string, top_content: string, remove_water: boolean = false, force_field: string = "martini22", itps?: string[], lipids?: any) {
    // let groups_to_del = 17;
    // if (lipids) {
    //   groups_to_del += lipids.length * 2;
    // }

    console.log("CREATE PDB WITH CONECT POUET")

    const force_fields = RadiusDatabase.getCompleteFilesForForceField(force_field)
    const pdb_or_gro_name = path.basename(pdb_or_gro_path)
    console.log(force_fields)
    const command = {
      exportVar: {
        "DEL_WATER_BOOL": remove_water ? "YES" : "NO",
        "INPUT_NAME": pdb_or_gro_name
      },
      inputs: {
        ...pathsToInputs([pdb_or_gro_path]),
        'input.top': str_to_stream(top_content),
        'run.mdp': CONECT_MDP_PATH,
        ...pathsToInputs(force_fields)
      }
    };

    if (itps) {
      command.inputs = { ...command.inputs, ...pathsToInputs(itps) }
    }

    try {
      //@ts-ignore
      const { stdout, jobFS } = await JMSurcouche.run('conect', command)
      console.log("JOB end", jobFS.job.id)
      const pdb_stream = await jobFS.readToStream("output-conect.pdb")
      const top_stream = await jobFS.readToStream('input.top')
      if (!pdb_stream.readable) {
        throw new Error(`PDB could not be created for an unknown reason (or more than 1 output exists)`);
      }

      console.log(await jobFS.list())
      if (!top_stream.readable) {
        throw new Error(`Top could not be copied for an unknown reason (or more than 1 top exists)`);
      }
      return { pdb_stream, top_stream };

    } catch (e) {
      logger.error("Error while job")
      console.log(e)
      throw new Error(`Error while job`);

    }

  }

  async createPdbWithConectFromStream(inputStream : Readable, inputType : "pdb" | "gro", top_content: string, remove_water: boolean = false, force_field: string = "martini22", itps?: {[name: string] : Readable}, lipids?: any) {
    // let groups_to_del = 17;
    // if (lipids) {
    //   groups_to_del += lipids.length * 2;
    // }

    console.log("CREATE PDB WITH CONECT POUET")

    const force_fields = RadiusDatabase.getCompleteFilesForForceField(force_field)
    const inputName = "input." + inputType
    const command = {
      exportVar: {
        "DEL_WATER_BOOL": remove_water ? "YES" : "NO",
        "INPUT_NAME": inputName
      },
      inputs: {
        [inputName] : inputStream,
        'input.top': str_to_stream(top_content),
        'run.mdp': CONECT_MDP_PATH,
        ...pathsToInputs(force_fields)
      }
    };

    if (itps) {
      command.inputs = { ...command.inputs, ...itps }
    }

    console.log(command); 

    try {
      //@ts-ignore
      const { stdout, jobFS } = await JMSurcouche.run('conect', command)
      console.log("JOB end", jobFS.job.id)
      const pdb_out = await jobFS.list("output-conect.pdb")
      const top_out = await jobFS.list('input.top')
      if (pdb_out.length !== 1) {
        throw new Error(`PDB could not be created for an unknown reason (or more than 1 output exists)`);
      }

      if (top_out.length !== 1) {
        throw new Error(`Top could not be copied for an unknown reason (or more than 1 top exists)`);
      }
      return { pdb: pdb_out[0], top: top_out[0] };

    } catch (e) {
      logger.error("Error while job")
      console.log(e)
      throw new Error(`Error while job`);

    }

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
  async createPdbWithConectWithoutWater(pdb_or_gro_filename: string, top_filename: string, base_directory: string, lipids?: any) {
    const pdb_water = await this.createPdbWithConect(pdb_or_gro_filename, top_filename, true, lipids);

    const pdb_no_w = base_directory + "/output-conect-no-w.pdb";
    const exists = await FsPromise.access(pdb_no_w, fs.constants.F_OK).then(() => true).catch(() => false);

    if (!exists) {
      throw new Error("PDB could not be created for an unknown reason. Check the files gromacs.stdout and gromacs.stderr in directory " + base_directory + ".");
    }

    return {
      water: pdb_water,
      no_water: pdb_no_w
    };
  }

  /**
   * Generate an object than contain "additionnal" bounds between atoms, that are created with the elastic network (option -elastic).
   * 
   * Specificy only ITP files that are generated by Martinize, not the force field itself!!
   * 
   * TODO: worker thread
   */
  async computeElasticNetworkBounds(top_file: string, itp_files: string[], workdir: string) {
    logger.verbose("[ELASTIC-BUILD] Constructing elastic network bonds.");

    const bounds: ElasticOrGoBounds[] = [];

    logger.debug("[ELASTIC-BUILD] Reading TOP+ITP files.");
    const top = await TopFile.read(top_file, itp_files);

    logger.debug(`[ELASTIC-BUILD] Available molecules: ${top.molecules
      .map(molecule => `${molecule.type} (${molecule.count} time${molecule.count > 1 ? 's' : ''})`)
      .join(', ')
      }`);

    // Incrementer for designating PDB line
    let i = 0;
    let elasticItps = [];
    let withoutElasticItps = [];
    for (const molecule of top.molecules) {
      logger.debug("[ELASTIC-BUILD] Reading molecule " + molecule.type + ".");

      // Get the number of atoms in a single chain of this molecule
      const itp = molecule.itp;
      const atom_count = itp.atoms.filter(line => line && !line.startsWith(';')).length;

      //const name = molecule.name; 
      //Write elastic bonds in an other itp file to avoid elastic bonds representation with ngl. Output connect will be computed without this new file, and then it will be included again. 
      const elastic_bonds = itp.getSubfield("bonds", "Rubber band")
      const elastic_itp_name = molecule.type + "_rubber_band.itp";
      const elastic_itp_path = workdir + "/" + elastic_itp_name

      const elastic_itp = new ItpFile();
      elastic_itp.appendField("bonds", elastic_bonds)

      fs.writeFileSync(elastic_itp_path, elastic_itp.toString())
      elasticItps.push(elastic_itp_path)

      //Delete elastic bonds from current itp
      const correctedItp = workdir + "/" + molecule.type + "_without_elastic.itp";
      itp.removeSubfield("bonds", "Rubber band");
      fs.writeFileSync(correctedItp, itp.toString()) //Write itp without rubber bands
      withoutElasticItps.push(correctedItp)


      itp.appendInclude(elastic_itp_name, "bonds");
      fs.writeFileSync(workdir + "/" + molecule.type + ".itp", itp.toString()) //Rewrite initial itp with include statement for rubber bands


      for (const band of elastic_bonds) {
        if (!band.startsWith(";")) {
          const [atom_from, atom_to,] = band.split(/\s+/g);
          bounds.push([
            Number(atom_from) + i,
            Number(atom_to) + i,
          ]);
        }
      }

      i += atom_count;

    }
    return {
      elastic_bounds: bounds,
      elastic_itps: elasticItps,
      itp_without_elastic: withoutElasticItps
    }

    /*for (const molecule of top.molecules) {
      logger.debug("[ELASTIC-BUILD] Reading molecule " + molecule.type + ".");

      // Get the number of atoms in a single chain of this molecule
      const itp = molecule.itp;

      const atom_count = itp.atoms.filter(line => line && !line.startsWith(';')).length;
      let chain_n = 0;

      // There is one ITP per chain
      for (let mol_count = 0; mol_count < molecule.count; mol_count++) {
        let should_read = false;
        chain_n++;
        
        logger.verbose(`[ELASTIC-BUILD] Chain ${chain_n} of molecule "${molecule.type}": ${itp.bonds.length} bond lines.`);

        for (const band of itp.bonds) {
          if (band.startsWith(';')) {
            // Find the elastic related comments
            if (
              // Todo verify if it is elastic bounds
              // band.startsWith('; Long elastic bonds for extended regions') ||
              // band.startsWith('; Short elastic bonds for extended regions') ||
              band.startsWith('; Rubber band')
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
    }*/


  }

  /**
   * Generate an object than contain "additionnal" bounds between atoms, that are created with the Go model (option -govs-include).
   * 
   * Specificy only ITP files that are generated by Martinize, not the force field itself!!
   * 
   * TODO: worker thread
   */
  async __UNSAFEcomputeGoModelBounds(top_file: string, itp_files: string[], remove_duplicates = true) {
    /*
     *  Même si il y a deux chaînes, en mode go elles seront dans une seule molécule (normalement).
     *
     *  1: Lire les noms des atoms Go dans [atoms]: {molecule_type}_{i}
     *  2: Lire les associations go_atom => real atom index dans [virtual_sitesn] après comment "; Virtual go site"
     *  3: Lire le fichier {molecule_type}_go-table-VirtGoSites.itp qui définit les liaisons go_atom <=> go_atom
     *  4: Convertir ces liaisons go <=> go en "additionnal" bonds entre real atoms
     */

    logger.verbose("[GO-VIRT-SITES] Reading system topology.");
    const top = await TopFile.read(top_file);
    const molecule_types: string[] = [];

    for (const molecule of top.getField('molecules', true)) {
      // molecule_0    1 (count is always 1 for Go model.)
      const [name,] = molecule.split(ItpFile.BLANK_REGEX);
      molecule_types.push(name);
    }

    // Here, we store bounds
    const bounds: ElasticOrGoBounds[] = [];
    const details: GoMoleculeDetails = {};

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
      const molecule_itp = await ItpFile.read(itp_files[molecule_itp_index]);
      const go_table = await ItpFile.read(itp_files[go_table_index]);

      // Read the molecule file
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Reading ITP files.`);

      // Count all atoms, used to increment atom counter at the end of loop
      const all_atom_count = molecule_itp.atoms.filter(line => line && !line.startsWith(';')).length;

      const prefix = molecule_type + '_';
      /** Link go atom index to real atom index. */
      const index_to_real: { [index: number]: number } = {};
      /** Link go atom name to go atom index. */
      const name_to_index: { [name: string]: number } = {};

      // WILL BE USEFUL WHEN BOUNDS ARE INSERTED/DELETED DYNAMICALLY
      /** Link go atom index to go atom name. */
      const index_to_name: { [index: number]: string } = {};
      /** Link real atom index to go atom index. */
      const real_to_index: { [index: number]: number } = {};

      details[molecule_type] = { index_to_real, name_to_index, index_to_name, real_to_index, count: all_atom_count };

      // Step 1: Find atoms that name start by "{molecule_type}_" in category "atoms"
      logger.debug(`[GO-VIRT-SITES] [${molecule_type}] Looking for virtual atoms.`);

      for (const atom_line of molecule_itp.atoms) {
        // Typical line is : 
        // 2575 molecule_0_9       9 LYS CA  2575    0 

        const [index, name,] = atom_line.split(ItpFile.BLANK_REGEX);

        if (name.startsWith(prefix)) {
          name_to_index[name] = Number(index);
          index_to_name[Number(index)] = name;
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
        real_to_index[Number(real_index)] = Number(go_index);
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
        const [name1, name2,] = line.split(ItpFile.BLANK_REGEX).filter(e => e);

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


    return { bounds, details };
  }

  protected async zip(dir: string) {
    const zip = new JSZip;

    for (const file of await FsPromise.readdir(dir)) {
      const stat = await FsPromise.stat(dir + "/" + file);

      // 10 MB max
      if (stat.size < 10 * 1024 * 1024 && !stat.isSymbolicLink() && stat.isFile()) {
        const name = file.endsWith('.stderr') || file.endsWith('.stdout') ? file + '.txt' : file;
        zip.file(name, await FsPromise.readFile(dir + "/" + file));
      }
    }

    return zip;
  }

  /**
   * Zip a directory.
   * 
   * Todo: worker thread
   * @param dir 
   */
  async zipDirectory(dir: string) {
    return (await this.zip(dir)).generateAsync({
      compression: "DEFLATE",
      compressionOptions: { level: 6 },
      type: "arraybuffer",
    });
  }

  async zipDirectoryString(dir: string) {
    return (await this.zip(dir)).generateAsync({
      compression: "DEFLATE",
      compressionOptions: { level: 6 },
      type: 'array',
    });
  }

}();
