import { TopFile } from 'itp-parser';
import fs, { promises as FsPromise } from 'fs';
import path from 'path';
import TmpDirHelper from '../TmpDirHelper';
import { LIPIDS_ROOT_DIR, INSANE_HACK_SCRIPT } from '../constants';
import { Martinizer } from './Martinizer';
import RadiusDatabase from '../Entities/RadiusDatabase';
import logger from '../logger';
import { Database } from '../Entities/CouchHelper';
import MoleculeOrganizer from '../MoleculeOrganizer';
import { ArrayValues } from '../helpers';
import { Lipid } from '../Entities/entities';
import ShellManager, { JobInputs, JMError } from './ShellManager';
import Errors, { ErrorType } from '../Errors';
import ItpFile from 'itp-parser-forked';

export const AvailablePbcStrings = ['hexagonal', 'rectangular', 'square', 'cubic', 'optimal', 'keep'] as const;
export type PbcString = ArrayValues<typeof AvailablePbcStrings>;

export const AvailableRotateTypes = ['random', 'princ', 'angle'] as const;
export type RotateString = ArrayValues<typeof AvailableRotateTypes>;

export interface InsaneSettings {
  pbc: PbcString;
  /** Box size: Must be an array of 3, 6 or 9 integers. */
  box: number[];
  area_per_lipid?: number;
  area_per_lipid_upper?: number;
  random_kick_size?: number;
  bead_distance?: number;
  center?: boolean;
  orient?: boolean;
  rotate?: RotateString;
  rotate_angle?: number;
  grid_spacing?: number;
  hydrophobic_ratio?: number;
  fudge?: number;
  shift_protein?: number;
  charge?: number;
  salt_concentration?: number;
  solvent_type?: string;
}

// @ts-ignore
const InsaneParamToCliArg: { [T in keyof InsaneSettings]: string } = {
  pbc: '-pbc',
  area_per_lipid: '-a',
  area_per_lipid_upper: '-au',
  random_kick_size: '-rand',
  bead_distance: '-bd',
  center: '-center',
  orient: '-orient',
  grid_spacing: '-od',
  hydrophobic_ratio: '-op',
  fudge: '-fudge',
  shift_protein: '-dm',
  charge: '-charge',
  salt_concentration: '-salt',
  solvent_type: '-sol',
};

export interface InsaneRunnerOptions {
  force_field: string, 
  molecule_pdb?: string, 
  molecule_top?: string,
  molecule_itps?: string[], 
  lipids?: LipidMap, 
  upper_leaflet?: LipidMap,
  settings?: Partial<InsaneSettings>,
}

type SimpleLipidMap = [string, number][];
export type LipidMap = SimpleLipidMap | string[];

export const MembraneBuilder = new class MembraneBuilder {
  readonly SUPPORTED_LIPIDS: { [prefix: string]: string[] } = {};

  constructor() {
    // Init the supported lipids
    const supported_prefixes = new Set(Object.values(RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION));

    for (const dir of fs.readdirSync(LIPIDS_ROOT_DIR)) {
      if (!supported_prefixes.has(dir)) {
        continue;
      }

      if (!(dir in this.SUPPORTED_LIPIDS)) {
        this.SUPPORTED_LIPIDS[dir] = [];
      }

      for (const itp of fs.readdirSync(LIPIDS_ROOT_DIR + dir)) {
        if (!itp.endsWith('.itp')) {
          continue;
        }

        const name = itp.slice(0, itp.indexOf('.itp'));
        this.SUPPORTED_LIPIDS[dir].push(name);
      }
    }
  }


  /**
   * Build a membrane using INSANE.
   * 
   * Steps:
   * - With the given molecule (PDB+ITPs), selected force field and the defined lipids in parameter, generate a GRO box with INSANE
   * - Include the required ITPs (lipid-specific + force field + xxx.itp of the given molecule) in a newly created topology file
   * - Generate a PDB with CONECT entries using `Martinizer.createPdbWithConect()` function
   * - Returns TOP, PDB and ITPs file locations
   */
  async run({ force_field, molecule_pdb, molecule_top, molecule_itps, lipids, upper_leaflet = [], settings = {} }: InsaneRunnerOptions) {
    let ff_location = RadiusDatabase.FORCE_FIELD_TO_FILE_NAME[force_field];
    const ff_prefix = RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field];
    if (!ff_location || !ff_prefix) {
      throw new Error("Unknown force field. Please select a good force field");
    }
    else if (typeof ff_location === 'string') {
      ff_location = [ff_location];
    }

    let lipid_param: SimpleLipidMap;
    let upper_lipid_param: SimpleLipidMap;

    if (lipids){
      if (!lipids.length) {
        throw new Error("You need at least one lipid to insert.");
      }
      
      // If string[], convert to [string, 1][]
      lipid_param = typeof lipids[0] === 'string' ? 
        (lipids as string[]).map(e => [e, 1]) : 
        lipids as SimpleLipidMap;
  
      // Same for upper leaflet
      upper_lipid_param = upper_leaflet.length && typeof upper_leaflet[0] === 'string' ? 
        (upper_leaflet as string[]).map(e => [e, 1]) : 
        upper_leaflet as SimpleLipidMap;
  
      /// WITH DATABASE
      // Download every lipid
      // const lipids_entities = await Database.lipid.getAndThrowIfMissing([...lipid_param, ...upper_lipid_param].map(e => e[0]), force_field);
      
      /// WITH FILES
      // Check if every lipid is supported
      if (
        !lipid_param.every(l => this.SUPPORTED_LIPIDS[ff_prefix].includes(l[0])) ||
        !upper_lipid_param.every(l => this.SUPPORTED_LIPIDS[ff_prefix].includes(l[0]))
      ) {
        // unsupported lipid
        throw new Error("Unsupported lipid.");
      }
    }


    // Get a tmp dir
    const workdir = await TmpDirHelper.get();

    logger.debug(`[INSANE] Starting a INSANE run in directory ${workdir}.`);

    // Register all options
    const options: InsaneSettings = Object.assign({
      pbc: 'square',
      box: [7, 7, 9],
    }, settings);

    logger.debug(`[INSANE] Options: ` + JSON.stringify(options, null, 2));

    // Check if box size is valid
    if (![3, 6, 9].includes(options.box.length)) {
      throw new Error("Box has unsupported number of dimensions.");
    }

    let command_line = '';
    // Build the command line with fixed options (or that requires treatment)
    if (molecule_pdb !== "") {
      command_line = `-f "${molecule_pdb}" `;
    }
      // Output files (system and topology) 
      command_line += "-o system.gro -p __insane.top " +
      // Box size
      `-box ${options.box.map(e => Math.trunc(e)).join(',')} `;
      // Add all the lipids
      if (lipids) {
        // Lower leaflet/both leaflets if -u is missing
        //@ts-ignore
        command_line += `-l ${lipid_param.map(e => `${e[0]}:${e[1]}`).join(' -l ')} ` +
        // Upper leaflet (if defined)
        //@ts-ignore
        (upper_lipid_param.length ? `-u ${upper_lipid_param.map(e => `${e[0]}:${e[1]}`).join(' -u ')} ` : "");
      }
      
    
    if (options.rotate) {
      if (options.rotate === 'angle') {
        command_line += `-rotate ${options.rotate_angle} `;
      }
      else if (AvailableRotateTypes.includes(options.rotate)) {
        command_line += `-rotate ${options.rotate} `;
      }
    }
    
    // Add every supported item
    for (const opt in options) {
      const o = opt as keyof InsaneSettings;
      if (!(opt in InsaneParamToCliArg) || options[o] === undefined) {
        // Unsupported or invalid option
        continue;
      }

      if (typeof options[o] === 'boolean') {
        command_line += `${InsaneParamToCliArg[o]} `;
      }
      else {
        command_line += `${InsaneParamToCliArg[o]} ${options[o]} `;
      }
    }

    logger.debug("[INSANE] Command line: " + command_line);
    logger.debug(`[INSANE] Running INSANE with given settings.`);

    let jobOpt:JobInputs = { 
      "exportVar" : {
          "basedir" : workdir,
          "insaneArgs" : command_line,
          "insaneHackBefore" : INSANE_HACK_SCRIPT.BEFORE,
          "insaneHackAfter" : INSANE_HACK_SCRIPT.AFTER, 
          "inputFile" : molecule_pdb as string
      },
      "inputs" : {}
    };   

    // Start insane
    try {
      await ShellManager.run('insane', ShellManager.mode == "jm" ? jobOpt : `${INSANE_HACK_SCRIPT.BEFORE} ${INSANE_HACK_SCRIPT.AFTER} ${molecule_pdb} ${command_line}`, workdir, "insane");
    } catch (e) {
      // Handle error and throw the right error
      console.error("ShellManager.run crash"); 
      console.error(e.stack); 
      if (e instanceof JMError) return Errors.throw(ErrorType.JMError, {error: e.message})
      throw new InsaneError('insane_crash', workdir, 'error' in e ? e.error.stack : e.stack);
    }

    // Create the new TOP file

    /*
    > merge includes of a full top file create with createTopFile() system+molecules and molecules from insane.top (except "Protein")
    > Inject ITPs of lipids (todo...) they should be available in server
    File with modifications: __prepared.top

    Compile with gromacs script
    ~/Prog/martinize-db/utils/create_conect_pdb.sh system.gro __prepared.top "/Users/alki/Prog/martinize-db/utils/run.mdp" --remove-water
    */

    
    logger.debug(`[INSANE] Creating TOP file.`);
    //If we have a _rubber_band.itp, don't include it in the top file because it's included in molecule_x.itp. Also need to have a molecule_x.itp without this included itp to compute pdb without elastic bonds in conect fields.

    const itps_for_top = molecule_itps?.filter(itp => !itp.includes("_rubber_band"))

    let itps_without_elastic: string[] = []; 
    if(itps_for_top && (itps_for_top?.length != molecule_itps?.length)) { //We have elastic bonds, so filter "#include" statements and write new *_without_elastic.itp
      for (const itp of itps_for_top){
        const readedItp = await ItpFile.read(itp); 
        const bonds = readedItp.getField("bonds")
        const bonds_without_included_elastic = bonds.filter(line => !(line.startsWith("#include") && line.includes("_rubber_band")))
        if (bonds.length !== bonds_without_included_elastic.length) logger.warn("[INSANE] It seems *_rubber_band.itp exists but is not included inside molecule_*.itp")
        readedItp.setField("bonds", bonds_without_included_elastic)
        const itpPath = workdir + "/" + readedItp.type + "_without_elastic.itp"
        fs.writeFileSync(itpPath, readedItp.toString())
        itps_without_elastic.push(itpPath); 
      }
    }

    let wo_elastic_top = undefined; 
    try {
      var { top: full_top } = await Martinizer.createTopFile(
        workdir,
        molecule_top,
        itps_for_top,
        force_field
      );
      if(itps_without_elastic.length > 0) {
        const { top } = await Martinizer.createTopFile(workdir, molecule_top, itps_without_elastic, force_field)
        wo_elastic_top = top; 
      }

    } catch (e) {
      throw new InsaneError('top_file_crash', workdir, e.stack);
    }



    
    logger.debug(`[INSANE] Reading built TOP file and INSANE generate TOP file.`);
    const insane_top = await TopFile.read(workdir + "/__insane.top");

    const molecule_full_top = await TopFile.read(full_top);
    const readed_wo_elastic_top = wo_elastic_top ? await TopFile.read(wo_elastic_top) : undefined


    if (lipids) {
      // Compile the top files
      // Includes are normally all resolved in molecule_full_top (with the force field !)
      // We need to includes also the lipids ITPs
      //@ts-ignore
      const lipids_itp_names = this.getUniqueLipids(lipid_param, upper_lipid_param).map(e => e + ".itp");

      // Add the includes at the end of headlines
      molecule_full_top.headlines.push(...lipids_itp_names.map(e => `#include "${e}"`));
      if(readed_wo_elastic_top) readed_wo_elastic_top.headlines.push(...lipids_itp_names.map(e => `#include "${e}"`))
    }

    // Compile the top files together
    logger.debug(`[INSANE] Writing prepared TOP file.`);
    
    const prepared_top = await this.writePreparedTopFile(workdir + "/__prepared.top", insane_top, molecule_full_top);
    const prepared_top_wo_elastic = readed_wo_elastic_top ? await this.writePreparedTopFile(workdir + "/__prepared-no-elastic.top", insane_top, readed_wo_elastic_top) : undefined

    // Create lipids ITP files in working dir.
    // FF(s) symlink has been created by createTopFile() method.
    // + symlink of the molecule ITPs (needed)
    logger.debug(`[INSANE] Creating files for lipids ITPs.`);

    if (lipids) {
      /// WITH DATABASE
      // await this.createLipidItpFiles(workdir, lipids_entities);
      /// WITH FILES
      //@ts-ignore
      await this.createLipidItpSymlinks(workdir, force_field, lipid_param, upper_lipid_param);
    }
    if (molecule_itps !== undefined) {
      for (const itp of new Set(molecule_itps)) {
        await FsPromise.symlink(itp, workdir + "/" + path.basename(itp));
      }
    }
    

    // Ok, all should be ready. Start gromacs!
    logger.debug(`[INSANE] Creating the CONECT-ed PDB with GROMACS.`);
    try {
      const to_use_top = prepared_top_wo_elastic ? prepared_top_wo_elastic : prepared_top
      const to_use_gro = molecule_pdb ? "system-insane-hack.gro" : "system.gro"
      var pdbs = await Martinizer.createPdbWithConectWithoutWater(workdir + "/" + to_use_gro, to_use_top, workdir, lipids);
    } catch (e) {
      throw new InsaneError('gromacs_crash', workdir, e.stack);
    }

    // Build the ITP list without the force field ones
    const itps_ff = RadiusDatabase.getFilesForForceField(force_field);

    logger.debug(`[INSANE] Building ITP list.`);
    const itps = [] as string[];
    for (const file of await FsPromise.readdir(workdir)) {
      if (itps_ff.includes(file)) {
        continue;
      }
      if (!file.endsWith('.itp')) {
        continue;
      }

      itps.push(workdir + "/" + file);
    }

    logger.debug(`[INSANE] Run seems to be ok :)`);

    return {
      itps,
      pdbs,
      top: prepared_top,
    };
  }

  async prepareRunWithDatabaseMolecule(id: string) {
    const molecule = await Database.molecule.get(id);

    // Extract ZIP file in a temporary directory
    const tmp_dir = await TmpDirHelper.get();

    const { pdb, top, itps } = await MoleculeOrganizer.extract(molecule.files, tmp_dir);

    return { 
      force_field: molecule.force_field, 
      top,
      pdb,
      itps, 
    };
  }

  protected getUniqueLipids(lower: SimpleLipidMap, upper: SimpleLipidMap) {
    const lipids = new Set<string>();

    for (const e of [...lower, ...upper]) {
      lipids.add(e[0]);
    }

    return [...lipids];
  }

  protected async writePreparedTopFile(filename: string, insane: TopFile, protein: TopFile) {
    const stream = fs.createWriteStream(filename);

    try {
      const headlines = protein.headlines;
      const system = protein.getField('system');
      const molecules_prot = protein.getField('molecules');
      const molecules_insane = insane.getField('molecules');

      stream.write(headlines.join('\n') + '\n\n');

      stream.write('[system]\n');
      stream.write(system.join('\n') + '\n\n');

      stream.write('[molecules]\n');
      stream.write(molecules_prot.join('\n') + '\n');
      // Insert everything except the protein one
      stream.write(molecules_insane.filter(e => !e.startsWith('Protein ') && !e.startsWith(';')).join('\n') + '\n');
    } finally {
      stream.close();
    }

    return filename;
  }
  
  protected async createLipidItpFiles(workdir: string, lipids: Lipid[]) {
    for (const lipid of lipids) {
      await FsPromise.writeFile(workdir + "/" + lipid.name + ".itp", lipid.itp);
    }
  }

  protected async createLipidItpSymlinks(workdir: string, force_field: string, lower: SimpleLipidMap, upper: SimpleLipidMap) {
    if (!RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field]) {
      throw new Error("Lipid ITP directory not found.")
    } 

    const lipid_itp_base_dir = LIPIDS_ROOT_DIR + RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[force_field] + "/";
    
    for (const lipid of new Set([...lower.map(e => e[0]), ...upper.map(e => e[0])])) {
      await FsPromise.symlink(lipid_itp_base_dir + lipid + ".itp", workdir + "/" + lipid + ".itp");
    }
  }
};

export class InsaneError extends Error {
  constructor(
    public message: 'insane_crash' | 'gromacs_crash' | 'top_file_crash', 
    public workdir: string, 
    public trace: string
  ) { super(message); }
}

export default MembraneBuilder;
