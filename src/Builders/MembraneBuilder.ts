import RadiusDatabase from '../Entities/RadiusDatabase';
import { promisify } from 'util';
import { exec } from 'child_process';
import TmpDirHelper from '../TmpDirHelper/TmpDirHelper';
import { PYTHON2_PATH, INSANE_PATH, LIPIDS_ROOT_DIR } from '../constants';
import { Martinizer } from './Martinizer';
import { TopFile } from '../ItpParser';
import fs, { promises as FsPromise } from 'fs';
import path from 'path';
import logger from '../logger';
import { Database } from '../Entities/CouchHelper';
import MoleculeOrganizer from '../MoleculeOrganizer';
import { ArrayValues } from '../helpers';

const ExecPromise = promisify(exec);

export const AvailablePbcStrings = ['hexagonal', 'rectangular', 'square', 'cubic', 'optimal', 'keep'] as const;
export type PbcString = ArrayValues<typeof AvailablePbcStrings>;


export interface InsaneSettings {
  pbc: PbcString;
  /** Box size: Must be an array of 3, 6 or 9 integers. */
  box: number[];
}

export interface InsaneRunnerOptions {
  force_field: string, 
  molecule_pdb: string, 
  molecule_top: string,
  molecule_itps: string[], 
  lipids: LipidMap, 
  upper_leaflet?: LipidMap,
  settings?: Partial<InsaneSettings>,
}

type SimpleLipidMap = [string, number][];
export type LipidMap = SimpleLipidMap | string[];

export const MembraneBuilder = new class MembraneBuilder {
  /*
   * todo
   * Supported lipids for INSANE.
   * 
   * ITP name must be {lipid_name}.itp
   * If ITP name == {itp_name}, ITP file must be under {server_root}/lipids/{martini_version}/{itp_name}
   * 
   * {martini_version} is automatically determined with force field name (see RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION).
   */
  readonly SUPPORTED_LIPIDS = [
    'DPPC',
    'DLPC',
  ];

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
    if (!ff_location) {
      throw new Error("Unknown force field. Please select a good force field");
    }
    else if (typeof ff_location === 'string') {
      ff_location = [ff_location];
    }

    if (!lipids.length) {
      throw new Error("You need at least one lipid to insert.");
    }
    
    // If string[], convert to [string, 1][]
    const lipid_param: SimpleLipidMap = typeof lipids[0] === 'string' ? 
      (lipids as string[]).map(e => [e, 1]) : 
      lipids as SimpleLipidMap;

    // Same for upper leaflet
    const upper_lipid_param: SimpleLipidMap = upper_leaflet.length && typeof upper_leaflet[0] === 'string' ? 
      (upper_leaflet as string[]).map(e => [e, 1]) : 
      upper_leaflet as SimpleLipidMap;

    // Check if every lipid is supported
    if (
      !lipid_param.every(l => this.SUPPORTED_LIPIDS.includes(l[0])) ||
      !upper_lipid_param.every(l => this.SUPPORTED_LIPIDS.includes(l[0]))
    ) {
      // unsupported lipid
      throw new Error("Unsupported lipid.");
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

    // Build the command line
    const command_line = `${PYTHON2_PATH} "${INSANE_PATH}" ` +
      // Input PDB file
      `-f "${molecule_pdb}" ` +
      // Output files (system and topology) 
      "-o system.gro -p __insane.top " +
      // Box type
      `-pbc ${options.pbc} ` +
      // Box size
      `-box ${options.box.map(e => Math.trunc(e)).join(',')} ` +
      // Solvant and insane fixed settings
      `-center -sol W ` + 
      // Add all the lipids
      // Lower leaflet/both leaflets if -u is missing
      `-l ${lipid_param.map(e => `${e[0]}:${e[1]}`).join(' -l ')} ` +
      // Upper leaflet (if defined)
      (upper_lipid_param.length ? `-u ${upper_lipid_param.map(e => `${e[0]}:${e[1]}`).join(' -u ')} ` : "") +
      // More options TODO...
      ``;

    logger.debug(`[INSANE] Running INSANE with given settings.`);

    const stdout_insane = fs.createWriteStream(workdir + '/insane.stdout');
    const stderr_insane = fs.createWriteStream(workdir + '/insane.stderr');

    // Start insane (todo catch errors)
    try {
      const process_promise = ExecPromise(
        command_line, 
        { cwd: workdir }
      );
  
      process_promise.child.stdout?.pipe(stdout_insane);
      process_promise.child.stderr?.pipe(stderr_insane);

      await process_promise;
    } finally {
      stdout_insane.close();
      stderr_insane.close();
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
    const { top: full_top } = await Martinizer.createTopFile(
      workdir,
      molecule_top,
      molecule_itps,
      force_field
    );

    logger.debug(`[INSANE] Reading built TOP file and INSANE generate TOP file.`);
    const insane_top = new TopFile(workdir + "/__insane.top", []);
    await insane_top.read();

    const molecule_full_top = new TopFile(full_top, []);
    await molecule_full_top.read();

    // Compile the top files
    // Includes are normally all resolved in molecule_full_top (with the force field !)
    // We need to includes also the lipids ITPs
    const lipids_itp_names = this.getUniqueLipids(lipid_param, upper_lipid_param).map(e => e + ".itp");

    // Add the includes at the end of headlines
    molecule_full_top.headlines.push(...lipids_itp_names.map(e => `#include "${e}"`));

    // Compile the top files together
    logger.debug(`[INSANE] Writing prepared TOP file.`);
    const prepared_top = await this.writePreparedTopFile(workdir + "/__prepared.top", insane_top, molecule_full_top);

    // Create symlinks of the lipids ITP in working dir.
    // FF(s) symlink has been created by createTopFile() method.
    // + symlink of the molecule ITPs (needed)
    logger.debug(`[INSANE] Creating symlinks for lipids ITPs.`);
    await this.createLipidItpSymlinks(workdir, force_field, lipid_param, upper_lipid_param);
    for (const itp of new Set(molecule_itps)) {
      await FsPromise.symlink(itp, workdir + "/" + path.basename(itp));
    }

    // Ok, all should be ready. Start gromacs!
    logger.debug(`[INSANE] Creating the CONECT-ed PDB with GROMACS.`);
    const pdbs = await Martinizer.createPdbWithConectWithoutWater(workdir + "/system.gro", prepared_top, workdir);

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

export default MembraneBuilder;
