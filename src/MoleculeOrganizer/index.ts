import fs, { promises as FsPromise } from 'fs';
import { MOLECULE_ROOT_DIR } from '../constants';
import logger from '../logger';
import { generateSnowflake, basenameWithoutExt } from '../helpers';
import JSZip from 'jszip';
import md5File from 'md5-file/promise';
import { Martinizer } from '../Builders/Martinizer';
import path from 'path';
import TmpDirHelper from '../TmpDirHelper';
// @ts-ignore 
import NodeStreamZip from 'node-stream-zip';
import Errors, { ErrorType } from '../Errors';
import { ConsoleTransportOptions } from 'winston/lib/winston/transports';
import { SimuFile } from '../routes/molecule/CreateMoleculeJson';

export const MoleculeOrganizer = new class MoleculeOrganizer {


  constructor() {

    try {
      fs.mkdirSync(MOLECULE_ROOT_DIR);
    } catch { }

  }

  /**
   * Get a ZIP, read by JSZip.
   * 
   * If the save doesn't exists, returns `undefined`.
   */
  async get(file_id: string): Promise<[JSZip, MoleculeSaveInfo] | undefined> {
    if (await this.exists(file_id)) {
      const zip_buffer = await FsPromise.readFile(this.getFilenameFor(file_id));
      const zip = await JSZip.loadAsync(zip_buffer);

      return [
        zip,
        (await this.getInfo(file_id))!
      ];
    }
  }

  async extract(id: string, in_directory: string) {
    const infos = await this.getInfo(id);
    if (!infos) {
      throw new Error("Save not found.");
    }

    const zip = new NodeStreamZip({
      file: this.getFilenameFor(id),
      storeEntries: true
    });

    await new Promise((resolve, reject) => {
      zip.on('ready', resolve);
      zip.on('error', reject);
    });

    function extract(input: string, output: string) {
      return new Promise((resolve, reject) => {
        zip.extract(input, in_directory + "/" + output, (err: any) => {
          if (err) reject(err);
          resolve();
        });
      });
    }

    // Extract top file, pdb file and itps
    await extract(infos.top.name, "molecule.top");
    await extract(infos.pdb.name, "molecule.pdb");

    for (const itp of infos.itp) {
      await extract(itp.name, itp.name);
    }

    // this throws ?
    zip.close();

    return {
      pdb: in_directory + '/molecule.pdb',
      top: in_directory + '/molecule.top',
      itps: infos.itp.map(e => `${in_directory}/${e.name}`),
    }
  }

  /**
   * Get the file infos attached to a save.
   */
  async getInfo(file_id: string): Promise<MoleculeSaveInfo | undefined> {
    if (await this.exists(file_id)) {
      return JSON.parse(await FsPromise.readFile(this.getInfoFilenameFor(file_id), "utf-8"));
    }
  }

  /**
   * Get the ZIP filename full path.
   *
   * You can use it in express in order to make the client download ZIP
   * ```ts
   * const filename = MoleculeOrganizer.getFilenameFor("139284920");
   * res.download(filename);
   * ```
   */
  getFilenameFor(file_id: string) {
    console.log( file_id)
    return MOLECULE_ROOT_DIR + file_id + ".zip";
  }

  protected getInfoFilenameFor(file_id: string) {
    return MOLECULE_ROOT_DIR + file_id + ".json";
  }

  /**
   * List saves available.
   * 
   * You can filter files you want by using a predicate on filename.
   */
  async list(predicate?: (file: string) => boolean): Promise<MoleculeSaveInfo[]> {
    const files = await FsPromise.readdir(MOLECULE_ROOT_DIR);

    let json_files = files.filter(f => f.endsWith('.json'));

    if (predicate) {
      json_files = json_files.filter(predicate);
    }

    // Todo in a worker ?
    return Promise.all(
      json_files.map(async f => JSON.parse(await FsPromise.readFile(MOLECULE_ROOT_DIR + f, "utf-8")))
    );
  }

  async existsMany(...file_ids: string[]) {
    for (let file_id of file_ids) {
      const a = await this.exists(file_id);
      if (!a) return a;
    }
    return true;
  }

  async exists(file_id: string) { // exists(...file_ids: string[])
    return FsPromise.access(this.getFilenameFor(file_id), fs.constants.F_OK).then(() => true).catch(() => false);
  }

  hash(file_id: string) {
    return md5File(this.getFilenameFor(file_id));
  }

  /**
   * Remove a save.
   */
  async remove(file_id: string) {
    try {
      await Promise.all([
        FsPromise.unlink(this.getFilenameFor(file_id)),
        FsPromise.unlink(this.getInfoFilenameFor(file_id))
      ]);
      logger.debug("Removed save ID #" + file_id);
    } catch { }
  }

  async removeAll() {
    const dir = await FsPromise.readdir(MOLECULE_ROOT_DIR);

    for (const file of dir) {
      await FsPromise.unlink(MOLECULE_ROOT_DIR + file);
    }
  }


  async createSymlinksInTmpDir(
    dir: string,
    pdb: Express.Multer.File | SimuFile,
    top: Express.Multer.File | SimuFile,
    itps: (Express.Multer.File | SimuFile)[],
    maps: (Express.Multer.File | SimuFile)[]
  ) {
    const pdb_name = path.basename(pdb.originalname);
    const full_pdb_name = dir + "/" + pdb_name;

    // Check if pdb has extension
    const pdb_ext = pdb_name.split('.')[1];
    if (!pdb_ext) {
      throw new Error("Uploaded PDB/GRO file must have an extension");
    }
    else if (pdb_ext !== 'pdb' && pdb_ext !== 'gro') {
      throw new Error("Uploaded PDB/GRO file must file extension '.pdb' or '.gro'.");
    }

    await FsPromise.symlink(pdb.path, dir + "/" + pdb_name);

    const top_name = generateSnowflake() + '.top';
    const full_top_name = dir + "/" + top_name;

    await FsPromise.symlink(top.path, dir + "/" + top_name);

    const full_itp_files: string[] = [];
    for (const file of itps) {
      const itp_name = basenameWithoutExt(file.originalname) + '.itp';
      full_itp_files.push(dir + "/" + itp_name);

      await FsPromise.symlink(file.path, dir + "/" + itp_name);
    }

    const full_map_files: string[] = [];
    for (const file of maps) {
      const map_name = basenameWithoutExt(file.originalname) + '.map';
      full_map_files.push(dir + "/" + map_name);

      await FsPromise.symlink(file.path, dir + "/" + map_name);
    }

    return {
      pdb: full_pdb_name,
      top: full_top_name,
      itps: full_itp_files,
      maps: full_map_files,
    };
  }

  /**
   * Create a JSZip object with the following data:
   * 
   * @param itps_path Path to ITP files
   * @param maps_path Path to MAP files
   * @param conect_pdb Path to CONECT-ed PDB
   * @param full_top Path to built TOP
   * @param top_name Original basename of the TOP
   * @param zip_destination_path ZIP path destination
   */
  async zipFromPaths(
    itps_path: string[],
    maps_path: string[],
    conect_pdb: string,
    full_top: string,
    top_name: string,
    zip_destination_path: string,
  ) {
    // Create ZIP
    const zip = new JSZip();

    const itp_files_info: FileSaveInfo[] = [];
    const map_files_info: FileSaveInfo[] = [];
    const targets = [itp_files_info, map_files_info];

    // Copy map and itps
    for (const item of [itps_path, maps_path]) {
      const target = targets.shift()!;

      for (const file of item) {
        const f_name = path.basename(file);
        const content = await FsPromise.readFile(file);
        zip.file(f_name, content);

        target.push({
          size: content.length,
          name: f_name,
        });
      }
    }

    const top_content = await FsPromise.readFile(full_top);
    zip.file(top_name, top_content);

    const pdb_content = await FsPromise.readFile(conect_pdb);
    zip.file(path.basename(conect_pdb), pdb_content);

    logger.debug("[MOLECULE-ORGANIZER] Saving in-memory ZIP file to disk.");

    await new Promise((resolve, reject) => {
      zip
        .generateNodeStream({ streamFiles: true, compression: "DEFLATE", compressionOptions: { level: 6 } })
        .pipe(fs.createWriteStream(zip_destination_path))
        .on('finish', () => {
          resolve();
        })
        .on('error', e => {
          reject(e);
        });
    });

    return {
      pdb_length: pdb_content.length,
      top_length: top_content.length,
      itp_files_info,
      map_files_info,
    }
  }

  /**
   * Check, compress and save a zipped version of the given ITP and PDB file.
   * 
   * Returns a save information.
   * 
   * TODO & verification of ITP+PDB
   * 
   * pdb_file can also be a GRO file, it doesn't matter. It is converted by GROMACS.
   */
  async save(
    itp_files: (Express.Multer.File | SimuFile)[],
    pdb_file: (Express.Multer.File | SimuFile),
    top_file: (Express.Multer.File | SimuFile),
    map_files: (Express.Multer.File | SimuFile)[],
    force_field: string, simple_force_field: boolean = true
  ): Promise<MoleculeSave> {
    logger.verbose("[MOLECULE-ORGANIZER] Saving upload " + itp_files.map(e => e.originalname).join(', ') + " and " + pdb_file.originalname);

    // TODO check ITP and PDB
    // ----------------------

    // Copy the files into a tmp dir
    const use_tmp_dir = await TmpDirHelper.get();

    logger.debug("[MOLECULE-ORGANIZER] Symlinking files into a temporary directory: " + use_tmp_dir + ".");

    // const {
    //   pdb: pdb_path,
    //   top: top_path,
    //   itps: itps_path,
    //   maps: maps_path
    // } = await this.createSymlinksInTmpDir(use_tmp_dir, pdb_file, top_file, itp_files, map_files);

    logger.debug("[MOLECULE-ORGANIZER] Creating extended TOP file for " + pdb_file.originalname + ". " + force_field);
    // Create the modified TOP and the modified pdb
    let full_top; 
    const itps_multer_path = itp_files.map(itp => itp.path)
    try {
      const ffForTop = simple_force_field ? "simple_" + force_field : force_field

      full_top = await Martinizer.createTopFileToString(top_file.path, itp_files.map(itp => itp.path), ffForTop);
      console.log("Return full top", full_top)
    } catch (e) {
      console.error(e)
      logger.warn("[MOLECULE-ORGANIZER] Unable to create extended TOP file. Maybe the ITPs are incorrects.");

      return Errors.throw(ErrorType.InvalidMoleculeFiles, {
        dir: use_tmp_dir,
        error: e,
      });
    }

    logger.debug("[MOLECULE-ORGANIZER] Extended TOP file created");

    let conectOutput: {pdb: string }; 
    logger.debug("[MOLECULE-ORGANIZER] Creating PDB with CONECT entries for " + pdb_file.originalname + ".");
    try {

      //var conectOutput = await Martinizer.createPdbWithConect(pdb_path, full_top, false, force_field, itps_path);

      const formatted_itp_paths: {[name: string]: string} = {}
      for (const itpPath of itps_multer_path){
        formatted_itp_paths[path.basename(itpPath)] = itpPath
      }

      conectOutput = await Martinizer.createPdbWithConectFromStream(pdb_file.path, "pdb", full_top, false, force_field, use_tmp_dir, formatted_itp_paths)

      console.log("CONECT OUTPUT", conectOutput)

    } catch (e) {
      logger.warn("[MOLECULE-ORGANIZER] Unable to create full PDB with GROMACS. Provided files might be incorrects.");

      return Errors.throw(ErrorType.InvalidMoleculeFiles, {
        dir: use_tmp_dir,
        error: e,
      });
    }

    // @ts-ignore
    logger.debug("[MOLECULE-ORGANIZER] CONECT-ed PDB created: " + path.basename(conectOutput.pdb) + ".");

    // Compressing and saving
    const save_id = generateSnowflake();
    const zip_name = this.getFilenameFor(save_id);
    const info_name = this.getInfoFilenameFor(save_id);

    const top_name = path.basename(top_file.path);
    const final_pdb = conectOutput.pdb
    const final_itps = await Promise.all(itp_files.map(async(itp) => {
      const symlink_path = use_tmp_dir + "/" + path.basename(itp.originalname)
      await FsPromise.symlink(itp.path, symlink_path)
      return symlink_path
      }))
    
    const final_maps = await Promise.all(map_files.map(async(map) => {
        const symlink_path = use_tmp_dir + "/" + path.basename(map.originalname)
        await FsPromise.symlink(map.path, symlink_path)
        return symlink_path
    }))

    console.log("final_itps", final_itps)
    
    const final_top = top_file.originalname.endsWith('.top') ?use_tmp_dir + "/" + path.basename(top_file.originalname) : use_tmp_dir + "/" + path.basename(top_file.originalname) + ".top"
    console.log("final_top", final_top)

    await FsPromise.symlink(top_file.path, final_top)    

    // Compress and get save data
    const {
      pdb_length,
      top_length,
      itp_files_info,
      map_files_info
    } = await this.zipFromPaths(
      final_itps,
      final_maps,
      conectOutput.pdb,
      final_top,
      top_name,
      zip_name
    );

    logger.debug("[MOLECULE-ORGANIZER] Computing MD5 hash of ZIP file.");

    // Calculate hash
    const hash = await md5File(zip_name);

    // TODO write better infos of the ZIP in the JSON
    // Send the save data
    const infos: MoleculeSaveInfo = {
      pdb: {
        size: pdb_length,
        // @ts-ignore
        name: path.basename(conectOutput.pdb)
      },
      top: {
        size: top_length,
        name: top_name,
      },
      itp: itp_files_info,
      map: map_files_info,
      hash,
      force_field
    };

    await FsPromise.writeFile(info_name, JSON.stringify(infos));

    logger.debug("[MOLECULE-ORGANIZER] ZIP has been created and saved with save ID #" + save_id + ".");

    return {
      id: save_id,
      name: zip_name,
      infos,
    };
  }
};

export default MoleculeOrganizer;

export interface FileSaveInfo {
  size: number;
  name: string;
}

export interface MoleculeSaveInfo {
  pdb: FileSaveInfo;
  itp: FileSaveInfo[];
  top: FileSaveInfo;
  map: FileSaveInfo[];
  hash: string;
  force_field: string;
}

export interface MoleculeSave {
  id: string;
  name: string;
  infos: MoleculeSaveInfo;
}
