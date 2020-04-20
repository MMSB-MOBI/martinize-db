import fs, { promises as FsPromise } from 'fs';
import { MOLECULE_ROOT_DIR } from '../constants';
import logger from '../logger';
import { generateSnowflake } from '../helpers';
import JSZip from 'jszip';
import md5File from 'md5-file/promise';
import os from 'os';
import { Martinizer } from '../Builders/Martinizer';
import path from 'path';
import TmpDirHelper from '../TmpDirHelper/TmpDirHelper';
// @ts-ignore 
import NodeStreamZip from 'node-stream-zip';

export const MoleculeOrganizer = new class MoleculeOrganizer {
  constructor() {
    try {
      fs.mkdirSync(MOLECULE_ROOT_DIR);
    } catch {}
  }

  /**
   * Get a ZIP, read by JSZip.
   * 
   * If the save doesn't exists, returns `undefined`.
   */
  async get(file_id: string) : Promise<[JSZip, MoleculeSaveInfo] | undefined> {
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
  async getInfo(file_id: string) : Promise<MoleculeSaveInfo | undefined> {
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
  async list(predicate?: (file: string) => boolean) : Promise<MoleculeSaveInfo[]> {
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

  async exists(file_id: string) {
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
    } catch {}
  }

  async removeAll() {
    const dir = await FsPromise.readdir(MOLECULE_ROOT_DIR);

    for (const file of dir) {
      await FsPromise.unlink(MOLECULE_ROOT_DIR + file);
    }
  }

  /**
   * Check, compress and save a zipped version of the given ITP and PDB file.
   * 
   * Returns a save information.
   * 
   * TODO make support of GRO files & verification of ITP+PDB
   */
  async save(itp_files: Express.Multer.File[], pdb_file: Express.Multer.File, top_file: Express.Multer.File, force_field: string) : Promise<MoleculeSave> {
    logger.verbose("Saving upload " + itp_files.map(e => e.originalname).join(', ') + " and " + pdb_file.originalname);

    // TODO check ITP and PDB
    // ----------------------

    logger.debug("Copying files into a temporary directory");
    // Copy the files into a tmp dir
    const use_tmp_dir = await TmpDirHelper.get();

    const pdb_name = path.basename(pdb_file.originalname);
    const full_pdb_name = use_tmp_dir + "/" + pdb_name;

    // Check if pdb has extension
    const pdb_ext = pdb_name.split('.')[1];
    if (!pdb_ext) {
      throw new Error("Uploaded PDB/GRO file must have an extension");
    }
    else if (pdb_ext !== 'pdb' && pdb_ext !== 'gro') {
      throw new Error("Uploaded PDB/GRO file must file extension '.pdb' or '.gro'.");
    }

    await FsPromise.copyFile(pdb_file.path, use_tmp_dir + "/" + pdb_name);

    const top_name = top_file.originalname.split('.')[0] + '.top';
    const full_top_name = use_tmp_dir + "/" + top_name;

    await FsPromise.copyFile(top_file.path, use_tmp_dir + "/" + top_name);
   
    const full_itp_files: string[] = [];
    for (const file of itp_files) {
      const itp_name = file.originalname.split('.')[0] + '.itp';
      full_itp_files.push(use_tmp_dir + "/" + itp_name);

      await FsPromise.copyFile(file.path, use_tmp_dir + "/" + itp_name);
    }

    logger.debug("Creating extended TOP file for " + pdb_file.originalname);
    // Create the modified TOP and the modified pdb
    const { top: full_top } = await Martinizer.createTopFile(use_tmp_dir, full_top_name, full_itp_files, force_field);
    logger.debug("Extended TOP file created: " + path.basename(full_top));

    logger.debug("Creating PDB with CONECT entries for " + pdb_file.originalname);
    const full_pdb = await Martinizer.createPdbWithConect(full_pdb_name, full_top, use_tmp_dir);
    logger.debug("CONECT-ed PDB created: " + path.basename(full_pdb));

    // Compressing and saving
    const save_id = generateSnowflake();
    const zip_name = this.getFilenameFor(save_id);
    const info_name = this.getInfoFilenameFor(save_id);

    // Create ZIP
    const zip = new JSZip();

    let itp_files_info: FileSaveInfo[] = [];
    for (const file of full_itp_files) {
      const itp_name = path.basename(file);
      const itp_content = await FsPromise.readFile(file);
      zip.file(itp_name, itp_content);

      itp_files_info.push({
        size: itp_content.length,
        name: itp_name,
      });
    }

    const top_content = await FsPromise.readFile(full_top);
    zip.file(path.basename(top_name), top_content);

    const pdb_content = await FsPromise.readFile(full_pdb);
    zip.file(path.basename(full_pdb), pdb_content);

    await new Promise((resolve, reject) => {
      zip
        .generateNodeStream({ streamFiles: true, compression: "DEFLATE", compressionOptions: { level: 6 } })
        .pipe(fs.createWriteStream(zip_name))
        .on('finish', () => {
          resolve();
        })
        .on('error', e => {
          reject(e);
        });
    });

    // Calculate hash
    const hash = await md5File(zip_name);


    // TODO write better infos of the ZIP in the JSON
    const infos: MoleculeSaveInfo = {
      pdb: {
        size: pdb_content.length,
        name: path.basename(full_pdb)
      },
      top: {
        size: top_content.length,
        name: path.basename(top_name),
      },
      itp: itp_files_info,
      hash,
      force_field
    };

    await FsPromise.writeFile(info_name, JSON.stringify(infos));

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
  hash: string;
  force_field: string;
}

export interface MoleculeSave {
  id: string;
  name: string;
  infos: MoleculeSaveInfo;
}
