import fs, { promises as FsPromise } from 'fs';
import { MOLECULE_ROOT_DIR } from '../constants';
import logger from '../logger';
import { getNameAndPathOfUploadedFile, generateSnowflake } from '../helpers';
import JSZip from 'jszip';
import md5File from 'md5-file/promise';

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
      const infos: MoleculeSaveInfo = JSON.parse(await FsPromise.readFile(this.getInfoFilenameFor(file_id), "utf-8"));

      const zip_buffer = await FsPromise.readFile(this.getFilenameFor(file_id));
      const zip = await JSZip.loadAsync(zip_buffer);
      
      return [
        zip,
        infos
      ];
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
      json_files = files.filter(predicate);
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
  async save(itp_files: Express.Multer.File[], pdb_file: Express.Multer.File, force_field: string) : Promise<MoleculeSave> {
    logger.debug("Saving ", itp_files, " and ", pdb_file);

    // TODO check ITP and PDB
    // ----------------------

    // Compressing and saving
    const save_id = generateSnowflake();
    const zip_name = this.getFilenameFor(save_id);
    const info_name = this.getInfoFilenameFor(save_id);

    const pdb_name = pdb_file.originalname.split('.')[0] + '.pdb';

    // Create ZIP
    const zip = new JSZip();

    let itp_files_info: FileSaveInfo[] = [];
    for (const file of itp_files) {
      const itp_name = file.originalname.split('.')[0] + '.itp';
      const itp_content = await FsPromise.readFile(file.path);
      zip.file(itp_name, itp_content);

      itp_files_info.push({
        size: itp_content.length,
        name: itp_name,
      });
    }

    const pdb_content = await FsPromise.readFile(pdb_file.path);
    zip.file(pdb_name, pdb_content);

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
        name: pdb_name
      },
      itp: itp_files_info,
      type: "pdb",
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
  pdb?: FileSaveInfo;
  itp: FileSaveInfo[];
  gro?: FileSaveInfo;
  type: "pdb" | "gro";
  hash: string;
  force_field: string;
}

export interface MoleculeSave {
  id: string;
  name: string;
  infos: MoleculeSaveInfo;
}
