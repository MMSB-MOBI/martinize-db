import fs, { promises as FsPromise } from 'fs';
import { MOLECULE_ROOT_DIR, UPLOAD_ROOT_DIR } from '../constants';
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
  
  /**
   * Remove a save.
   */
  async remove(file_id: string) {
    try {
      await Promise.all([
        FsPromise.unlink(this.getFilenameFor(file_id)),
        FsPromise.unlink(this.getInfoFilenameFor(file_id))
      ]);
    } catch {}
  }

  /**
   * Check, compress and save a zipped version of the given ITP and PDB file.
   * 
   * Returns a save information.
   * 
   * TODO make support of GRO files
   */
  async save(itp_file: string, pdb_file: string) : Promise<MoleculeSave> {
    logger.debug("Saving ", itp_file, " and ", pdb_file);

    const [itp_fn, itp_fp] = getNameAndPathOfUploadedFile(itp_file);
    const [pdb_fn, pdb_fp] = getNameAndPathOfUploadedFile(pdb_file);

    // TODO check ITP and PDB
    // ----------------------

    // Compressing and saving
    const save_id = generateSnowflake();
    const zip_name = this.getFilenameFor(save_id);
    const info_name = this.getInfoFilenameFor(save_id);
    const pdb_name = pdb_fn.split('.')[0] + '.pdb';
    const itp_name = itp_fn.split('.')[0] + '.itp';

    const zip = new JSZip();
    const itp_content = await FsPromise.readFile(itp_fp);
    zip.file(itp_name, itp_content);

    const pdb_content = await FsPromise.readFile(pdb_fp);
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
      itp: {
        size: itp_content.length,
        name: itp_name
      },
      type: "pdb",
      hash
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
  itp: FileSaveInfo;
  gro?: FileSaveInfo;
  type: "pdb" | "gro";
  hash: string;
}

export interface MoleculeSave {
  id: string;
  name: string;
  infos: MoleculeSaveInfo;
}
