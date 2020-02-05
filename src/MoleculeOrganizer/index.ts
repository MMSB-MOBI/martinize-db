import fs, { promises as FsPromise } from 'fs';
import { MOLECULE_ROOT_DIR, UPLOAD_ROOT_DIR } from '../constants';
import logger from '../logger';
import { getNameAndPathOfUploadedFile, generateSnowflake } from '../helpers';
import JSZip from 'jszip';

export default new class MoleculeOrganizer {
  constructor() {
    try {
      fs.mkdirSync(MOLECULE_ROOT_DIR);
    } catch {}
  }

  /**
   * Get a ZIP, read by JSZip
   * @param file_id 
   */
  async get(file_id: string) : Promise<[JSZip, MoleculeSaveInfo] | undefined> {
    if (await this.exists(file_id)) {
      const infos: MoleculeSaveInfo = JSON.parse(await FsPromise.readFile(MOLECULE_ROOT_DIR + file_id + ".json", "utf-8"));

      const zip_buffer = await FsPromise.readFile(MOLECULE_ROOT_DIR + file_id + ".zip");
      const zip = await JSZip.loadAsync(zip_buffer);
      
      return [
        zip,
        infos
      ];
    }
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
    return FsPromise.access(file_id + ".zip", fs.constants.F_OK).then(() => true).catch(() => false);
  }
  
  /**
   * Remove a save.
   */
  async remove(file_id: string) {
    try {
      await Promise.all([
        FsPromise.unlink(MOLECULE_ROOT_DIR + file_id + ".zip"),
        FsPromise.unlink(MOLECULE_ROOT_DIR + file_id + ".json")
      ]);
    } catch {}
  }

  /**
   * Check, compress and save a zipped version of the given ITP and PDB file.
   * 
   * Returns a save information.
   */
  async save(itp_file: string, pdb_file: string) : Promise<MoleculeSave> {
    logger.debug("Saving ", itp_file, " and ", pdb_file);

    const [itp_fn, itp_fp] = getNameAndPathOfUploadedFile(itp_file);
    const [pdb_fn, pdb_fp] = getNameAndPathOfUploadedFile(pdb_file);

    // TODO check ITP and PDB
    // ----------------------

    // Compressing and saving
    const save_id = generateSnowflake();
    const root_name = MOLECULE_ROOT_DIR + save_id;
    const zip_name = root_name + ".zip";
    const info_name = root_name + ".json";
    const pdb_name = pdb_fn.split('.')[0] + '.pdb';
    const itp_name = itp_fn.split('.')[0] + '.itp';

    const zip = new JSZip();
    const itp_content = await FsPromise.readFile(itp_fp);
    zip.file(itp_name, itp_content);

    const pdb_content = await FsPromise.readFile(pdb_fp);
    zip.file(pdb_name, pdb_content);

    await new Promise((resolve, reject) => {
      zip
        .generateNodeStream({ streamFiles: true })
        .pipe(fs.createWriteStream(zip_name))
        .on('finish', () => {
          resolve();
        })
        .on('error', e => {
          reject(e);
        });
    });

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
    };

    await FsPromise.writeFile(info_name, JSON.stringify(infos));

    return {
      id: save_id,
      name: zip_name,
      infos,
    };
  }
};

export interface FileSaveInfo {
  size: number;
  name: string;
}

export interface MoleculeSaveInfo {
  pdb?: FileSaveInfo;
  itp: FileSaveInfo;
  gro?: FileSaveInfo;
  type: "pdb" | "gro";
}

export interface MoleculeSave {
  id: string;
  name: string;
  infos: MoleculeSaveInfo;
}
