import { VanDerWaalsRadius } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import nano = require("nano");
import { promises as FsPromise } from 'fs';
import { exec } from 'child_process';

export default class RadiusDatabase extends AbstractDatabase<VanDerWaalsRadius> {
  protected static readonly FORCE_FIELD_TO_FILE_NAME: { [ff: string]: string | string[] } = {
    martini304: 'martini_v3.0.4.itp',
    elnedyn22p: 'martini_v2.2P.itp',
    elnedyn22: 'martini_v2.2.itp',
    elnedyn: 'martini_v2.2.itp',
    martini22: 'martini_v2.2.itp',
    martini22p: 'martini_v2.2P.itp',
    martini23: 'martini_v2.3P.itp',
  };

  /**
   * Create a new set a custom atoms for {force_field}.
   */
  async create(force_field: string) {
    const atoms: { [name: string]: number } = {};
    
    if (!(force_field in RadiusDatabase.FORCE_FIELD_TO_FILE_NAME)) {
      throw new Error("Force field " + force_field + " is not referenced.");
    }

    const itp_files = RadiusDatabase.FORCE_FIELD_TO_FILE_NAME[force_field];

    for (const file of (Array.isArray(itp_files) ? itp_files : [itp_files])) {
      const [out, ] = await new Promise((resolve, reject) => {
        exec(`grep -P "^\\s*(?'g1'\\S+)\\s+(\\k{g1})\\b" "${file}" | sort | cut -d' ' -f-2,4- | uniq`, (err, stdout, stderr) => {
          if (err) {
            reject(err);
            return;
          }
          resolve([stdout, stderr]);
        });
      }) as [string, string];

      for (const line of out.split('\n')) {
        const l = line.trim();
        if (!l) {
          continue;
        }

        const [residue, , r1, r2] = l.split(/\s+/);
        atoms[residue] = Number(r1) + Number(r2);
      }
    }

    // Atoms is okay

    const exists = await this.exists(force_field);
    if (exists) {
      await this.delete(await this.get(force_field));
    }

    // Insert the thing
    const radius: VanDerWaalsRadius = {
      id: force_field,
      atoms
    };

    const ok = await this.save(radius);

    if (!ok.ok) {
      throw new Error("Unable to insert force field.");
    }

    return radius;
  }

  protected async createOrCheckExistance(force_field: string) {
    const exists = await this.exists(force_field);
    if (!exists) {
      return this.create(force_field);
    }
    return this.get(force_field);
  }

  /**
   * For a given molecule, get the related radius for every atom linked to its ITPs.
   * 
   * @param itp_files ITP **content** (not their filename!).
   */
  async getRadius(force_field: string, itp_files: string[]) {
    const atoms = (await this.createOrCheckExistance(force_field)).atoms;
    const subset: typeof atoms = {};

    for (const file of itp_files) {
      let seen_atoms = false;

      for (const line of file.split('\n')) {
        if (!seen_atoms && line !== '[ atoms ]') {
          continue;
        }
        seen_atoms = true;

        if (!line.trim()) {
          break;
          // atoms are read.
        }

        // Typical line is: 1 Q5  1 GLY BB   1   1
        const residue = line.split(/\s+/)[1];

        if (residue in atoms) {
          subset[residue] = atoms[residue];
        }
      }
    }

    return subset;
  }
}
