import { VanDerWaalsRadius } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import nano = require("nano");
import fs, { promises as FsPromise } from 'fs';
import readline from 'readline';
import { FORCE_FIELD_DIR } from "../constants";

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
    const files = Array.isArray(itp_files) ? itp_files : [itp_files];

    for (const file of files) {
      const rl = readline.createInterface({
        input: fs.createReadStream(FORCE_FIELD_DIR + file),
        crlfDelay: Infinity, // crlfDelay option to recognize all instances of CR LF in file as a single line break.
      });

      for await (const line of rl) {
        const l = line.trim();
        
        if (!l || !l.match(/^(?<g1>\S+)\s+(\k<g1>)\b/)) {
          continue;
        }

        const [residue, , r1, r2] = l.split(/\s+/);
        if (residue && r1 && r2) {
          atoms[residue] = Number(r1) + Number(r2);
        }
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
