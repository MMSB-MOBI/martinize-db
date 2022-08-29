import { VanDerWaalsRadius } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import fs, { promises as FsPromise } from 'fs';
import readline from 'readline';
import { FORCE_FIELD_DIR } from "../constants";

export default class RadiusDatabase extends AbstractDatabase<VanDerWaalsRadius> {
  static readonly FORCE_FIELD_TO_FILE_NAME: { [ff: string]: string[] } = {
    martini3001 : ['martini_v3.0.0.itp', 'martini_v3.0.0_ions_v1.itp', 'martini_v3.0.0_solvents_v1.itp'],
    elnedyn22p: ['martini_v2.2P.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    elnedyn22: ['martini_v2.2.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    elnedyn: ['martini_v2.2.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    martini22: ['martini_v2.2.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    simple_martini22: ['martini_v2.2.itp'], 
    martini22p: ['martini_v2.2P.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    simple_martini22p: ['martini_v2.2P.itp'],
    martini23p: ['martini_v2.3P.itp', 'martini_v2.0_ions.itp', 'martini_v2.0_solvents.itp'],
    martini23_CNP: ['martini_v2.3_CNP.itp', 'martini_v2.0_ions.itp'],
    simple_martini23_CNP : ['martini_v2.3_CNP.itp']
  };

  static readonly FORCE_FIELD_TO_MARTINI_VERSION: { [ff: string]: string } = {
    martini3001: '3_0',
    elnedyn22p: '2_2',
    elnedyn22: '2_2',
    elnedyn: '2_2',
    martini22: '2_2',
    simple_martini22 : '2_2', 
    martini22p: '2_2',
    simple_martini22p : '2_2',
    martini23p: '2_2',
    martini23_CNP : '2_2'
  };

  static getFilesForForceField(name: string) {
    if (!(name in this.FORCE_FIELD_TO_FILE_NAME)) {
      return [];
    }

    const files = this.FORCE_FIELD_TO_FILE_NAME[name];
    if (typeof files === 'string') {
      return [files];
    }
    return files;
  }

  static getCompleteFilesForForceField(name : string) {
    if (!(name in this.FORCE_FIELD_TO_FILE_NAME)) {
      return [];
    }

    const files = this.FORCE_FIELD_TO_FILE_NAME[name];

    return files.map(f => FORCE_FIELD_DIR + "/" + f) 
  }

  /*static getFilesForForceFieldObj(name : string) : {[file_name : string] : string} {
    const filesObj : {[file_name : string] : string}  = {}
    if (!(name in this.FORCE_FIELD_TO_FILE_NAME)) {
      return {};
    }

    const files = this.FORCE_FIELD_TO_FILE_NAME[name];
    for (const f of files){
      filesObj[f] = FORCE_FIELD_DIR + "/" + f
    }
    return filesObj;
  }*/

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

  protected async getAtomMapAndRadius(force_field: string, itp_files: (string | NodeJS.ReadableStream)[]) {
    const atoms = (await this.createOrCheckExistance(force_field)).atoms;
    const subset: typeof atoms = {};
    const index_to_residue: { [index: number]: string } = {};
    let i = 1;

    // ITP files must be in TOP file include order !!
    for (const file of itp_files) {
      let seen_atoms = false;

      const rl = readline.createInterface({
        input: typeof file === 'string' ? fs.createReadStream(file) : file,
        crlfDelay: Infinity, // crlfDelay option to recognize all instances of CR LF in file as a single line break.
      });

      for await (const line of rl) {
        if (!seen_atoms && line !== '[ atoms ]') {
          continue;
        }
        seen_atoms = true;
        if (line === '[ atoms ]') {
          continue;
        }

        if (!line.trim()) {
          break;
          // atoms are read.
        }

        // Typical line is: 1 Q5  1 GLY BB   1   1
        const splitted_line = line.split(/\s+/).filter(e => e);
        let residue = splitted_line[1];

        if (residue.startsWith('molecule_0_')) {
          // Go virt site
          // Original name 'CA'
          residue = 'CA';
        }

        index_to_residue[i] = residue;
        i++;

        if (residue in atoms) {
          subset[residue] = atoms[residue];
        }
      }

      rl.close();
    }

    return [subset, index_to_residue] as [typeof atoms, { [index: number]: string }];
  }

  async transformPdb(pdb: string | NodeJS.ReadableStream, itp_files: (string | NodeJS.ReadableStream)[], force_field: string) {
    // Read the pdb and transform it
    const lines: string[] = [];
    
    const [radius, indexes] = await this.getAtomMapAndRadius(force_field, itp_files);

    const rl = readline.createInterface({
      input: typeof pdb === 'string' ? fs.createReadStream(pdb) : pdb,
      crlfDelay: Infinity, // crlfDelay option to recognize all instances of CR LF in file as a single line break.
    });

    for await (const line of rl) {
      if (!line.startsWith('ATOM')) {
        lines.push(line);
        continue;
      }

      const atom_number = Number(line.slice(7, 12));
      
      // Atom name needs to be changed.
      if (atom_number in indexes) {
        lines.push(
          line.slice(0, 12) +
          indexes[atom_number].padEnd(4, ' ') +
          line.slice(16)
        );
      }
      else {
        lines.push(line);
      }
    }

    return [
      lines.join('\n'),
      radius
    ] as [string, typeof radius];
  }

  /**
   * For a given molecule, get the related radius for every atom linked to its ITPs.
   * 
   * @param itp_files ITP **filename**.
   */
  async getRadius(force_field: string, itp_files: (string | NodeJS.ReadableStream)[]) {
    return (await this.getAtomMapAndRadius(force_field, itp_files))[0];
  }
}
