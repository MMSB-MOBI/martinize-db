import { exec } from 'child_process';
import fs, { promises as FsPromise } from 'fs';
import os from 'os';
import Errors, { ErrorType } from '../Errors';


export interface MartinizeSettings {
  /** PDB file path */
  input: string;
  
  /** Ignore residues */
  ignore?: string[];

  /** Ignore hydrogens */
  ignh?: boolean;

  /** Force field */
  ff: string;

  /** Position restrains */
  position: 'none' | 'all' | 'backbone';

  /** Position restrain force const */
  posref_fc?: number;

  /** Use collagen */
  collagen?: boolean;

  /** Use dihedral */
  dihedral?: boolean;

  /** Elastic bounds */
  elastic?: boolean;
  /** Elastic force const */
  ef?: number;
  /** Elastic lower bound */
  el?: number;
  /** Elastic upper bound */
  eu?: number;
  /** Elastic decay alpha */
  ea?: number;
  /** Elastic decay power */
  ep?: number;
  /** Elastic remover minimum force */
  em?: number;
  /** List of bead names for elastic bound (comma separated in martinize) */
  eb?: string[];

  /** Use govs */
  use_go_virtual_sites?: boolean;

  /** Set neutral termini */
  neutral_termini?: boolean;
  /** Apply side chains corrections */
  side_chain_fix?: boolean;
  /** Cystein bounds */
  cystein_bridge?: string;
}

export const Martinizer = new class Martinizer {
  async run(settings: Partial<MartinizeSettings>) {
    const full: MartinizeSettings = Object.assign({}, {
      input: '',
      ff: 'martini22',
      position: 'none'
    }, settings);

    if (!full.input.trim()) {
      throw new Error("Invalid input file");
    }

    const basename = full.input;
    const with_ext = basename + '.pdb';

    // Check dssp ps
    let command_line = "martinize2 -f " + with_ext + " -x output.pdb -dssp mkdssp -ff " + full.ff + " -p " + full.position + " ";

    if (full.ignore) {
      command_line += " " + full.ignore.join(',');
    }
    if (full.ignh) {
      command_line += " -ignh";
    }

    if (full.posref_fc) {
      command_line += " -pf " + full.posref_fc.toString();
    }
    if (full.collagen) {
      command_line += " -collagen ";
    }
    if (full.dihedral) {
      command_line += " -ed ";
    }
    if (full.elastic) {
      command_line += " -elastic ";
    }
    if (full.ef) {
      command_line += " -ef " + full.ef.toString();
    }
    if (full.el) {
      command_line += " -el " + full.el.toString();
    }
    if (full.eu) {
      command_line += " -eu " + full.eu.toString();
    }
    if (full.ea) {
      command_line += " -ea " + full.ea.toString();
    }
    if (full.ep) {
      command_line += " -ep " + full.ep.toString();
    }
    if (full.em) {
      command_line += " -em " + full.em.toString();
    }
    if (full.eb) {
      command_line += " -eb " + full.eb.toString();
    }
    if (full.use_go_virtual_sites) {
      command_line += " -govs-include ";
    }
    if (full.neutral_termini) {
      command_line += " -nt ";
    }
    if (full.side_chain_fix) {
      command_line += " -scfix ";
    }
    if (full.cystein_bridge) {
      command_line += " -cys " + full.cystein_bridge;
    }

    await FsPromise.rename(basename, with_ext);

    try {
      // Run the command line
      const tmp_dir = os.tmpdir();
      const dir = await FsPromise.mkdtemp(tmp_dir + "/");
  
      const exists = await FsPromise.access(dir, fs.constants.F_OK).then(() => true).catch(() => false);
      if (!exists) {
        await FsPromise.mkdir(dir, { recursive: true });
      }
  
      await new Promise((resolve, reject) => {
        exec(command_line, { cwd: dir }, (err, stdout, stderr) => {
          if (err) {
            reject(err);
            return;
          }
          resolve([stdout, stderr]);
        });
      }) as [string, string];
  
      // Scan for files in result dir
      let pdb_file: string = dir + "/output.pdb";
      const itp_files: string[] = [];
  
      const exists_pdb = await FsPromise.access(pdb_file, fs.constants.F_OK).then(() => true).catch(() => false);
      if (!exists_pdb) {
        return Errors.throw(ErrorType.Server, { error: "Unable to find created pdb." });
      }

      // If go mode, we should compute map + run a python script to refresh ITPs files.
      if (settings.use_go_virtual_sites) {
        // Must create the go sites
        const [out, ] = await new Promise((resolve, reject) => {
          exec('cut -f1 -d \' \' output.pdb | grep -c ATOM', { cwd: dir }, (err, stdout, stderr) => {
            if (err) {
              reject(err);
              return;
            }
            resolve([stdout, stderr]);
          });
        }) as [string, string];
        
        // GET THE MAP FILE FROM A CUSTOM WAY.
        // For now, this throw an error.
        throw new Error("No way to get the MAP file.");

        // Ensure to have the script at the correct place...
        await new Promise((resolve, reject) => {
          exec(`python create_goVirt.py -s output.pdb -f MOLECULE.map --Natoms ${out.trim()} --moltype molecule_0`, { cwd: dir }, (err, stdout, stderr) => {
            if (err) {
              reject(err);
              return;
            }
            resolve([stdout, stderr]);
          });
        }) as [string, string];

        // Ok, ITP file refreshed.
      }

      for (const file of await FsPromise.readdir(dir)) {
        if (file.endsWith('.itp')) {
          itp_files.push(dir + "/" + file);
        }
      }

      return {
        pdb: pdb_file,
        itps: itp_files
      };
    } finally {
      await FsPromise.rename(with_ext, basename);
    }
  }
}();
