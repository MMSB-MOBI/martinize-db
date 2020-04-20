import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles } from '../../helpers';
import Uploader from '../Uploader';
import { promises as FsPromise } from 'fs';
import MembraneBuilder, { LipidMap, PbcString, AvailablePbcStrings, InsaneSettings } from '../../Builders/MembraneBuilder';
import { SettingsJson } from '../../types';
import { SETTINGS_FILE } from '../../constants';
import Errors, { ErrorType } from '../../Errors';
import TmpDirHelper from '../../TmpDirHelper/TmpDirHelper';
import path from 'path';

const MembraneBuilderRouter = Router();

// Middleware that wipe uploaded files after request
MembraneBuilderRouter.use((req, res, next) => {
  function after() {
    // Response is sended
    cleanMulterFiles(req);
    res.removeListener('finish', after);
  }

  res.once('finish', after);
  next();
});

function detectType(ext: string) {
  switch (ext) {
    case 'itp': return 'chemical/x-include-topology';
    case 'top': return 'chemical/x-topology';
    case 'pdb': return 'chemical/x-pdb';
  }
  return '';
}

async function getFormattedFile(file: string) {
  const name = path.basename(file);
  const type = detectType(file.slice(file.indexOf('.') + 1));

  return {
    name,
    type,
    content: await FsPromise.readFile(file, 'utf-8'),
  };
}

function checkPbc(pbc: string) : pbc is PbcString {
  return AvailablePbcStrings.includes(pbc as any);
}

/**
 * Create a new molecule.
 * 
 * Expect a enctype `multipart/form-data`, because you should have two files in it.
 * 
 * Expect file `itp=itp_file`, `pdb=pdb_file`, `gro=gro_file`.
 * You must specify an ITP, and at least one PDB or one GRO file.
 * 
 * 
 */
MembraneBuilderRouter.post('/', Uploader.fields([
  { name: 'itp', maxCount: 99 }, 
  { name: 'top', maxCount: 1 },
  { name: 'pdb', maxCount: 1 },
]), (req, res) => {
  (async () => {
    // Init
    const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));
    const files = req.files as { [fieldname: string]: Express.Multer.File[] };
    const molecule_id = req.body.from_id as string | undefined;
    const { pbc, box, force_field, lipids: lipids_str, upper_leaflet: upper_leaflet_str } = req.body;

    const opts: Partial<InsaneSettings> = {};

    // Test force field
    
    if (!settings.force_fields.includes(force_field)) {
      return Errors.throw(ErrorType.InvalidForceField);
    }
    
    const molecule_entries = {
      molecule_pdb: "",
      molecule_top: "",
      molecule_itps: [] as string[],
    };

    if (molecule_id) {
      // from molecule id
      const { pdb, itps, top } = await MembraneBuilder.prepareRunWithDatabaseMolecule(molecule_id);
      molecule_entries.molecule_itps = itps;
      molecule_entries.molecule_pdb = pdb;
      molecule_entries.molecule_top = top;
    }
    else {
      // except them from files
      if (!files || !files.itp || !files.top || !files.pdb) {
        return Errors.throw(ErrorType.MissingFiles);
      }
      if (!files.itp.length || !files.top.length || !files.pdb.length) {
        return Errors.throw(ErrorType.MissingFiles);
      }

      // todo test file size !
      const tmp_dir = await TmpDirHelper.get();

      // Creating a symlink for pdb/top
      await FsPromise.symlink(files.top[0].path, tmp_dir + "/full.top");
      molecule_entries.molecule_top = tmp_dir + "/full.top";
      await FsPromise.symlink(files.pdb[0].path, tmp_dir + "/output.pdb");
      molecule_entries.molecule_pdb = tmp_dir + "/output.pdb";
      
      // Symlink for itps
      for (const itp of files.itp) {
        let name = itp.originalname;
        if (!name.endsWith('.itp')) {
          name += ".itp";
        }
        await FsPromise.symlink(itp.path, tmp_dir + "/" + name);
        molecule_entries.molecule_itps.push(tmp_dir + "/" + name);
      }

      // ready !
    }

    if (!lipids_str) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    // Parse lipid str
    const lipids = (lipids_str as string).split(',').map(e => { 
      const res = e.split(':'); 
      
      if (res.length > 1) {
        return [res[0], parseInt(res[1], 10)] as [string, number];
      }
      return [res[0], 1] as [string, number];
    });

    // If upper leaflet, parse it
    let upper_leaflet: LipidMap = [];
    if (upper_leaflet_str) {
      upper_leaflet = (upper_leaflet_str as string).split(',').map(e => { 
        const res = e.split(':'); 
        
        if (res.length > 1) {
          return [res[0], parseInt(res[1], 10)] as [string, number];
        }
        return [res[0], 1] as [string, number];
      });
    }

    // Parse settings
    if (pbc && checkPbc(pbc)) {
      opts.pbc = pbc;
    }
    if (box) {
      const items = (box as string).split(',').map(e => parseInt(e, 10));

      if (!items.every(e => !isNaN(e) && e > 0)) {
        return Errors.throw(ErrorType.Format);
      }
      opts.box = items;
    }
    // todo other settings

    const { pdbs: { water, no_water }, top, itps } = await MembraneBuilder.run({
      force_field,
      lipids,
      upper_leaflet,
      ...molecule_entries,
      settings: opts,
    });

    res.json({
      water: await getFormattedFile(water),
      no_water: await getFormattedFile(no_water),
      top: await getFormattedFile(top),
      itps: await Promise.all(itps.map(i => getFormattedFile(i))),
    });
  })().catch(errorCatcher(res));
});

MembraneBuilderRouter.all('/', methodNotAllowed(['POST']));

export default MembraneBuilderRouter;
