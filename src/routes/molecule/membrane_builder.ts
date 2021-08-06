import { Router } from 'express';
import { methodNotAllowed, errorCatcher, cleanMulterFiles } from '../../helpers';
import Uploader from '../Uploader';
import { promises as FsPromise } from 'fs';
import MembraneBuilder, { LipidMap, PbcString, AvailablePbcStrings, InsaneSettings, InsaneError } from '../../Builders/MembraneBuilder';
import { SettingsJson } from '../../types';
import { SETTINGS_FILE } from '../../constants';
import Errors, { ErrorType } from '../../Errors';
import TmpDirHelper from '../../TmpDirHelper';
import path from 'path';
import { Martinizer } from '../../Builders/Martinizer';
import logger from '../../logger';
import { inspect } from 'util';

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

function isBoolTrue(v: string) {
  return v === "true" || v === "t" || v === "1";
}

// Items automatically coerced to numbers when presents
const VALID_BODY_ITEMS = [
  'area_per_lipid', 
  'area_per_lipid_upper',
  'random_kick_size',
  'bead_distance',
  'grid_spacing',
  'hydrophobic_ratio',
  'fudge',
  'shift_protein',
] as const;


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
    function convertOrThrow(str: string) : number {
      const n = Number(str);
      if (isNaN(n)) {
        return Errors.throw(ErrorType.Format);
      }

      return n;
    }

    //console.log(req.body)

    // Init
    const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));
    const files = req.files as { [fieldname: string]: Express.Multer.File[] };
    const molecule_id = req.body.from_id as string | undefined;
    const { 
      pbc, 
      box, 
      lipids: lipids_str, 
      upper_leaflet: upper_leaflet_str 
    } = req.body;
    let force_field: string = req.body.force_field;

    const opts: Partial<InsaneSettings> = {};

    // Test force field
    
    const molecule_entries = {
      molecule_pdb: "",
      molecule_top: "",
      molecule_itps: [] as string[],
    };

    if (molecule_id) {
      // from molecule id
      const { pdb, itps, top, force_field: ff } = await MembraneBuilder.prepareRunWithDatabaseMolecule(molecule_id);
      molecule_entries.molecule_itps = itps;
      molecule_entries.molecule_pdb = pdb;
      molecule_entries.molecule_top = top;
      force_field = ff;
    }
    else {
      if (!settings.force_fields.includes(force_field)) {
        return Errors.throw(ErrorType.InvalidForceField);
      }
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

    if (!lipids_str || !force_field) {
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

      if (!items.every(e => !isNaN(e) && e >= 0)) {
        return Errors.throw(ErrorType.Format);
      }
      opts.box = items;
    }

    // Auto convert all items that are numbers
    for (const item of VALID_BODY_ITEMS) {
      if (item in req.body && req.body[item]) {
        opts[item] = convertOrThrow(req.body[item]);
      } 
    }

    // Handle rotate
    if (req.body.rotate && req.body.rotate !== 'none') {
      if (req.body.rotate === 'angle') {
        opts.rotate_angle = convertOrThrow(req.body.rotate_angle);
        opts.rotate = 'angle';
      }
      else {
        opts.rotate = req.body.rotate;
      }
    }

    // Handle booleans
    if (isBoolTrue(req.body.center)) {
      opts.center = true;
    }
    if (isBoolTrue(req.body.orient)) {
      opts.orient = true;
    }

    try {
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
    } catch (e) {
      logger.error('[INSANE] Insane run failed.');
      logger.error(e); 

      if (e instanceof InsaneError) {
        const dir = e.workdir;

        res
          .status(400)
          .json({
            error: e.message,
            trace: e.trace,
            zip: await Martinizer.zipDirectoryString(dir)
          });
      }
      else {
        logger.error(inspect(e));
        throw e;
      }
    }
  })().catch(errorCatcher(res));
});

MembraneBuilderRouter.all('/', methodNotAllowed(['POST']));

export default MembraneBuilderRouter;
