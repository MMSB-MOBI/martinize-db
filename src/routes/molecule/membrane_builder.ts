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
import { plainToInstance } from 'class-transformer';
import { ClientInsaneSettingsDto, FileDto } from './membrane_builder.dto';
import { validateOrReject } from 'class-validator';
import { AvailableForceFields } from '../types';

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
    
    // Init
    console.log("files", req.files);
    const validatedParams = plainToInstance(ClientInsaneSettingsDto, req.body);

    //Validate file names
    const files = req.files as { [fieldname: string]: Express.Multer.File[] };
    if(Object.keys(files).length > 0) {
      if(files.pdb.length > 1) return Errors.throw(ErrorType.TooManyFiles)
      if(files.top.length > 1) return Errors.throw(ErrorType.TooManyFiles)
      const validatedPdb = plainToInstance(FileDto, files.pdb[0])
      const validatedTop = plainToInstance(FileDto, files.top[0])
      const validatedItps = files.itp.map(itp => plainToInstance(FileDto, itp))
  
      try {
        await validateOrReject(validatedParams); 
        await validateOrReject(validatedPdb); 
        await validateOrReject(validatedTop); 
        await Promise.all(validatedItps.map(dto => validateOrReject(dto)))
      } catch(e) {
        res.status(400).json({ error: true, statusCode: 400, errorCode: 'PARAMS_VALIDATION_ERROR', e })
        return; 
      }
    }


    //const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));
     // Maybe security check
    const molecule_id = validatedParams.from_id
    const { 
      pbc, 
      box, 
      lipids: lipids_str, 
      upper_leaflet: upper_leaflet_str 
    } = validatedParams
    let force_field = validatedParams.force_field


    const opts: Partial<InsaneSettings> = {};

    // Test force field
    
    const molecule_entries = {
      molecule_pdb: "",
      molecule_top: "",
      molecule_itps: [] as string[],
    }; // m

    if (molecule_id) {
      // from molecule id
      const { pdb, itps, top, force_field: ff } = await MembraneBuilder.prepareRunWithDatabaseMolecule(molecule_id.toString());
      molecule_entries.molecule_itps = itps;
      molecule_entries.molecule_pdb = pdb;
      molecule_entries.molecule_top = top;
      force_field = ff as AvailableForceFields;
    } // m
    else {
      if (validatedParams.molecule_added) {
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
      
    }

    let upper_leaflet: LipidMap = [];
    let lipids = undefined;
    if (validatedParams.lipids_added) {
      if (!lipids_str || !force_field) {
        return Errors.throw(ErrorType.MissingParameters);
      } // l  

      // Parse lipid str
      lipids = (lipids_str as string).split(',').map(e => { 
        const res = e.split(':'); 
        
        if (res.length > 1) {
          return [res[0], parseInt(res[1], 10)] as [string, number];
        }
        return [res[0], 1] as [string, number];
      }); // l

      // If upper leaflet, parse it
      if (upper_leaflet_str) {
        upper_leaflet = (upper_leaflet_str as string).split(',').map(e => { 
          const res = e.split(':'); 
          
          if (res.length > 1) {
            return [res[0], parseInt(res[1], 10)] as [string, number];
          }
          return [res[0], 1] as [string, number];
        });
      } // l
    }

    // Parse settings
    opts.pbc = pbc;
    const items = (box as string).split(',').map(e => parseInt(e, 10));
    if (!items.every(e => !isNaN(e) && e >= 0)) {
      return Errors.throw(ErrorType.Format);
    }
    opts.box = items;
    
    // Handle rotate
    if (validatedParams.rotate !== 'none') {
      if (validatedParams.rotate === 'angle') {
        opts.rotate_angle = validatedParams.rotate_angle;
        opts.rotate = 'angle';
      }
      else {
        opts.rotate = validatedParams.rotate;
      }
    }

    if (validatedParams.molecule_added) {
      if (validatedParams.center) {
        opts.center = true;
      }
    }
    if (validatedParams.lipids_added && validatedParams.molecule_added) {
      if (validatedParams.orient) {
        opts.orient = true;
      }
    }

    opts.salt_concentration = validatedParams.salt_concentration;
    if(validatedParams.charge !== 0){
      opts.charge = validatedParams.charge;
    }
    opts.solvent_type = validatedParams.solvent_type;
    
    if(!force_field){
      return Errors.throw(ErrorType.MissingParameters)
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
        // @ts-ignore
        water: await getFormattedFile(water),
        no_water: await getFormattedFile(no_water),
        top: await getFormattedFile(top),
        itps: await Promise.all(itps.map(i => getFormattedFile(i))),
      });
    } catch (e) {
      console.log(e)
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
  })().catch(
    errorCatcher(res)
  );
});

MembraneBuilderRouter.all('/', methodNotAllowed(['POST']));

export default MembraneBuilderRouter;
