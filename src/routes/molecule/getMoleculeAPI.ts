import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { BaseMolecule } from '../../Entities/entities';
// @ts-ignore 
import StreamZip from 'node-stream-zip';
import JSZip from 'jszip';
import { MangoQuery } from 'nano';
import logger from '../../logger';
import { inspect } from 'util';
// Get a pdb from a file ID
const GetMoleculeAPI = Router();

export async function getReadableStream(name: string, zip: StreamZip) {
  return new Promise((resolve, reject) => {
    zip.stream(name, (err: any, stm: any) => {
      if (err)
        reject(err);
      resolve(stm);
    });
  }) as Promise<NodeJS.ReadableStream>
}

function readStreamsPromises(streams: NodeJS.ReadableStream[]): Promise<string>[] {
  let readedStream = "";
  return streams.map(stream => new Promise((resolve, reject) => {
    stream.on('data', (chunk: string) => {
      readedStream += chunk
    })
    stream.on('end', () => resolve(readedStream))
    stream.on('error', reject)
  }))
}


//If format isnt provided give the last update of this model 
GetMoleculeAPI.get('/:forcefield/:id.:format?/:version?', (req, res) => {
  (async () => {
    logger.info(`[GetMoleculeAPI] ${inspect(req.params)}`);
    const selectruc: MangoQuery = { selector: { alias: req.params.id, force_field: req.params.forcefield } }
    if (req.params.version) {
      selectruc.selector["version"] = req.params.version
    }

    const molcouch = await Database.molecule.find(selectruc)


    // File does not exists
    if (molcouch.length === 0) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const molecule = await MoleculeOrganizer.getInfo(molcouch[0].files);

    logger.info(`[GetMoleculeAPI] molcouch[0].files ${molcouch[0].files}`);
    const zip = new StreamZip({
      file: MoleculeOrganizer.getFilenameFor(molcouch[0].files),
      storeEntries: true
    });
    
    await new Promise((resolve, reject) => {
      zip.on('ready', ()=> resolve);
      zip.on('error', reject);
    });
    
    if (req.params.format === "itp") {
      const itp_streams = await Promise.all(molecule!.itp.map(itp_file => getReadableStream(itp_file.name, zip)))
      const itp_finals = await Promise.all(readStreamsPromises(itp_streams))
      logger.info(`[GetMoleculeAPI]::itp Sending ${itp_finals[0]}`);
      res.send(itp_finals[0]);
    }
    else if (req.params.format === "pdb") {
      const pdb_stream = await getReadableStream(molecule!.pdb!.name, zip)
      const pdb_final = await Promise.all(readStreamsPromises([pdb_stream]))
      logger.info(`[GetMoleculeAPI]::pdb Sending ${pdb_final[0]}`);
      res.send(pdb_final[0]);
    }
    else if (req.params.format === "gro") {
      if (molecule!.gro) {
        const gro_stream = await getReadableStream(molecule!.gro!.name, zip)
        const gro_final = await Promise.all(readStreamsPromises([gro_stream]))
        logger.info(`[GetMoleculeAPI]::pdb Sending ${gro_final[0]}`);
        res.send(gro_final[0]);
      }
      else {
        logger.error(`[GetMoleculeAPI]::gro \".gro not found\"`);
        res.status(404).json({ "error": ".gro not found." });
      }

    }
    else if ((req.params.format === undefined) || (req.params.format === "zip")) {
      const filename = MoleculeOrganizer.getFilenameFor(molcouch[0].files);
      res.download(filename);
    }
    else {
      logger.error(`[GetMoleculeAPI] Format .${req.params.format} unkown.`);
      res.send({ "error": "Format ." + req.params.format + " unkown." })
    }
  })().catch(errorCatcher(res));
});

export default GetMoleculeAPI;