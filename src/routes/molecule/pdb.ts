import { Router } from 'express';
import { errorCatcher } from '../../helpers';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import Errors, { ErrorType } from '../../Errors';
import { Database } from '../../Entities/CouchHelper';
import { BaseMolecule } from '../../Entities/entities';
// @ts-ignore 
import NodeStreamZip from 'node-stream-zip';

// Get a pdb from a file ID
const PdbGetterRouter = Router();

async function getReadableStream(name: string, zip: NodeStreamZip) {
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

PdbGetterRouter.get('/:id', (req, res) => {
  (async () => {
    // Find/get file by ID
    const molecule = await MoleculeOrganizer.getInfo(req.params.id);

    console.log("/molecule/representation", molecule?.itp)
    // File does not exists or does not have a pdb file attached
    if (!molecule || !molecule.itp.length) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    // Find the molecule to find which force field is used.
    let molecule_data: BaseMolecule;
    const holders = await Database.molecule.find({ selector: { files: req.params.id } });
    if (holders.length) {
      molecule_data = holders[0];
    }
    else {
      const holders_stashed = await Database.stashed.find({ selector: { files: req.params.id } });
      molecule_data = holders_stashed[0];
    }

    if (!molecule_data) {
      return Errors.throw(ErrorType.MoleculeNotFound);
    }

    const itps: NodeJS.ReadableStream[] = [];

    const zip = new NodeStreamZip({
      file: MoleculeOrganizer.getFilenameFor(req.params.id),
      storeEntries: true
    });

    await new Promise((resolve, reject) => {
      zip.on('ready', resolve);
      zip.on('error', reject);
    });
    
    // Read every ITP in a stream
    // for (const itp of molecule.itp) {
    //   const stream = await new Promise((resolve, reject) => {
    //     zip.stream(itp.name, (err: any, stm: any) => {
    //       if (err)
    //         reject(err);
    //       resolve(stm);
    //     });
    //   }) as NodeJS.ReadableStream;
      
    //   itps.push(stream);
    // }

    // Read the PDB in a stream
    const pdb_stream = await getReadableStream(molecule.pdb!.name, zip)
    const top_stream = await getReadableStream(molecule.top!.name, zip)

    const itp_streams = await Promise.all(molecule.itp.map(itp_file => getReadableStream(itp_file.name, zip)))

    //console.log(itp_streams)

    //const itp_streams = await Promise.all( new Promise(()))   
    // Generate the needed radius
    // const [pdb_final, radius] = await Database.radius.transformPdb(pdb_stream, itps, molecule_data.force_field);


    const pdb_final = await Promise.all(readStreamsPromises([pdb_stream]))

    const top_final = await Promise.all(readStreamsPromises([top_stream]))

    const itp_finals = await Promise.all(readStreamsPromises(itp_streams))


    res.json({
      radius: {},
      pdb: pdb_final[0], 
      top: top_final[0],
      itps : itp_finals
    });
  })().catch(errorCatcher(res));
});

export default PdbGetterRouter;
