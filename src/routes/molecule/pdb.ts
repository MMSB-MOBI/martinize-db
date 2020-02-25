import { Router } from 'express';
import { errorCatcher, sleep } from '../../helpers';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import Errors, { ErrorType } from '../../Errors';
import NodeStream from 'stream';

// Get a pdb from a file ID
const PdbGetterRouter = Router();

PdbGetterRouter.get('/:id.pdb', (req, res) => {
  (async () => {
    // Find/get file by ID
    const molecule = await MoleculeOrganizer.get(req.params.id);

    // File does not exists or does not have a pdb file attached
    if (!molecule || !molecule[1].pdb) {
      return Errors.throw(ErrorType.ElementNotFound);
    }

    const [zip, save_info] = molecule;

    // Extract the pdb file from the internal ZIP file
    const file_as_buffer = await zip.file(save_info.pdb!.name).async("nodebuffer");

    // Create the write stream
    const write_stream = new NodeStream.PassThrough;

    res.set('Content-Type', 'chemical/x-pdb');
    // Must set for download progress estimation
    res.set('Content-Length', file_as_buffer.length.toString()); 

    // Attach write stream to request
    write_stream.pipe(res);

    // Send slowly the PDB to client
    const PARTS_NUMBER = 5;
    const PART_LENGTH = Math.ceil(file_as_buffer.length / PARTS_NUMBER);
    
    for (let i = 0; i < PARTS_NUMBER; i++) {
      const start = i * PART_LENGTH;
      const part = file_as_buffer.slice(start, start + PART_LENGTH);

      write_stream.push(part);
      await sleep(350);
    }

    // End stream then req
    write_stream.end("");
  })().catch(errorCatcher(res));
});

export default PdbGetterRouter;
