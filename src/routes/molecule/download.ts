import { Router } from 'express';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import MoleculeOrganizer from '../../MoleculeOrganizer';
import JSZip from 'jszip';
var fs = require("fs");


const DownloadMoleculeRouter = Router();



DownloadMoleculeRouter.get('/', (req, res) => {
  (async () => {
    if(req.query.id){
      const { id, filename } = req.query as Record<string, string>;
      if (!id) {
        return Errors.throw(ErrorType.MissingParameters);
      }
  
      try {
        BigInt(id);
      } catch {
        return Errors.throw(ErrorType.Format);
      }
  
      if (!(await MoleculeOrganizer.exists(id))) {
        return Errors.throw(ErrorType.ElementNotFound);
      }
  
      // todo check filename
      res.download(MoleculeOrganizer.getFilenameFor(id), filename);
    }
    else {
      //@ts-ignore
      let molecule;
      const zip = new JSZip();
      
      for(let champ in req.query){
        //@ts-ignore
        let molecules = req.query[champ].split(/;,|;/);
        let mol: string;
        for(mol of molecules){
          mol = mol.replace(/'/gi, '"')
          console.log(mol)
          mol === '' ? undefined : molecule = JSON.parse(mol);

          molecule.filenameFor = MoleculeOrganizer.getFilenameFor(molecule.id);

          zip.file(molecule.filenameFor);
        }
      }
    zip.generateAsync({type: "binarystring"}).then((z) => res.download(z, "molecules.zip"))
    }
  })().catch(errorCatcher(res));
});

DownloadMoleculeRouter.all('/', methodNotAllowed('GET'))

export default DownloadMoleculeRouter;
