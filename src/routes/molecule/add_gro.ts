import { readFileSync } from "fs";
import md5File from "md5-file/promise";
import { cpuUsage, off } from "process";
import { Martinizer } from "../../Builders/Martinizer";
import { Database } from "../../Entities/CouchHelper";
import MoleculeOrganizer from "../../MoleculeOrganizer";
import { InfosJson } from "./CreateMoleculeJson";
import StreamZip from 'node-stream-zip';
import { getReadableStream } from "./getMoleculeAPI";
import TmpDirHelper from "../../TmpDirHelper";

export const addGroFromBatch = async (batch : InfosJson[]) => {
    const tmpDir = await TmpDirHelper.get()
    for (const molecule of batch){
        await addOneGro(molecule, tmpDir)
    }
}

const addOneGro = async(molecule : InfosJson, tmpDir : string) => {
    for (const version of molecule.versions) {
        console.log('#', molecule.alias, version.force_field)
        const query = {selector : {alias : molecule.alias, force_field:version.force_field}}
        const inDb = await Database.molecule.find(query)
        console.log(version)
        const currentItpMd5 = await md5File(version.itp.path)
        console.log("itp md5", currentItpMd5)

        //const { pdb_stream } = await Martinizer.createPdbWithConect(version.gro.path, readFileSync(version.top.path).toString(), false, version.force_field)

        //console.log(pdb_stream); 
        
        for (const version2 of inDb) {
            console.log('# Entry', version2.alias, version2.force_field, version2.version)
            const filesInfo = await MoleculeOrganizer.getInfo(inDb[0].files)
            console.log(filesInfo)

            const zipName = await MoleculeOrganizer.getFilenameFor(inDb[0].files)
            
            const zip = new StreamZip.async({ file: zipName });
  
            await zip.extract(filesInfo!.itp[0].name, tmpDir + '/test.itp')
            await zip.close()

            const dbItpMd5 = await md5File(tmpDir + '/test.itp')

            console.log('db itp md5', dbItpMd5)

            

        }

        
        
        

    }    
    
}