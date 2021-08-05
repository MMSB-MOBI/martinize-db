import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path from 'path';


interface MartinizeFiles {

}

export const HistoryOrganizer = new class HistoryOrganizer{
    constructor(){
        logger.info("HistoryOrganizer constructor")
        try{
            fs.mkdirSync(HISTORY_ROOT_DIR)
        }catch{}
    }

    async save(job_id: string, files: string[]):Promise<any> {
        logger.debug("HistoryOrganizer save")
        const dirPath = HISTORY_ROOT_DIR + "/" + job_id
        fs.mkdirSync(dirPath)
        files.forEach(file => {
            logger.silly(`copy ${file}`)
            const name = path.basename(file); 
            fs.copyFileSync(file, dirPath + "/" + name)
        })
    }
}

export default HistoryOrganizer;