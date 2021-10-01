import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path from 'path';
import { Database } from '../Entities/CouchHelper';
import { Job } from "../Entities/entities"


interface MartinizeFiles {

}

export const HistoryOrganizer = new class HistoryOrganizer{
    constructor(){
        logger.info("HistoryOrganizer constructor")
        try{
            fs.mkdirSync(HISTORY_ROOT_DIR)
        }catch{}
        
    }

    async _saveToFileSystem(job_id: string, files: string[]):Promise<any> {
        logger.debug("HistoryOrganizer save")
        const dirPath = HISTORY_ROOT_DIR + "/" + job_id
        fs.mkdirSync(dirPath)
        files.forEach(file => {
            logger.silly(`copy ${file}`)
            const name = path.basename(file); 
            fs.copyFileSync(file, dirPath + "/" + name)
        })
    }

    async save(jobInfos: Job, files: string[], userId: string){
        await this._saveToFileSystem(jobInfos.id, files)
        await Database.job.addToJob(jobInfos)
        await Database.history.addToHistory(userId, jobInfos.id)
    }

    async getHistory(userId: string){
        const jobIds = await Database.history.getAllJobs(userId)
        return await Database.job.getJobsDetails(jobIds)
        
    }
}

export default HistoryOrganizer;