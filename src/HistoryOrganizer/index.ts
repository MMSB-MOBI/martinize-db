import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path, { resolve } from 'path';
import { Database } from '../Entities/CouchHelper';
import { getFormattedFile } from "../helpers"; 

interface MartinizeFileNames {
    all_atom : string; 
    coarse_grained: string; 
    itp_files : string[]; 
    top_file : string; 
}

export const HistoryOrganizer = new class HistoryOrganizer{
    constructor(){
        logger.info("HistoryOrganizer constructor")
        try{
            fs.mkdirSync(HISTORY_ROOT_DIR)
        }catch{}
        
    }

    async _saveToFileSystem(job_id: string, files: string[]):Promise<any> {
        const dirPath = HISTORY_ROOT_DIR + "/" + job_id
        fs.mkdirSync(dirPath)
        files.forEach(file => {
            logger.silly(`copy ${file}`)
            const name = path.basename(file); 
            fs.copyFileSync(file, dirPath + "/" + name)
        })
    }

    async _deleteFromFileSystem(jobId: string){
        logger.debug(`Delete ${jobId} from file system`)
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        await FsPromise.rmdir(dirPath, {recursive: true}); 
    }

    async _saveToCouch(job: any){
        const jobDoc = {id : job.jobId, ...job}
        await Database.job.addToJob(jobDoc)
        await Database.history.addToHistory(job.userId, job.jobId)
    }

    async saveToHistory(job : any, files: string[]){
        return await Promise.all([this._saveToFileSystem(job.jobId, files), this._saveToCouch(job)]) 
    }

    async getHistory(userId: string){
        const jobIds = await Database.history.getAllJobs(userId)
        return await Database.job.getJobsDetails(jobIds, userId) 
    }

    async getJob(jobId: string) {
        return await Database.job.get(jobId)
    }

    async readFiles(jobId: string, files:MartinizeFileNames){
        const readedFiles = {
            all_atom : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.all_atom}`),
            top_file : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.top_file}`),
            coarse_grained : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.coarse_grained}`),
            itp_files : await Promise.all(files.itp_files.map(i => getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${i}`)))
        }

        return readedFiles
    }
    
    async deleteJob(jobId : string){
        const job = await Database.job.get(jobId)
        const user = job.userId
        const resp = await Database.job.delete(job)
        await Database.history.deleteJobs(user, [jobId]); 
        await this._deleteFromFileSystem(jobId); 
        return {'deleted' : jobId}
    }
    
}


export default HistoryOrganizer;