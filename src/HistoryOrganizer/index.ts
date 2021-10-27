import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path, { resolve } from 'path';
import { Database } from '../Entities/CouchHelper';
import { Job } from '../Entities/entities'
import { getFormattedFile } from "../helpers"; 
import { isCouchNotFound, notFoundOnFileSystem } from '../Errors';
import { JobFilesNames, JobReadedFiles } from '../types'


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

    async updateJobInFileSystem(jobId: string, itp_files:Express.Multer.File[]){
        const jobDir = HISTORY_ROOT_DIR + "/" + jobId; 
        if (!fs.existsSync(jobDir)){
            throw new Error("Job directory doesn't exist")
        }
        logger.debug(`move new itp files to job directory ${jobId}`)
        return await Promise.all(itp_files.map(async(file) => {
            await FsPromise.rename(file.path, jobDir + "/" + file.originalname)
        }))

    }

    async updateJobForSavedBonds(jobId: string, itp_files: Express.Multer.File[]){
        const fileNames = itp_files.map(file => file.originalname)
        return await Database.job.updateManuallySavedBonds(jobId, fileNames); 
    }

    async deleteFromFileSystem(jobId: string){
        logger.debug(`Delete ${jobId} from file system`)
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        await FsPromise.rmdir(dirPath, {recursive: true}); 
    }

    async _deleteFromFileSystemIfExists(jobId : string){
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        try {
            await FsPromise.rmdir(dirPath, {recursive:true})
            logger.debug(`${jobId} deleted from file system`)
        }catch(e){
           if (notFoundOnFileSystem(e)) logger.debug(`job ${jobId} doesn't exist on file system, no deletion`)
           else throw(e)
        }
    }

    async deleteFromCouch(jobId: string){
        const job = await Database.job.get(jobId)
        const user = job.userId
        return await Promise.all([Database.job.delete(job), Database.history.deleteJobs(user, [jobId])]); 
    }

    async _deleteFromCouchIfExists(jobId : string){
        try {
            const job = await Database.job.get(jobId); 
            const user = job.userId
            await Promise.all([Database.job.delete(job), Database.history.deleteJobs(user, [jobId])])
            logger.debug(`${jobId} deleted from couch`)
        } catch(e) {
            if (isCouchNotFound(e)) logger.debug(`job ${jobId} doesn't exist on couch, no deletion`)
            else throw(e)
        }
    }

    async _saveToCouch(job: any){
        const jobDoc = {id : job.jobId, ...job}
        await Database.job.addToJob(jobDoc)
        await Database.history.addToHistory(job.userId, job.jobId)
    }

    async saveToHistory(job : any, files: string[]){
        return new Promise((res, rej) => {
            this._saveToFileSystem(job.jobId, files).then(() => {
                this._saveToCouch(job).then(() => res(job.jobId)).catch(e => {
                    this.deleteFromFileSystem(job.jobId)
                    rej(e)
                })
            }).catch(e => rej(e))
        })
        
    }

    async getHistory(userId: string){
        const jobIds = await Database.history.getAllJobs(userId)
        return await Database.job.getJobsDetails(jobIds, userId) 
    }

    async getJob(jobId: string) : Promise<Job> {
        return new Promise( async (res, rej) => {
            try {
                const job = await Database.job.get(jobId)
                res(job)
            } catch(e){
                if(isCouchNotFound(e)) { 
                    this._deleteFromFileSystemIfExists(jobId)
                    rej("not_found")
                }
                else rej(e)
            }
        })
        
    }

    async readFiles(jobId: string, files:JobFilesNames): Promise<JobReadedFiles>{
        return new Promise(async (res, rej) => {
            try {
                const readedFiles = {
                    all_atom : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.all_atom}`),
                    top_file : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.top_file}`),
                    coarse_grained : await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.coarse_grained}`),
                    itp_files : await Promise.all(files.itp_files.map(i => getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${i}`)))
                }
                res(readedFiles)
            } catch(e) {
                if (notFoundOnFileSystem(e)){
                    this._deleteFromCouchIfExists(jobId)
                    rej("not_found")
                } 
                else rej(e)
            }
        })
    }
    
    async deleteJob(jobId : string){
        return Promise.all([this.deleteFromCouch(jobId), this.deleteFromFileSystem(jobId)])

    }

    
    
}


export default HistoryOrganizer;