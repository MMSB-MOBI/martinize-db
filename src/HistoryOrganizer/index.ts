import logger from '../logger';
import { HISTORY_ROOT_DIR } from '../constants';
import fs, { promises as FsPromise } from 'fs';
import path, { resolve } from 'path';
import { Database } from '../Entities/CouchHelper';
import { Job } from '../Entities/entities'
import { getFormattedFile } from "../helpers"; 
import { isCouchNotFound, notFoundOnFileSystem } from '../Errors';
import { JobFilesNames, JobReadedFiles } from '../types'

function getUserJobObject(jobsDoc:Job[]){
    let obj : {[userId: string]: string[]} = {}
    jobsDoc.forEach(job => {
        if (!obj.hasOwnProperty(job.userId)) obj[job.userId] = [job.id]
        else obj[job.userId].push(job.id)
    })
    return obj
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
            logger.verbose(`copy ${file}`)
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

    async updateJobForSavedBonds(jobId: string, itp_files_names: string[][]){
        return await Database.job.updateManuallySavedBonds(jobId, itp_files_names); 
    }

    async deleteFromFileSystem(jobId: string){
        logger.debug(`Delete ${jobId} from file system`)
        const dirPath = HISTORY_ROOT_DIR + "/" + jobId
        await FsPromise.rmdir(dirPath, {recursive: true}); 
    }

    async deleteMultipleFromFileSystem(jobIds : string []){
        logger.debug(`Delete ${jobIds} from file system`)
        return await Promise.all(jobIds.map(id => FsPromise.rmdir(HISTORY_ROOT_DIR + "/" + id, {recursive:true})))
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

    async deleteMultipleFromCouch(jobIds: string[]){
        logger.debug(`Delete multipe ${jobIds} from couch`)
        const jobs = await Database.job.bulkGet(jobIds)
        let notFoundIdx: number[] = []; 
        const filteredJobs = jobs.filter((job, idx) => {
            if (job === null){
                notFoundIdx.push(idx); 
                return false
            } 
            return true
        })
       
        const usersRelatedToJobs = getUserJobObject(filteredJobs)
        let promises = jobs.map(job => Database.job.delete(job))
        for (const [user, jobIds] of Object.entries(usersRelatedToJobs)){
            promises.push(Database.history.deleteJobs(user, jobIds))
        } 
        return await Promise.all(promises); 
        
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
                    itp_files : await Promise.all(files.itp_files.map(async mol_itp => await Promise.all(mol_itp.map(i => getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${i}`))))), 
                    warnings: await getFormattedFile(`${HISTORY_ROOT_DIR}/${jobId}/${files.warnings}`),
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
    
    async deleteJobs(jobIds : string[]){    
        return Promise.all([this.deleteMultipleFromCouch(jobIds), this.deleteMultipleFromFileSystem(jobIds)])

    }

    
    
}


export default HistoryOrganizer;