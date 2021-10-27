import AbstractDatabase from "./AbstractDatabase"
import logger from "../logger";
import { Database } from "./CouchHelper";
import { Job } from './entities'


export default class JobDatabase extends AbstractDatabase<Job> {

    async addToJob(jobInfos : any){
        logger.debug(`JobDatabase : add job ${jobInfos.id}`)
        const exists = await this.exists(jobInfos.id); 
        if (!exists) {
            logger.debug("JobDatabase : create new job entry")
            return this.save(jobInfos)
        }
        else {
            logger.error("This job already exists. Should not happen")
            throw new Error("Job already exists")
        }
    }

    async getJobsDetails(jobIds : string[], userId: string){
            const jobsDetails = await this.bulkGet(jobIds)
            const jobsFound = jobsDetails.filter(job => job !== null).map(job => job.id)
            const notFound = jobIds.filter(x => !jobsFound.includes(x));
            if (notFound.length > 0) {
                logger.warn(`job(s) ${notFound} not found in job database`)
                //Clean
                await Database.history.deleteJobs(userId, notFound);
                
            }
            return jobsDetails.filter(job => job !== null); 
        }


    async updateManuallySavedBonds(id: string, newItpFiles: string[]) {
        const updateFnc = (doc : Job) => {

            const newItps = newItpFiles.filter(itp => ! doc.files.itp_files.includes(itp))
            doc.files.itp_files = [...doc.files.itp_files, ...newItps]

            if(doc.manual_bonds_edition) return doc
            doc.manual_bonds_edition = true 
            return doc
        }
        return await this.update(id, updateFnc) 
    }
        

}