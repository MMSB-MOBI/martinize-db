import AbstractDatabase from "./AbstractDatabase"
import logger from "../logger";


export default class JobDatabase extends AbstractDatabase<any> {

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

    async getJobsDetails(jobIds : string[]){
        const jobsDetails = await this.bulkGet(jobIds)
        return jobsDetails
    }
}