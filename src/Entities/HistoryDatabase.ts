import { History } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import logger from "../logger";
const util = require('util');


export default class HistoryDatabase extends AbstractDatabase<History>{

  async addToHistory(user_id: string, jobId : string){
    logger.debug(`HistoryDatabase : user ${user_id} and job ${jobId}`)

    const exists = await this.exists(user_id);
    if (!exists){
      logger.debug("HistoryDatabase : create new user history")
      const doc: History = {
        id: user_id, 
        job_ids: [jobId]
      }
      return this.save(doc)
    }
    else {
      logger.debug("HistoryDatabase : add job to already existing history")
      const currentHistory = await this.get(user_id)
      currentHistory.job_ids.push(jobId)
      return this.save(currentHistory)
      
    }
  }

  async getAllJobs(user_id: string) {
    logger.debug(`History database : get all jobs for user ${user_id}`)
    const jobs : History = await this.get(user_id)
    logger.debug(`${jobs.job_ids.length} jobs found`)
    return jobs.job_ids
  }

}


