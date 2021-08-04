import { History } from "./entities";
import AbstractDatabase from "./AbstractDatabase";
import logger from "../logger";
const util = require('util');

export default class HistoryDatabase extends AbstractDatabase<History>{

  async addToHistory(user_id: string, job_id:string){
    logger.debug(`HistoryDatabase : user ${user_id}`)
    const exists = await this.exists(user_id);
    if (!exists){
      logger.debug("HistoryDatabase : create new user history")
      const doc: History = {
        id: user_id, 
        job_ids: [job_id]
      }
      return this.save(doc)
    }
    else {
      logger.debug("HistoryDatabase : add job to already existing history")
      const currentHistory = await this.get(user_id)
      currentHistory.job_ids.push(job_id)
      
      return this.save(currentHistory)
      
    }
  }

}
