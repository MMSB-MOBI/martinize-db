
import { Router } from 'express'; 
import { cpuUsage } from 'process';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher, sendError } from '../../helpers'
import Errors, { ApiError, ErrorType } from "../../Errors";
import logger from '../../logger';


const ListHistoryRouter = Router(); 

ListHistoryRouter.get('/', async (req, res) => {
    const user = req.user?.user_id; 
    if(user) {
        HistoryOrganizer.getHistory(user).then(jobs => res.json(jobs.reverse()))
        .catch(e => {
            logger.error("Error while get user history", e)
            console.log(e)
            if (e.error && e.error === "not_found" && e.reason && e.reason === "deleted") sendError(Errors.make(ErrorType.HistoryNotFound), res)
            else sendError(Errors.make(ErrorType.Server), res)
            
        })
    }
    else sendError(Errors.make(ErrorType.UserNotProvided), res)
        
})

export default ListHistoryRouter