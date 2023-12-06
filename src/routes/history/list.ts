
import { Router } from 'express'; 
import { cpuUsage } from 'process';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher, sendError } from '../../helpers'
import Errors, { ApiError, ErrorType } from "../../Errors";
import logger from '../../logger';
import { inspect } from 'util';

const ListHistoryRouter = Router(); 

ListHistoryRouter.get('/', async (req, res) => {
    logger.info("[Router:history::list] Hello !!");
    const user = req.query?.user; 
    if(user) {
        HistoryOrganizer.getHistory(user as string).then(jobs => res.json(jobs.reverse()))
        .catch(e => {
            logger.error(`[Router:history::list] Error while get user history : ${e}`);
            if (e.error && e.error === "not_found") 
                sendError(Errors.make(ErrorType.HistoryNotFound), res);
            else 
                sendError(Errors.make(ErrorType.Server), res);
            
        })
    }
    else {
        logger.error(`[Router:history::list] UserNotProvided , \"req.user?.user_id\" missing from\n${inspect(req)}`);
        sendError(Errors.make(ErrorType.UserNotProvided), res)
    }
        
})

export default ListHistoryRouter