
import { Router } from 'express'; 
import { cpuUsage } from 'process';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher, sendError } from '../../helpers'
import Errors, { ApiError, ErrorType } from "../../Errors";


const ListHistoryRouter = Router(); 

ListHistoryRouter.get('/', async (req, res) => {
    const user = req.user?.user_id; 
    if(user) {
        HistoryOrganizer.getHistory(user).then(jobs => res.json(jobs))
        .catch(e => {
            if (e.error && e.error === "not_found") sendError(Errors.make(ErrorType.HistoryNotFound), res)
            else sendError(Errors.make(ErrorType.Server), res)
            
        })
    }
    else sendError(Errors.make(ErrorType.UserNotProvided), res)
        
})

export default ListHistoryRouter