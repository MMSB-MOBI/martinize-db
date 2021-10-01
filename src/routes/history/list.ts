
import { Router } from 'express'; 
import { cpuUsage } from 'process';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher } from '../../helpers'

const ListHistoryRouter = Router(); 

ListHistoryRouter.get('/', (req, res) => {
    (async () => {
        const user = req.user?.user_id; 
        if(user) {
            const jobs = await HistoryOrganizer.getHistory(user); 
            console.log("jobs")
            console.log(jobs)
            res.json(jobs)
        }
        else throw new Error("no user given in history/list request")
        
    })().catch(errorCatcher(res))
})

export default ListHistoryRouter