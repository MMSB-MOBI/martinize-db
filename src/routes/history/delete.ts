import { Router } from 'express';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher } from '../../helpers'

const DeleteHistoryRouter = Router()

DeleteHistoryRouter.get('/', (req, res) => {
    (async () => {
        if (!req.query.jobIds){
            throw new Error("server receives no jobs")
        }
        const jobIds = (req.query.jobIds as string).split(",")
        const resp = await HistoryOrganizer.deleteJobs(jobIds)
        res.json(resp); 

    })().catch(errorCatcher(res))
})

export default DeleteHistoryRouter