import { Router } from 'express';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher } from '../../helpers'

const DeleteHistoryRouter = Router()

DeleteHistoryRouter.get('/', (req, res) => {
    (async () => {
        const jobId = req.query.jobId as string
        if (jobId) {
            const resp = await HistoryOrganizer.deleteJob(jobId)
            res.json(resp); 
        }
        else throw new Error("server receives no job")


    })().catch(errorCatcher(res))
})

export default DeleteHistoryRouter