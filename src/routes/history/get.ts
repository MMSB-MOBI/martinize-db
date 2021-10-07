import { Router } from 'express';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher } from '../../helpers'

const GetHistoryRouter = Router()

GetHistoryRouter.get('/', (req, res) => {
    (async () => {
        console.log("POST GET")
        console.log(req)
        const jobId = req.query.jobId as string
        if (jobId) {
            const job = await HistoryOrganizer.getJob(jobId);
            //job.files = await HistoryOrganizer.readFiles(jobId, job.files); 

            job.files = await HistoryOrganizer.readFiles(jobId, job.files)
            
            res.json(job)
        }
        else throw new Error("server receives no job")


    })().catch(errorCatcher(res))
})

export default GetHistoryRouter