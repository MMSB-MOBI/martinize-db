import { Router } from 'express';
import { sendError } from '../../helpers';
import HistoryOrganizer from '../../HistoryOrganizer'
import Errors, { ApiError, ErrorType } from "../../Errors";

const GetHistoryRouter = Router()

GetHistoryRouter.get('/', async (req, res) => {
        const jobId = req.query.jobId as string
        if (jobId) {
          try {
            const job = await HistoryOrganizer.getJob(jobId)
            const files = await HistoryOrganizer.readFiles(jobId, job.files)
            job.files = files
            res.json(job)
          }catch(e){
            console.error(e)
            if (e.code && e.code === "ENOENT") {
                await HistoryOrganizer.deleteJob(jobId)
                sendError(Errors.make(ErrorType.HistoryFilesNotFound), res)
            }
            else sendError(Errors.make(ErrorType.Server, e), res)
          }
        }
        else sendError(Errors.make(ErrorType.JobNotProvided), res)
    })

export default GetHistoryRouter