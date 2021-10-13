import { Router } from 'express';
import { sendError } from '../../helpers';
import HistoryOrganizer from '../../HistoryOrganizer'
import Errors, { ApiError, ErrorType } from "../../Errors";

const GetHistoryRouter = Router()

GetHistoryRouter.get('/', async (req, res) => {
        const jobId = req.query.jobId as string
        if (jobId) {
           HistoryOrganizer.getJob(jobId).then(async(job) => {
                HistoryOrganizer.readFiles(jobId, job.files)
                    .then(files => {
                        job.files = files
                        res.json(job)
                    }).catch(e => {
                        if (e.code && e.code === "ENOENT") sendError(Errors.make(ErrorType.HistoryFilesNotFound), res)
                        else sendError(Errors.make(ErrorType.Server, e), res)
                    })
                
           })
           .catch(e => sendError(Errors.make(ErrorType.Server,e), res));
        }
        else sendError(Errors.make(ErrorType.JobNotProvided), res)


    })

export default GetHistoryRouter