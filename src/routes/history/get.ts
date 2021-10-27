import { Router } from 'express';
import { sendError } from '../../helpers';
import HistoryOrganizer from '../../HistoryOrganizer'
import Errors, { ApiError, ErrorType } from "../../Errors";
import { ReadedJob } from '../../types';

const GetHistoryRouter = Router()

GetHistoryRouter.get('/', async (req, res) => {
        const jobId = req.query.jobId as string
        if (jobId) {
          HistoryOrganizer.getJob(jobId).then(job => {
            const {files, ...jobBase} = job
            HistoryOrganizer.readFiles(jobId, files).then(readedFiles => {
              const readedJob : ReadedJob = {
                ...jobBase, 
                files : readedFiles
              }
              res.json(readedJob)
            }).catch(e => {
              if (e === "not_found") sendError(Errors.make(ErrorType.HistoryFilesNotFound), res)
              else sendError(Errors.make(ErrorType.Server, e), res)
            })
          }).catch(e => {
            if (e === "not_found") sendError(Errors.make(ErrorType.JobNotFound), res)
            else sendError(Errors.make(ErrorType.Server, e), res)
          })
        }
        else sendError(Errors.make(ErrorType.JobNotProvided), res)
    })

export default GetHistoryRouter