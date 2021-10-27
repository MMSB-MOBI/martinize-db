import { Router } from 'express'
import { Job } from '../../Entities/entities'
import Uploader from '../Uploader';
import multer from 'multer';
import HistoryOrganizer from '../../HistoryOrganizer';
import { sendError, errorCatcher } from '../../helpers';
import Errors, { ApiError, ErrorType } from "../../Errors";

const UpdateHistoryRouter = Router(); 

UpdateHistoryRouter.post('/', Uploader.array('files',99), (req, res) => {
    (async () => {
        const jobId = req.body.jobId
        const files = req.files as Express.Multer.File[]

        if (!jobId || !files){
            return Errors.throw(ErrorType.MissingParameters);
        }

        await HistoryOrganizer.updateJobInFileSystem(jobId, files)
        await HistoryOrganizer.updateJobForSavedBonds(jobId, files)

        res.json({jobId: 'updated'}); 

    })().catch(errorCatcher(res))
})

export default UpdateHistoryRouter