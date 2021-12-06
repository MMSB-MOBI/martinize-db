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
        const jobId = req.body.jobId as string; 
        const files = req.files as Express.Multer.File[]    
        const molIdxs = req.body.molIdxs as number[]; 

        if (!jobId || !files || !molIdxs ){
            return Errors.throw(ErrorType.MissingParameters);
        }
        
        const organizedFileNames: string[][]  = []
        for (const [index, value] of molIdxs.entries()){
            if(!organizedFileNames[value]) organizedFileNames[value] = []
            organizedFileNames[value].push(files[index].originalname)
        }
        
        await HistoryOrganizer.updateJobInFileSystem(jobId, files)
        await HistoryOrganizer.updateJobForSavedBonds(jobId, organizedFileNames)

        res.json({jobId: 'updated'}); 

    })().catch(errorCatcher(res))
})

export default UpdateHistoryRouter