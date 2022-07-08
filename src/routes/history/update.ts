import { Router } from 'express'
import Uploader from '../Uploader';
import HistoryOrganizer from '../../HistoryOrganizer';
import { errorCatcher } from '../../helpers';
import Errors, { ErrorType } from "../../Errors";

const UpdateHistoryRouter = Router(); 

UpdateHistoryRouter.post('/', Uploader.array('files',99), (req, res) => {
    (async () => {
        const jobId = req.body.jobId as string; 
        const files = req.files as Express.Multer.File[]    
        const molIdxs = req.body.molIdxs as number[]; 
        const editionComment = req.body.editionComment as string; 

        if (!jobId || !files || !molIdxs){
            return Errors.throw(ErrorType.MissingParameters);
        }
        
        const organizedFileNames: string[][]  = []
        for (const [index, value] of molIdxs.entries()){
            if(!organizedFileNames[value]) organizedFileNames[value] = []
            organizedFileNames[value].push(files[index].originalname)
        }
        
        const newId = await HistoryOrganizer.updateJobInFileSystem(jobId, files)
        await HistoryOrganizer.updateJobAndCreateANewOne(jobId, newId, organizedFileNames, editionComment)

        res.json({jobId: 'updated'}); 

    })().catch(errorCatcher(res))
})

export default UpdateHistoryRouter