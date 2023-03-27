
import { Router } from 'express';
import { cpuUsage } from 'process';
import HistoryOrganizer from '../../HistoryOrganizer'
import { errorCatcher, sendError } from '../../helpers'
import Errors, { ApiError, ErrorType } from "../../Errors";
import logger from '../../logger';
import { ReadedJob } from '../../types';


const GetItpHistoryRouter = Router();

// 'src/routes/history/get.ts'
GetItpHistoryRouter.get('/:job_id', async (req, res) => {
  const job_id = req.params.job_id
  console.log("ITP", job_id)
  const job_details = await HistoryOrganizer.getJob(job_id); // Trust user id exists, not check on user ownership
  // raise exception on empty results ^^
  const { files } = job_details;
  HistoryOrganizer.readFiles(job_id, files).then(readedFiles => {

    res.json( {itp : readedFiles.itp_files[0][0].content , gro : readedFiles.gro })
  }).catch(e => {
    if (e === "not_found") sendError(Errors.make(ErrorType.HistoryFilesNotFound), res)
    else sendError(Errors.make(ErrorType.Server, e), res)
  })
})

export default GetItpHistoryRouter