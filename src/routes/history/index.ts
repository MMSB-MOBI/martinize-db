import { Router } from 'express';
import ListHistoryRouter from './list'; 
import GetHistoryRouter from './get'
import DeleteHistoryRouter from './delete'
import UpdateHistoryRouter from './update'
import GetItpHistoryRouter from './itp';
 

const HistoryRouter = Router();

HistoryRouter.use('/list', ListHistoryRouter); 
HistoryRouter.use('/get', GetHistoryRouter); 
HistoryRouter.use('/delete', DeleteHistoryRouter); 
HistoryRouter.use('/update', UpdateHistoryRouter); 
HistoryRouter.use('/itp', GetItpHistoryRouter); 
 
export default HistoryRouter; 

