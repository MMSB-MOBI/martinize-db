import { Router } from 'express';
import ListHistoryRouter from './list'; 
import GetHistoryRouter from './get'
import DeleteHistoryRouter from './delete'
import UpdateHistoryRouter from './update'

const HistoryRouter = Router();

HistoryRouter.use('/list', ListHistoryRouter); 
HistoryRouter.use('/get', GetHistoryRouter); 
HistoryRouter.use('/delete', DeleteHistoryRouter); 
HistoryRouter.use('/update', UpdateHistoryRouter); 

export default HistoryRouter; 

