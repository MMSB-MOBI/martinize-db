import { Router } from 'express';
import ListHistoryRouter from './list'; 
import GetHistoryRouter from './get'
import DeleteHistoryRouter from './delete'

const HistoryRouter = Router();

HistoryRouter.use('/list', ListHistoryRouter); 
HistoryRouter.use('/get', GetHistoryRouter); 
HistoryRouter.use('/delete', DeleteHistoryRouter); 

export default HistoryRouter; 

