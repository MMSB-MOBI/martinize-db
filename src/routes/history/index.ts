import { Router } from 'express';
import ListHistoryRouter from './list'; 
import GetHistoryRouter from './get'

const HistoryRouter = Router();

HistoryRouter.use('/list', ListHistoryRouter); 
HistoryRouter.use('/get', GetHistoryRouter); 

export default HistoryRouter; 

