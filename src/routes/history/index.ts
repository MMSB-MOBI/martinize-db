import { Router } from 'express';
import ListHistoryRouter from './list'; 

const HistoryRouter = Router();

HistoryRouter.use('/list', ListHistoryRouter); 

export default HistoryRouter; 

