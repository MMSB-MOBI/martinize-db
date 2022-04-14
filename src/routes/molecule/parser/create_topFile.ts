import {ItpFile, TopFile} from 'itp-parser'
const fs = require('fs');
import * as path from 'path';
import Errors, { ErrorType } from '../../../Errors';
import logger from '../../../logger';
import { dirRecursive } from './parser_files';
import { rtrim } from '../../../helpers';


/**
 * Create a top file for an itp using the TopFIle interface
 * @param itp_dir - path to the itp file
 * @param outdir - path to the moelcule files where the top file will be saved
 * @returns 
 */
export const create_top_file = function(itp: string, outdir: string){
    const top = TopFile.readFromString('');
    
    top.headlines.push(`#include "martini.itp"`);
    top.headlines.push(`#include "${path.basename(itp)}"`);

    top.setField('system', ['This is an auto generated system']);
    top.appendFieldLine('molecules', ';moleculetype\tcount');

    const itp_file = ItpFile.readFromString(fs.readFileSync(itp, 'utf-8'));
    const mol_name = itp_file.name
    const itp_file_name = rtrim(path.basename(itp), ".itp")
    

    if (mol_name != '') {
        top.appendFieldLine('molecules', `${mol_name}\t1`);
    }
    else {
        return Errors.throw(ErrorType.IncorrectItpName);
    }
    
    const content = top.toString();
    fs.writeFileSync(`${outdir}/${itp_file_name}.top`, content);

}


/**
 * Create a top file for all the itp for the molecules in a directory when it does not exist
 * @param indir - Directory containing the molecules directories
 */
export const create_top_in_dir = (indir:string) => {
    const exceptions : {[key:string]: string[]} = {'itp number' : [], 'other': []}
    let exceptionsOccurs = false
    for (const itpDir of dirRecursive(indir)){
        try {
            const files = fs.readdirSync(itpDir)
            const itp = files.filter((f: string) => f.endsWith('.itp'))
            const top = files.filter((f: string) => f.endsWith('.top'))

            if (top.length > 1) {
                exceptions['top number'].push(itpDir)
                exceptionsOccurs = true
                continue
            } 
            else if (top.length === 1) {
                const topPath = `${itpDir}/${top[0]}`
                logger.debug(`mv already existing top ${top[0]}`)
                fs.renameSync(topPath, `${topPath}.ori`)
            }
            if (itp.length !== 1) {
                exceptions['itp number'].push(itpDir)
                exceptionsOccurs = true
                continue
            }
            const itpFile = `${itpDir}/${itp[0]}`
            create_top_file(itpFile, itpDir)
        } catch(e) {
            logger.warn(e)
            exceptions['other'].push(itpDir)
            exceptionsOccurs = true
        }
           

    }

    if(exceptionsOccurs){
        logger.warn("TOP CAN'T BE CREATED FOR THE FOLLOWING MOLECULES")
        for (const why in exceptions){
            logger.warn(`# ${why}:`)
            console.log(exceptions[why].join("\n"))
        }
    }


}

