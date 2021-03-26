const { ItpFile, TopFile } = require('itp-parser');
const fs = require('fs');
import * as path from 'path';
const dree = require('dree');
/*

export const create_top = async function(itp_dir: string, outdir: string){

    const top = TopFile.readFromString('');
    const itp = itp_dir;

    //top.setField(TopFile.HEADLINE_KEY, []);
    
    top.headlines.push(`#include "martini.itp"`);
    top.headlines.push(`#include "${path.basename(itp)}"`);


    top.setField('system', ['This is an auto generated system']);
    top.appendFieldLine('molecules', ';moleculetype\tcount');

    let itp_file = ItpFile.readFromString(fs.readFileSync(itp, 'utf-8'));
    top.appendFieldLine('molecules', `${itp_file.name}\t1`);

    //console.log(itp_file.name);


    const content = top.toString();

    let version = path.basename(itp).split('_').slice(-1).join().split('.').slice(0, -1).join();

    fs.writeFileSync(`${outdir}/${itp_file.name}_${version}.top`, content);

}


export const create_top_in_dir = function(indir:string){

    const options = {
        depth:10,
        extensions: ['itp', 'top'],
        exclude: /0\d\.top/
    };

    const tree = dree.scan(indir, options, (element)=> {
        
        if(element.path.match('.itp$')){
            create_top(element.path, element.path.split('/').slice(0, -1).join('/'));
        }
        

    }, (element) => {
            element.children.forEach(child => {
                if(child.name.match('top')){
                    fs.rename(child.path, child.path.split('/').slice(0, -1).join('/')+'/'+child.path.split('/').slice(-2, -1).join('/')+'_01.top', function(err) {
                        if ( err ) console.log('ERROR: ' + err);
                    });

                }
            });
            
    });
}

*/