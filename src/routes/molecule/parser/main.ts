/*const neodoc = require('neodoc');

const args = neodoc.run(`
usage: main.js DIRECTORY

options:
    -h --help  Show the help message

Before running the script make sure that the name of your itp files follow this pattern : 
    martini_(version of martini parameters in v+number format)_(name of the molecule)_(version of the molecule).itp
            OR
    martini_(version of martini parameters in v+number format)_(name of the molecule)_(version of the molecule)_(protonation state [p+number if positive, m+number if negative, z0 for zero]).itp
Make sure that no other "_" are present in the filename.
`);
*/

// 
import * as pf from "./parser_files";
import * as ct from './create_topFile';
import { CreateMoleculeFromJson } from "../CreateMoleculeJson";

/*
console.log("\nWriting the Top files");
ct.create_top_in_dir(args["DIRECTORY"]);
console.log("Top files written in the top_files directory\n");

console.log("Start of parsing");
let dict_molecules = pf.parser_files(args["DIRECTORY"]);
console.log("End of parsing\n");

*/
//console.log("Writing the Json files");
let batch = pf.parser_files(process.argv[2]);
//console.log("Json files written in the json_files directory\n");

CreateMoleculeFromJson(batch);
console.log('ok');