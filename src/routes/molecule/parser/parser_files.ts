import { Dree } from "dree";
import ItpFile from "itp-parser";
import { relative } from "path";
import { Excel } from "../../../cli/molecule_cli";
import Errors, { ErrorType } from "../../../Errors";
import logger from "../../../logger";
import { GoTerms } from "../../../types";
import { InfosJson, SimuFile } from "../CreateMoleculeJson";

const dree = require('dree');
const fs = require('fs');


export var MOLECULE : InfosJson = {
    versions: [],
    name: '',
    alias: '',
    category: 0,
    create_way: '',
    directory: '',
    comments: '',
    citation: '',
    top: [],
    map: [],
    gro: {originalname: '', path: '', size:0}
}




/**
 * Parse the molecules files names and return the informations in a Json Object
 * @param path - Directory where the molecules files are found
 * @param type - Type of molecule inserted
 * @returns - Object containing the molecules informations
 */
export const parser_files = function(path : string) : InfosJson[] { //, type : keyof typeof GoTerms[]) : InfosJson[] {
    
    // Options of the scan function : only read direct children files and directories of the root 
    // and apply fileCallback only itp end gro files
    const options = {
        depth:2,
        extensions: ['itp', 'gro', 'pdb', 'top', 'map']
    };
    
    /**
     * Callback function of dree.scan when a file is given as argument
     * @param element - File (dree format)
     */
    const fileCallback = function (element : Dree) {
    
        // Parsing and insertion of the itp file in the Json
        if (element.name.match("itp$") != null) {
            let itp_info = element.name.split("_");
            let sizeFile = 0;
            if (element.size) {
                sizeFile = Number(element.size.split(' ')[0])*100;
            }

            //let molecule_itp: ItpFile;
            let infos_itp = Object.create(null);
            console.log("test");
            //ItpFile.read(element.relativePath).then((file) => {
                let file = fs.readFileSync(element.path, {encoding:'utf8', flag:'r'});
                let file_splitten = file.split('\n')
                file_splitten.shift();
                //molecule_itp = file;
                let comments = []; //molecule_itp!.getField('headlines');
                for(let line of file_splitten) {
                    if (line.startsWith(';')) {
                        comments.push(line);
                    } else {
                        break;
                    }
                };
                console.log('comments')

                let champ = '';
                for (let line2 of comments) {
                    if (line2.startsWith('; ')) {
                        champ = line2;
                        infos_itp[line2] = '';
                    } else if (line2.startsWith(';\t')) {
                        infos_itp[champ] += line2;
                    }
                }
                console.log(infos_itp);
                infos_itp.keys.forEach((element: string) => {
                    if (element === '; Name') {
                        MOLECULE.name = infos_itp['; Name:'].split(';\t')[1];
                    } else if (element === '; Categories:') {
                        // TODO maybe check if it is a goname
                        //@ts-ignore
                        MOLECULE.category = infos_itp['; Categories:'].split(';\t')[1].split(',');
                    } else if (element === '; Reference(s):') {
                        MOLECULE.citation = infos_itp['; Reference(s):'].split(';\t')[1];
                    } else if (element !== '') {
                        MOLECULE.comments += element + '\n' + infos_itp[element];
                    }
                });
                console.log(MOLECULE);
            

            

            // Case with protonation
            if (element.name.match("^martini_v\.+_[A-Z]?[A-Z0-9]+_\\d{2}_[zpm]\\d.itp$") != null) {

                let tmp = {"number":itp_info[3].split('.')[0], 
                    "itp": {
                        originalname: element.name, 
                        path: element.path,
                        size: sizeFile
                    },
                    "force_field":itp_info[1], 
                    "protonation":itp_info[4].split(".")[0]
                };


                MOLECULE.versions.splice(Number(tmp.number), 0, tmp);

                
            }
            // Cases without protonation
            
            else  if (element.name.match("^martini_v\.+_[A-Z]?[A-Z0-9]+_\\d{2}.itp$") != null) {
                let tmp = {"number":itp_info[3].split('.')[0], 
                    "itp": {
                        originalname: element.name, 
                        path: element.path,
                        size: sizeFile
                    }, 
                    "force_field":itp_info[1]
                };


                MOLECULE.versions.splice(Number(tmp.number), 0, tmp);

            }

            // Error syntax
            else {
                return Errors.throw(ErrorType.IncorrectItpName);
            }
        } 
        

        // Insertion of the gro file in the Json
        else if(element.name.match("gro$") != null || element.name.match("pdb$") != null) {

            let sizeFile = 0;
            if (element.size) {
                sizeFile = Number(element.size.split(' ')[0])*100;
            }

            let gro = {
                originalname: element.name, 
                path: element.path,
                size: sizeFile
            };

            MOLECULE.gro = gro;
            
        }

        // Insertion of the map file if exists in the Json
        else if (element.name.match("map$") != null) {

            let sizeFile = 0;
            if (element.size) {
                sizeFile = Number(element.size.split(' ')[0])*100;
            }

            let map = {
                originalname: element.name,
                path: element.path,
                size: sizeFile
            }
            MOLECULE.map.push(map);
        }
        
        // Insertion of the top file in the Json
        else {

            let sizeFile = 0;
            if (element.size) {
                sizeFile = Number(element.size.split(' ')[0])*100;
            }

            let top = {
                version: element.name.split('_')[1].split('.top')[0],
                infos: {
                    originalname: element.name, 
                    path: element.path,
                    size: sizeFile
                }
            };
            MOLECULE.top.splice(Number(top.version), 0, top);
        }
        
    };
    

    /**
     * Callback function of dree.scan when a folder is given as argument
     * @param element - Folder (dree format)
     */
    const dirCallback = function (element: Dree) {
            

        // Initialize a molecule with its name
        let name_molecule = element.path.split('/')[element.path.split('/').length -1];
        MOLECULE = {
            versions: [],
            name: name_molecule,
            alias: name_molecule,
            category: 0,
            create_way: 'hand',
            comments: '',
            citation: '',
            directory: element.path,
            top: [],
            map: [],
            gro: {originalname: '', path: '', size:0}
        }

        // If the directory is not the current folder (prevent it from being scanned again), 
        // scan all the children directories
        if (!element.relativePath.includes('.')) {
            if (element.children){

                // If the directory contains molecule files
                if(element.children[0].type == 'file') {
                    const tree = dree.scan(element.path, options, fileCallback);
    
                }

                // If the moelcule files are in a subdirectory
                else {
                    const tree = dree.scan(element.path, options, ()=>{}, dirCallback);
                }
    
                //TODO gestion des erreurs
                // check if no molecule file is missing after parsing the directory
                if (element.children[0].type == 'file') {
                    if (MOLECULE.versions.length == 0 || MOLECULE.gro.originalname === ''){
                        let err = Errors.make(ErrorType.MissingFiles);
                        logger.warn(name_molecule + ' : ' + err.data.message);
                        Excel.text += name_molecule+ ',X,,,\n';
                    }
                    else if (MOLECULE.top.length == 0) {
                        return Errors.throw(ErrorType.MissingTopFiles);
                    }
                    else {
                        batch.push(MOLECULE);
                    }
                }
            }
        }
        
    };
    

    var batch : InfosJson[] = [];

    // Launch the dree.scan function with the root directory
    const tree = dree.scan(path, options, ()=>{}, dirCallback);

    return(batch);

}