import { Dree } from "dree";
import MoleculeOrganizer from "../../../MoleculeOrganizer";
import { GoTerms } from "../../../types";
import { InfosJson, SimuFile } from "../CreateMoleculeJson";

const dree = require('dree');
const fs = require('fs');


interface TopFile {
    version: string,
    infos: SimuFile
}


var molecule : InfosJson = {
    versions: [],
    name: '',
    alias: '',
    category: 'lipids',
    create_way: '',
    directory: '',
    top: [],
    map: [],
    gro: {originalname: '', path: '', size:0}

}




/**
 * Parse the molecules files names and return the informations in a dictionary
 * @param path - Diretory where the molecules files are found
 * @param dict_molecules - Dictionnary containing the data
 * @returns - Dictionary containing the molecules informations
 */
export const parser_files = function(path : string) : InfosJson[] {
    
    // Options of the scan function : only read direct children files and directories of the root 
    // and apply fileCallback only itp end gro files
    const options = {
        depth:1,
        extensions: ['itp', 'gro', 'pdb', 'top', 'map']
    };
    
    /**
     * Callback function of dree.scan when a file is given as argument
     * @param element - File
     * @param stat
     */
    const fileCallback = function (element : Dree) {
    
        // Parsing and insertion of the itp file in the dictionary
        if (element.name.match("itp$") != null) {
            let itp_info = element.name.split("_");
            let s = 0;
            if (element.size) {
                s = Number(element.size.split(' ')[0])*100;
            }
            

            // Case with protonation
            if (element.name.match("^martini_v\.+_[A-Z]?[A-Z0-9]+_\\d{2}_[zpm]\\d.itp$") != null) {

                let tmp = {"number":itp_info[3].split('.')[0], 
                    "itp": {
                        originalname: element.name, 
                        path: element.path,
                        size: s
                    },
                    "force_field":itp_info[1], 
                    "protonation":itp_info[4].split(".")[0]
                };

                
                //dict_molecules[name_molecule]["versions"].splice(tmp['number'], 0, tmp);

                molecule.versions.splice(Number(tmp.number), 0, tmp);

                
            }
            // Cases without protonation
            
            else {
                let tmp = {"number":itp_info[3].split('.')[0], 
                    "itp": {
                        originalname: element.name, 
                        path: element.path,
                        size: s
                    }, 
                    "force_field":itp_info[1]
                };

                //dict_molecules[name_molecule]["versions"].splice(tmp['number'], 0, tmp);

                molecule.versions.splice(Number(tmp.number), 0, tmp);

            }
        } 
        

        // Insertion of the gro file in the dictionary
        else if(element.name.match("gro$") != null || element.name.match("pdb$") != null) {

            let s = 0;
            if (element.size) {
                s = Number(element.size.split(' ')[0])*100;
            }

            let gro = {
                originalname: element.name, 
                path: element.path,
                size: s
            };

            molecule.gro = gro;
            
        }

        else if (element.name.match("map$") != null) {

            let s = 0;
            if (element.size) {
                s = Number(element.size.split(' ')[0])*100;
            }

            let map = {
                originalname: element.name,
                path: element.path,
                size: s
            }
            molecule.map.push(map);
        }
        
        else {

            let s = 0;
            if (element.size) {
                s = Number(element.size.split(' ')[0])*100;
            }

            let top = {
                version: element.name.split('_')[1].split('.top')[0],
                infos: {
                    originalname: element.name, 
                    path: element.path,
                    size: s
                }
            };
            molecule.top.splice(Number(top.version), 0, top);
        }
        
    };
    
    /**
     * Callback function of dree.scan when a folder is given as argument
     * @param element - Folder
     * @param stat 
     */
    const dirCallback = function (element: Dree) {
            
        // Il the directory is not the current folder (prevent it from being scanned again), 
        // scan all the children directories
        let name_molecule = element.path.split('/')[element.path.split('/').length -1];
        if (!element.relativePath.includes('.')){
            molecule = {
                versions: [],
                name: name_molecule,
                alias: name_molecule,
                category: 'lipids',
                create_way: 'hand',
                directory: element.path,
                top: [],
                map: [],
                gro: {originalname: '', path: '', size:0}
    
            }
            const tree = dree.scan(element.path, options, fileCallback, dirCallback);
            if (molecule.versions.length > 0) {
                batch.push(molecule);
            }
        }

        
        
    };
    

    var batch : InfosJson[] = [];

    // Launch the dree.scan function with the root directory
    const tree = dree.scan(path, options, fileCallback, dirCallback);

    //console.log(batch);

    return(batch);
}


// molecule load /home/achopin/Documents/database/martini-molecule-repository/martini2_lipids_test/Glycosphingolipids
