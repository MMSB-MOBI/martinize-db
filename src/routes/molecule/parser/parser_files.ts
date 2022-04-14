import { Dree } from "dree";
import ItpFile from "itp-parser";
import { relative } from "path";
import { Excel } from "../../../cli/molecule_cli";
import Errors, { ErrorType } from "../../../Errors";
import logger from "../../../logger";
import { CategoryTree, GoTerms, SettingsJson } from "../../../types";
import { InfosJson, SimuFile } from "../CreateMoleculeJson";
import { SETTINGS_FILE } from "../../../constants";
import { exec } from "child_process";
import path from 'path'
import { rtrim } from "../../../helpers";

const dree = require('dree');
const fs = require('fs');

const SETTINGS : SettingsJson = JSON.parse(fs.readFileSync(SETTINGS_FILE, 'utf-8'))

const AVAILABLE_CAT_DIR = Object.values(SETTINGS.category_tree).map(catObj => catObj.dir)

const AVAILABLE_FF = SETTINGS.force_fields

export const decodeCategory = (cat: string) => {
    return SETTINGS.category_tree[cat].name
}

const inverseCatTree = (catTree: CategoryTree) : {[key: string]: string} => {
    const newObj: any = {}
    for (const cat in catTree){
        newObj[catTree[cat].dir] = cat
    }

    return newObj
}

const CORR_CAT_TREE = inverseCatTree(SETTINGS.category_tree)

interface ItpParsed {
    name: string; 
    alias : string;
    category: string; 
    forceField : string; 
    version : string; 
    references? : string; 
    cmdLine?: string; 
    comments? : string;
}


const getDirInside = (path : string) => {
    return fs.readdirSync(path).filter((dir: string) => fs.lstatSync(`${path}/${dir}`).isDirectory())
}

export function* dirRecursive (path: string) : Generator<string> {
    const subdir = getDirInside(path)
    if(subdir.length === 0){
        yield path
    }
    for (const dir of subdir){
        yield* dirRecursive(`${path}/${dir}`)
    }
}

const createItpx = (input_itp: string, output_itp : string, category : string, ff: string, version: string, erase: boolean) => {
    const regexLipidsName = new RegExp('^;;;;;; Martini lipid topology for (?<name>.+),')
    const regexLipidsName2 = new RegExp('^;;;;;; Martini lipid topology for (?<name>.+)')
    const regexCmdLine = new RegExp('^;.+Args are: (?<cmdline>.+)')

    const itp = ItpFile.readFromString(fs.readFileSync(input_itp, 'utf-8'))
    if(erase) fs.writeFileSync(output_itp, itp.toString()) //Copy the original itp

    const alias = itp.name
    const fileName = rtrim(path.basename(input_itp), ".itp")

    let lipidsRegexMatch;

    if(itp.headlines.length > 0){
        lipidsRegexMatch = itp.headlines[0].match(regexLipidsName)
        if (! lipidsRegexMatch) {
            lipidsRegexMatch = itp.headlines[0].match(regexLipidsName2)
        }
    }
    const cmdLineRegexMatch = itp.headlines.length > 1 ? itp.headlines[1].match(regexCmdLine) : undefined

    const name = lipidsRegexMatch?.groups?.name ?? fileName
    const cmdline = cmdLineRegexMatch?.groups?.cmdline ?? undefined

    //Add comments to headlines
    itp.headlines.push('; Category', `;   ${category}`, ';', '; Name', `;   ${name}`, ';', '; Alias', `;   ${alias}`, ';', '; Force field', `;   ${ff}`, ';', '; Version', `;   ${version}`, ';')
    if(cmdline){
        itp.headlines.push('; Command line', `;   ${cmdline}`, ';')
    }

    const new_itp = itp.toString();
    const itpToWrite = erase ? input_itp : output_itp
    fs.writeFileSync(itpToWrite, new_itp)

}


export const completeItpFiles = (path: string, erase = false) => {
    const exceptions : {[key:string]:string[]} = {'ffNotFound': [], 'categoryNotFound': [], 'itpNumber': [], 'other': []}
    const categories = getDirInside(path)
    for (const cat of categories){
        const catPath = `${path}/${cat}`
        if (!(AVAILABLE_CAT_DIR.includes(cat))){
            exceptions['categoryNotFound'].push(catPath)
            continue
        }
        
        const ffVersion = getDirInside(catPath)
        for (const ff of ffVersion) {
            const ffPath = `${catPath}/${ff}`
            if(!(AVAILABLE_FF.includes(ff))){
                exceptions['ffNotFound'].push(ffPath)
                continue
            }
            const molecules = getDirInside(ffPath)
            for (const mol of molecules){
                const molPath = `${ffPath}/${mol}`
                const files = fs.readdirSync(molPath)
                const itp = files.filter((f: string) => f.endsWith('.itp'))
                if (itp.length !== 1){
                    console.error('ERR : not 1 itp')
                    exceptions['itpNumber'].push(molPath)
                    continue
                }
                const itp_file = itp[0]
                const itpPath = `${molPath}/${itp_file}`
                //console.log('itp', itp, 'gro', gro)
                const outItp = erase ? `${itpPath}.ori` : `${itpPath}.new`
                createItpx(itpPath, outItp, cat, ff, '1.0', erase)
            }
        }
           
    }
    logger.warn("CAN'T CREATE ITPS FOR FILES IN FOLLOWING DIRECTORIES :")
    for (const why in exceptions) {
        logger.warn(`# ${why} :`)
        console.log(exceptions[why].join("\n"))
    }
}


/**
 * Parse the molecules files names and return the informations in a Json Object
 * @param path - Directory where the molecules files are found
 * @param type - Type of molecule inserted
 * @returns - Object containing the molecules informations
 */

/*export const _parser_files = function(path : string) : InfosJson[] { //, type : keyof typeof GoTerms[]) : InfosJson[] {
    
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
    /*const fileCallback = function (element : Dree) {
    
        // Parsing and insertion of the itp file in the Json
        if (element.name.match("itp$") != null) {
            let itp_info = element.name.split("_");
            let sizeFile = 0;
            if (element.size) {
                sizeFile = Number(element.size.split(' ')[0])*100;
            }

            let infos_itp = Object.create(null);
                let file = fs.readFileSync(element.path, {encoding:'utf8', flag:'r'});
                let file_splitten = file.split('\n')
                file_splitten.shift();
                let comments = []; 
                for(let line of file_splitten) {
                    line = line.trim(); 
                    if (line.startsWith(';')) {
                        line = line.trim();
                        comments.push(line);
                    } else {
                        break;
                    }
                };

                let comments_parsed = '';
                let citation_version = '';
                let champ = '';
                for (let line2 of comments) {
                    if (line2.match(/^; \S+/)) {
                        champ = line2;
                        infos_itp[line2] = '';
                    } else {
                        infos_itp[champ] += line2;
                    }
                }
                Object.keys(infos_itp).forEach((element: string) => {
                    element = element.trim(); 
                    if (element === '; name:') {
                        MOLECULE.name = infos_itp['; name:'].split(';\t')[1];
                    } else if (element === '; Category:') {
                        let category = infos_itp['; Category:'].split(';\t')[1].split(',');
                        for(let cat of category){
                            cat = cat.trim();
                            if (cat.localeCompare('lipids', undefined, { sensitivity: 'base' }) === 0){
                                //@ts-ignore
                                MOLECULE.category.includes('MC:0005') ? undefined : MOLECULE.category.push('MC:0005');
                            }
                            else if (cat.localeCompare('sugars', undefined, { sensitivity: 'base' }) === 0){
                                //@ts-ignore
                                MOLECULE.category.includes('MC:0001') ? undefined : MOLECULE.category.push('MC:0001');
                            }
                            else if (cat.localeCompare('proteins', undefined, { sensitivity: 'base' }) === 0){
                                //@ts-ignore
                                MOLECULE.category.includes('MC:0002') ? undefined : MOLECULE.category.push('MC:0002');
                            }  
                            else if (cat.localeCompare('polymers', undefined, { sensitivity: 'base' }) === 0){
                                //@ts-ignore
                                MOLECULE.category.includes('MC:0003') ? undefined : MOLECULE.category.push('MC:0003');
                            }
                            else if (cat.localeCompare('amino acids', undefined, { sensitivity: 'base' }) === 0){
                                //@ts-ignore
                                MOLECULE.category.includes('MC:0004') ? undefined :  MOLECULE.category.push('MC:0004');
                            }
                            else {
                                logger.warn(cat.localeCompare('lipids', 'en', { sensitivity: 'base' }));
                                logger.warn("Category " + cat + " not recognized");
                                throw new Error("Category " + cat + " not recognized");
                            }
                        }
                        
                    } else if (element === '; Reference(s):') {
                        citation_version = infos_itp['; Reference(s):'];
                    } else if (element !== '') {
                        comments_parsed += '\n' + element + '\n' + infos_itp[element];
                    }
                });

            

            // Case with protonation
            if (element.name.match("^martini_v\.+_[A-Z]?[A-Z0-9]+_\\d{2}_[zpm]\\d.itp$") != null) {

                let tmp = {"number":itp_info[3].split('.')[0], 
                    "itp": {
                        originalname: element.name, 
                        path: element.path,
                        size: sizeFile
                    },
                    "force_field":itp_info[1], 
                    "protonation":itp_info[4].split(".")[0],
                    "comments":comments_parsed,
                    "citation":citation_version,
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
                    "force_field":itp_info[1],
                    "comments":comments_parsed,
                    "citation":citation_version,
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
    /*const dirCallback = function (element: Dree) {
            

        // Initialize a molecule with its name
        let name_molecule = element.path.split('/')[element.path.split('/').length -1];
        MOLECULE = {
            versions: [],
            name: '',
            alias: name_molecule,
            //@ts-ignore
            category: [],
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

    console.log("FIRST MOL", batch[0].versions[0].itp)

    return(batch);

}*/

export const parser_files = (path: string) => { 
    const mandatoryHeadlineFields = {name : "; Name", alias : "; Alias", category : '; Category', forceField : '; Force field', version : '; Version'}
    const optionalHeadlineField = {references : "; Reference(s):", cmdLine: "; Command line", comments: "0"}
    const exceptions : {[key: string]: string[]} = {'itpNumber': [], 'groNumber': [], 'topNumber':[], 'missingField': [], 'other' : []}
    let exceptionsOccurs = false
    const batch : InfosJson[] = []
    const itpDirIterator = dirRecursive(path)
    for (const itpDir of itpDirIterator){
        console.log(itpDir)
        const files = fs.readdirSync(itpDir)
        const itp = files.filter((f: string) => f.endsWith('.itp'))
        const gro = files.filter((f : string) => f.endsWith('.gro'))
        const maps = files.filter((f : string) => f.endsWith('.map'))
        const top = files.filter((f: string) => f.endsWith('.top'))

        if (itp.length !== 1) {
            exceptionsOccurs = true
            exceptions.itpNumber.push(itpDir)
            continue
        }
        if (gro.length !== 1){
            exceptionsOccurs = true
            exceptions.groNumber.push(itpDir)
            continue
        }
        if(top.length !== 1){
            exceptionsOccurs = true
            exceptions.topNumber.push(itpDir)
            continue
        }

  

        try {
            const itpFile = `${itpDir}/${itp[0]}`
            const groFile = `${itpDir}/${gro[0]}`
            const topFile = `${itpDir}/${top[0]}`
            const itpInfos : ItpParsed = parseItp(itpFile, mandatoryHeadlineFields, optionalHeadlineField)
            
            
            const molecule: InfosJson = {
                versions : [{
                    number: itpInfos.version,
                    itp : {originalname : itp[0], path: itpFile, size : 0},
                    force_field : itpInfos.forceField,
                    comments : itpInfos.comments ?? '',
                    citation : itpInfos.references ?? ''
                }], 
                name : itpInfos.name, 
                alias : itpInfos.alias,
                category: [CORR_CAT_TREE[itpInfos.category]],
                create_way : "hand",
                directory : itpDir, 
                top : [{version : itpInfos.version, infos : {originalname : top[0], path: topFile, size: 0}}],
                map : maps.map((m:string) => ({originalname: m, path: `${itpDir}/${m}`, size: 0})),
                gro : {originalname : gro[0], path: groFile, size : 0}
            }
            batch.push(molecule)

        } catch(e) {
            //Catch the missing field error
            exceptionsOccurs = true
            exceptions['other'].push(itpDir)
            continue
        }
     
        
    }

    if(exceptionsOccurs){
        logger.warn("THESE MOLECULES CAN'T BE LOADED")
        for (const why in exceptions) {
            logger.warn(`# ${why}`)
            console.log(exceptions[why].join("\n"))
        }
    }

    return batch

}

const parseItp = (file : string, mandatoryHeadlineFields : {[key:string]: string}, optionalHeadlineFields: {[key: string]: string}) => {
    const itp = ItpFile.readFromString(fs.readFileSync(file, 'utf-8'))

    const headlines = itp.headlines

    const symbolToParse = [...Object.values(mandatoryHeadlineFields), ...Object.values(optionalHeadlineFields)]

    const parsedHeadlines = parseHeadlines(headlines, symbolToParse)

    const variables : any = {}

    for (const field in optionalHeadlineFields){
        const term = optionalHeadlineFields[field]
        variables[field] = term in parsedHeadlines ? flatItpField(parsedHeadlines[term]) : ''
    }

    for(const field in mandatoryHeadlineFields) {
        const term = mandatoryHeadlineFields[field]
        if (!(term in parsedHeadlines)){
            throw new Error(`No ${field} field in itp`)
        }
        variables[field] = flatItpField(parsedHeadlines[term])
    }

    return variables as ItpParsed
}

const parseHeadlines = (headlines : string[], keys: string[]) => {
    let parsed : {[key: string]: string[]} = {'0':[]}
    let currentParsed = '0'
    for(const line of headlines){
        let delThisEmptyLine = false; 
        if(line === ";" && currentParsed !== '0'){
            delThisEmptyLine = true
            currentParsed = '0'
        }
        if(keys.includes(line)){
            currentParsed = line
            
            if(line in parsed){
                console.error(`${line} already parsed, should not happen`)
                continue
            }
            parsed[currentParsed] = []
           
        } else {
            if(!delThisEmptyLine) parsed[currentParsed].push(line)
        }
        
    }
    return parsed
}

const flatItpField = (field: string[]) => {
    const flatArray = []
    for (const line of field){
        const trimmedLine = line[0] === ";" ? line.substring(1).trim() : line.trim()
        flatArray.push(trimmedLine)
    }
    return flatArray.join('\n')
}