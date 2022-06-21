import ItpFile from "itp-parser";
import logger from "../../../logger";
import { CategoryTree, SettingsJson } from "../../../types";
import { InfosJson, VersionItp } from "../CreateMoleculeJson";
import { SETTINGS_FILE } from "../../../constants";
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


export const parser_files = (path: string) => { 
    const mandatoryHeadlineFields = {name : "; Name", alias : "; Alias", category : '; Category', forceField : '; Force field', version : '; Version'}
    const optionalHeadlineField = {references : "; Reference(s):", cmdLine: "; Command line", comments: "0"}
    const exceptions : {[key: string]: string[]} = {'itpNumber': [], 'groNumber': [], 'topNumber':[], 'missingField': [], 'other' : []}
    let exceptionsOccurs = false
    const batch : InfosJson[] = []
    const itpDirIterator = dirRecursive(path)
    for (const itpDir of itpDirIterator){
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
            

            //Check if already in batch
            
            const alreadyInBatch = batch.find(mol => mol.alias === itpInfos.alias)
            const version : VersionItp = {
                directory : itpDir, 
                number: itpInfos.version,
                itp : {originalname : itp[0], path: itpFile, size : 0},
                force_field : itpInfos.forceField,
                comments : itpInfos.comments ?? '',
                citation : itpInfos.references ?? '',
                command_line: itpInfos.cmdLine ?? '',
                gro : {originalname : gro[0], path: groFile, size : 0},
                map : maps.map((m:string) => ({originalname: m, path: `${itpDir}/${m}`, size: 0})),
                top : {originalname : top[0], path: topFile, size : 0}
            }
            if(alreadyInBatch){
                logger.debug(`${itpInfos.alias} is already loaded, add a new version`)
                alreadyInBatch.versions.push(version)
            }
            else {
                const molecule: InfosJson = {
                    versions : [version], 
                    name : itpInfos.name, 
                    alias : itpInfos.alias,
                    category: [CORR_CAT_TREE[itpInfos.category]],
                    create_way : "hand",
                    directory : itpDir, 
                
                }
                batch.push(molecule)
            }
            

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
            logger.warn(`# ${why} ${exceptions[why].length}`)
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