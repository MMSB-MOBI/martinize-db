import { CONECT_PDB_PATH, CREATE_GO_PATH, CREATE_MAP_PATH, INSANE_PATH, MARTINIZE_PATH, MARTINIZE_VENV, RUN_POLYPLY_PATH, INSANE_VENV, POLYPLY_VENV, JOB_MANAGER_SETTINGS, MINIMIZEPDB, SLURM_PROFILES, INIT_POLYPLY_PATH, CREATE_MAP_RCSU_PATH } from '../constants';
import { ArrayValues, generateSnowflake } from '../helpers';

const SupportedScripts = ['insane', 'conect', 'convert', 'go_virt', 'ccmap', 'martinize', 'polyply', 'get_residue_avaible', 'map_rcsu'] as const;
export type SupportedScript = ArrayValues<typeof SupportedScripts>;
import jmClient from 'ms-jobmanager'
import logger from '../logger';
import path from 'path';
import { Readable } from 'stream';
import { ErrorType } from '../Errors';


export interface JobInputs {
    exportVar: { [key: string]: string },
    inputs: { [key: string]: any },
    modules? : string[],
    script? : string,
    jobProfile? : string,
    sysSettingsKey? : string, 
  };


type RunMode = 'server' | 'local';
const DEFAULT_RUN_MODE = 'local'

const SCRIPTS: { [scriptName in SupportedScript]: string } = {
    'conect': CONECT_PDB_PATH,
    'convert': MINIMIZEPDB,
    'go_virt': CREATE_GO_PATH,
    'ccmap': CREATE_MAP_PATH,
    'insane': INSANE_PATH,
    'martinize': MARTINIZE_PATH,
    'polyply': RUN_POLYPLY_PATH,
    'get_residue_avaible': INIT_POLYPLY_PATH,
    'map_rcsu': CREATE_MAP_RCSU_PATH
};

const SERVER_MODULES: { [scriptName in SupportedScript]: string[] } = {
    'conect': ['gromacs'],
    'convert': ['gromacs'],
    'go_virt': ['mad-utils'],
    'ccmap': ['mad-utils'],
    'insane': ['insane'],
    'martinize': ['martinize2'],
    'polyply': ['polyply'],
    'get_residue_avaible': ['polyply'],
    'map_rcsu': ['rcsu']
}

const LOCAL_CONFIG: { [scriptName in SupportedScript]: { venv?: string, modules?: string[] } } = {
    'conect': { modules: ['gromacs/2020.7'] },
    'convert': { modules: ['gromacs/2020.7'] },
    'go_virt': { venv: MARTINIZE_VENV },
    'ccmap': { venv: MARTINIZE_VENV },
    'insane': { venv: INSANE_VENV },
    'martinize': { venv: MARTINIZE_VENV },
    'polyply': { venv: POLYPLY_VENV },
    'get_residue_avaible': { venv: POLYPLY_VENV },
    'map_rcsu': {}
}

export default new class JMSurcouche {
    public mode: RunMode = DEFAULT_RUN_MODE;
    public id = "";
    async run(what_to_launch: SupportedScript, args: JobInputs) {
        if (this.id === "") this.id = generateSnowflake()
        args.script = SCRIPTS[what_to_launch]

        if (this.mode === 'local') {
            const venv = LOCAL_CONFIG[what_to_launch].venv
            const modules = LOCAL_CONFIG[what_to_launch].modules
            if (venv) args.exportVar = { ...args.exportVar, venv }
            if (modules) args.modules = modules

        } else {
            args.modules = SERVER_MODULES[what_to_launch]
            args.jobProfile = SLURM_PROFILES.JOB_PROFILE,
            args.sysSettingsKey = SLURM_PROFILES.SYS_SETTINGS
        } 

        logger.debug('Launch job : ' + what_to_launch + ' with mode ' + this.mode)
        
        jmClient.start(JOB_MANAGER_SETTINGS.address, JOB_MANAGER_SETTINGS.port)

        return await jmClient.pushFS(args)

    }
}

export const pathsToInputs = (paths : string[]) : {[fileName : string ] : string} => {
    const inputs: {[fileName : string ] : string} = {}
    for (const p of paths) {
        inputs[path.basename(p)] = p
    }
    return inputs
} 

export const str_to_stream = (str: string) => {
    const ma_stream: Readable = new Readable();
    ma_stream.push(str);
    ma_stream.push(null)
    return ma_stream
}

