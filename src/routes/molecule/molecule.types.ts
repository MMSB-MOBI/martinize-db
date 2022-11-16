import { ClientSettingsMartinize } from "./molecule.dto";

type AvailableJobType = "martinize" | "polyply"

interface JobToSave {
    jobId: string
    userId: string
    type: AvailableJobType
    date: string //to type best
}

interface MartinizeFilesPath {
    all_atom : string
    coarse_grained: string
    itp_files : string[][]
    top_file: string
    warnings: string;
    
}


export interface MartinizeJobToSave extends JobToSave {
    name?: string; //pdb name
    settings: ClientSettingsMartinize
    files: MartinizeFilesPath
    radius : {[bead: string]: number}
}

export interface PolyplyJobToSave extends JobToSave {
    name?: string; //pdb name
    settings: ClientSettingsMartinize
     
     
}