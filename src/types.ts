import { JobBase } from './Entities/entities' 

export const GoTerms = {
  "MC:0001": "<Go Term Regular Name>",
  "MC:0005": "",
  "MC:0002":"",
  "MC:0003":"",
  "MC:0004":"",
};

export type UserRole = "admin" | "curator";

export interface JSONWebTokenPartial {
  /** Issued at */
  iat: string;
  /** Expiration (timestamp) */
  exp: string;
  /** Issuer */
  iss: string;
  /** ID */
  jti: string;
}

export interface TokenPayload {
  user_id: string, 
  created_at: string;
}

export type JSONWebToken = JSONWebTokenPartial & TokenPayload;

export interface SettingsJson {
  force_fields: string[];
  force_fields_info: ForceFielsdInfo
  create_way: { [wayId: string]: string };
  category_tree: CategoryTree;
}

interface ForceFielsdInfo{
  [ff_name: string]: {
    polarizable: boolean
  };
}

export interface CategoryTree {
  [go_id: string]: {
    children: CategoryTree,
    name: string;
    dir : string; 
  };
}

export interface JobFilesNames { // just job files names or job files content
  all_atom: string
  coarse_grained: string
  itp_files: string[][]
  top_file: string
  warnings: string
}

export interface JobReadedFiles { // just job files names or job files content
  all_atom: ReadedFile
  coarse_grained: ReadedFile
  itp_files: ReadedFile[][]
  top_file: ReadedFile
  warnings: ReadedFile
}

export interface ReadedFile {
  content: string; 
  type : string; 
  name : string; 
}

export interface JobSettings {
    builder_mode : "classic" | "go" | "elastic"; 
    ff : AvailableForceFields; 
    advanced : boolean; 
    commandline : string; 
    cter : "COOH-ter"; 
    nter : "NH2-ter";
    sc_fix : boolean; 
    position : "backbone" | "all" | "none"
    cystein_bridge : "none" | "auto"
    elastic? : boolean; 
    use_go? : boolean; 
    ea? : number; 
    ef? : number; 
    el? : number; 
    em? : number; 
    ep? : number; 
    eu? : number; 
}

export interface ReadedJob extends JobBase {
  files : JobReadedFiles
}

type AvailableForceFields = "martini3001" | "elnedyn22" | "elnedyn22p" | "elnedyn" | "martini22" | "martini22p"
