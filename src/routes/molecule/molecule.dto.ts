import { Transform } from "class-transformer";
import { IsAlphanumeric, IsBoolean, IsDateString, IsDefined, IsIn, IsNumber, IsNumberString, IsOptional, IsString, Matches } from "class-validator";
import { AvailableForceFields, FORCE_FIELDS, castStringTrueFalseToBoolean } from "../types";


const POSITION = ['none', 'all', 'backbone'] as const;
type AvailablePositions = typeof POSITION[number]

const CTER = ["COOH-ter"] as const; 
type AvailableCter = typeof CTER[number]

const NTER = ["NH2-ter"] as const; 
type AvailableNter = typeof NTER[number]

const CYSTEIN_BRIDGE = ['none', 'auto'] as const; 
type AvailableCysteinBridge = typeof CYSTEIN_BRIDGE[number]

const BUILDER_MODE = ["classic", "elastic", "go"] as const; 
type AvailableBuilderMode = typeof BUILDER_MODE[number]

const JOB_TYPE = ["martinize"] as const; 
type AvailableJobType = typeof JOB_TYPE[number]

class ClientSettings {
    @IsDefined()
    @IsNumberString()
    user_id!: string; 
}

export class ClientSettingsMartinize extends ClientSettings {
    @IsDefined()
    @IsIn(FORCE_FIELDS)
    ff!: AvailableForceFields

    @IsDefined()
    @IsIn(POSITION)
    position!: AvailablePositions

    @IsDefined()
    @IsIn(CTER)
    cter!: AvailableCter

    @IsDefined()
    @IsIn(NTER)
    nter!: AvailableNter

    @IsDefined()
    @IsBoolean()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    sc_fix!:boolean

    @IsDefined()
    @IsIn(CYSTEIN_BRIDGE)
    cystein_bridge!: AvailableCysteinBridge

    @IsOptional()
    @IsBoolean()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    elastic?: boolean

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    ef?: number; 

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    el?: number; 

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    eu?: number; 

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    ea?: number; 

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    ep?: number; 

    @IsOptional()
    @IsNumber()
    @Transform(({ value }) => Number(value))
    em?: number; 

    @IsOptional()
    @IsBoolean()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    use_go?: boolean

    @IsDefined()
    @IsIn(BUILDER_MODE)
    builder_mode!: AvailableBuilderMode

    @IsDefined()
    @IsBoolean()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    send_mail!: boolean

    @IsOptional()
    @IsString()
    @Matches('\^[A-Z0-9a-z.\-_]*$')
    pdb_name?: string
}

export class JobToSave {
    @IsDefined()
    @IsString()
    @IsAlphanumeric()
    jobId!: string

    @IsDefined()
    @IsNumberString()
    userId!: string

    @IsDefined()
    @IsIn(JOB_TYPE)
    type!: AvailableJobType
}