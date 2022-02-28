import { Transform } from "class-transformer";
import { IsNotEmpty, IsString, IsIn, IsBoolean, IsOptional, Matches, IsNumber } from "class-validator";
import { CpuInfo } from "os";
import { AvailableForceFields, FORCE_FIELDS, castStringTrueFalseToBoolean } from "../types";

const BOX_TYPE = ["hexagonal", "rectangular", "square", "cubic", "optimal", "keep"] as const
type AvailableBoxType = typeof BOX_TYPE[number]

const ROTATE = ["none", "random", "princ", "angle"] as const
type AvailableRotateMode = typeof ROTATE[number]

const SOLVENT = ["W", "PW"] as const
type AvailableSolvent = typeof SOLVENT[number]


export class ClientInsaneSettingsDto {
    @IsBoolean()
    @IsNotEmpty()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    lipids_added!: boolean; 

    @IsBoolean()
    @IsNotEmpty()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    molecule_added!: boolean; 

    @IsString()
    @IsOptional()
    @Matches('^([A-Z]+:[0-9]+,*)*$') //Match examples : DIPC:1 | DIPC:1, | DIPC:1,DLPC:1 and so on. Should also check if it's available lipids. Maybe not here.
    lipids?: string; 

    @IsString()
    @IsOptional()
    @Matches('^([A-Z]+:[0-9]+,*)*$') //Match examples : DIPC:1 | DIPC:1, | DIPC:1,DLPC:1 and so on. Should also check if it's available lipids. Maybe not here.
    upper_leaflet?: string; 

    @IsString()
    @IsOptional()
    @IsIn(FORCE_FIELDS)
    force_field?: AvailableForceFields;

    @IsString()
    @IsOptional()
    @Matches('^[0-9]$')
    from_id?:string

    @IsString()
    @IsNotEmpty()
    @Matches('^([0-9]+,[0-9]+,[0-9]+,*){1,3}$') // Match 3 int separated by comma between 0 and 3 times
    box!: string; 

    @IsString()
    @IsNotEmpty()
    @IsIn(BOX_TYPE)
    pbc!: AvailableBoxType; 

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    area_per_lipid!:number

    @IsNumber()
    @IsOptional()
    @Transform(({ value }) => Number(value))
    area_per_lipid_upper?:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    bead_distance!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    fudge!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    grid_spacing!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    hydrophobic_ratio!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    random_kick_size!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    shift_protein!:number

    @IsString()
    @IsNotEmpty()
    @IsIn(ROTATE)
    rotate!: AvailableRotateMode;

    @IsNumber()
    @IsOptional()
    @Transform(({ value }) => Number(value))
    rotate_angle?:number

    @IsBoolean()
    @IsOptional()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    center?: boolean; 

    @IsBoolean()
    @IsOptional()
    @Transform(({ value }) => castStringTrueFalseToBoolean(value))
    orient?: boolean; 

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    salt_concentration!:number

    @IsNumber()
    @IsNotEmpty()
    @Transform(({ value }) => Number(value))
    charge!:number

    @IsString()
    @IsNotEmpty()
    @IsIn(SOLVENT)
    solvent_type!: AvailableSolvent;
}

export class FileDto{
    @IsString()
    @IsNotEmpty()
    @Matches('^[A-Z0-9a-z._-]*$')
    originalname!: string; 


}
