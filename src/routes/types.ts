export const FORCE_FIELDS = [
    "elnedyn",
    "elnedyn22",
    "elnedyn22p",
    "martini22",
    "martini22p",
    "martini3001"
] as const
export type AvailableForceFields = typeof FORCE_FIELDS[number]

export function castStringTrueFalseToBoolean(val: string) : boolean | undefined {
    //Don't work with 0 and 1
    if(val.toLowerCase() === "true") return true
    if(val.toLowerCase() === "false") return false
}
