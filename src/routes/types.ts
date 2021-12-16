export const FORCE_FIELDS = [
    "elnedyn",
    "elnedyn22",
    "elnedyn22p",
    "martini22",
    "martini22p",
    "martini3001"
] as const
export type AvailableForceFields = typeof FORCE_FIELDS[number]