import Errors, {ErrorType} from '../../Errors'; 

//REPLACE THIS BY WEBSOCKET GATEWAY

export const FORCE_FIELDS = [
    "elnedyn",
    "elnedyn22",
    "elnedyn22p",
    "martini22",
    "martini22p",
    "martini3001"
] as const
export type AvailableForceFields = typeof FORCE_FIELDS[number]

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

export interface ClientSettings {
    ff : string, 
    position : string; 
    cter : string
    nter : string 
    sc_fix : string; 
    cystein_bridge : string; 
    elastic? : string; 
    ef? : string; 
    el? : string; 
    eu? : string; 
    ea? : string; 
    ep? : string; 
    em? : string; 
    use_go? : string; 
    builder_mode : string; 
    send_mail : string; 
    user_id? : string; 
    pdb_name?: string; 
}

export interface ClientSettingsValidate {
  ff : AvailableForceFields,
  position : AvailablePositions,
  cter : AvailableCter, 
  nter : AvailableNter, 
  side_chain_fix? : boolean, 
  cystein_bridge : AvailableCysteinBridge,
  elastic? : boolean
  ef? : number; 
  el? : number; 
  eu? : number; 
  ea? : number; 
  ep? : number; 
  em? : number; 
  use_go ? : boolean; 
  builder_mode : AvailableBuilderMode; 
  send_mail : Boolean;
  user_id? : number;  
  pdb_name?: string ; //string that ends with .pdb
}


export function validateClientSettings(settings : ClientSettings) : ClientSettingsValidate{

  let validated = {} as ClientSettingsValidate; 

  if(settings.ff){
    const ff = settings.ff as AvailableForceFields
    if(FORCE_FIELDS.includes(ff)) validated.ff = ff
    else return Errors.throw(ErrorType.Format)
  } 
  else{
    return Errors.throw(ErrorType.MissingParameters)
  }

  if(settings.position) {
    const position = settings.position as AvailablePositions
    if(POSITION.includes(position)) validated.position = position
    else throw Errors.throw(ErrorType.Format)
  } else return Errors.throw(ErrorType.MissingParameters)

  if(settings.cter) {
    const cter = settings.cter as AvailableCter
    if(CTER.includes(cter)) validated.cter = cter
    else throw Errors.throw(ErrorType.Format)
  } else return Errors.throw(ErrorType.MissingParameters)

  if(settings.nter) {
    const nter = settings.nter as AvailableNter
    if(NTER.includes(nter)) validated.nter = nter
    else return Errors.throw(ErrorType.Format)
  } else return Errors.throw(ErrorType.MissingParameters)

  if(checkBooleanString(settings.sc_fix)) validated.side_chain_fix = true; 

  if(settings.cystein_bridge){
    const cystein_bridge = settings.cystein_bridge as AvailableCysteinBridge
    if(CYSTEIN_BRIDGE.includes(cystein_bridge)) validated.cystein_bridge = cystein_bridge
    else return Errors.throw(ErrorType.Format)
  } else return Errors.throw(ErrorType.MissingParameters)

  if(settings.elastic === "true") validated.elastic = true; 

  if(settings.ef) validated.ef = numberOrError(settings.ef)
  if(settings.el) validated.el = numberOrError(settings.el)
  if(settings.eu) validated.eu = numberOrError(settings.eu)
  if(settings.ea) validated.ef = numberOrError(settings.ea)
  if(settings.ep) validated.ef = numberOrError(settings.ep)

  if(settings.use_go === "true") validated.use_go = true

  if(settings.builder_mode){
    const builder_mode = settings.builder_mode as AvailableBuilderMode
    if(BUILDER_MODE.includes(builder_mode)) validated.builder_mode = builder_mode
    else return Errors.throw(ErrorType.Format)
  } else return Errors.throw(ErrorType.MissingParameters)

  if(checkBooleanString(settings.send_mail)) validated.send_mail = true
  else validated.send_mail = false

  if(settings.user_id) validated.user_id = numberOrError(settings.user_id)

  if(settings.pdb_name){
    if(settings.pdb_name.endsWith(".pdb")) validated.pdb_name = settings.pdb_name
    else return Errors.throw(ErrorType.Format)
  }

  return validated

}

function checkBooleanString(val?: string) {
  if (val){
    if(isBooleanString(val)){
      return val === "true" ? true : false
    } else Errors.throw(ErrorType.Format)
    if(val === "true") return true
    else Errors.throw(ErrorType.Format)
  }
  else Errors.throw(ErrorType.MissingParameters)
}

function isBooleanString(val : string): boolean {
  return val === "true" || val === "false"
}

function numberOrError(num: string) {
  const pos = Number(num);
  if (isNaN(pos) || pos < 0) {
    return Errors.throw(ErrorType.Format);
  }
  return pos;
}

