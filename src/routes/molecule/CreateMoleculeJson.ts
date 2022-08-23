const fs = require('fs');
const dree = require('dree');

import { MoleculeChecker } from './MoleculeChecker';
import { Molecule, MoleculeVersion } from '../../Entities/entities';
import { Database } from '../../Entities/CouchHelper';
import Errors, { ErrorType } from '../../Errors';
import logger from '../../logger';
import { CONNECTED_USER_CLI } from '../../cli/user_cli';





/**
 * One itp file version written in the Json file
 */
export interface VersionItp extends MinVersionItp {
  directory: string; 
  itp: SimuFile,
  top: SimuFile,
  map: SimuFile[],
  gro : SimuFile
  protonation? : string,
  comments: string,
  citation: string,
  command_line: string
  inserted?: boolean;
};

export interface MinVersionItp {
  number: string,
  force_field: string
}

/**
 * All the informations on the molecule stored in the Json file used to simulate a request
 */
export interface InfosJson {
  versions: VersionItp[],
  /* Name of the molecule */
  name: string,
  alias: string,
  category: string[],
  create_way: string,
  directory: string,
};

/**
 * Simulate the files informations contained in the request sent to the server
 */
export interface SimuFile {
  originalname: string;
  path: string;
  size: number
}

/**
 * Simulate the informations contained in a request sent to the server to add a molecule to the database
 */
export interface SimuRequest{
  full_user: {
    id: string,
    role: string
  };
  body: {
    name: string,
    alias: string,
    smiles: string,
    version: string,
    category: string[],
    command_line: string,
    comments: string,
    create_way: string,
    force_field: string,
    validation: string,
    citation: string,
    parent: string | null
    tree_id? : string
  };
  files: {
    itp: SimuFile[],
    pdb: SimuFile[],
    top: SimuFile[] | [],
    map: SimuFile[] | []
  };
}

interface InsertionRecap {
  inserted : {[molecule_alias : string] : MoleculeInsert }
  not_inserted : {[reason: string]: {[molecule_alias : string] : MoleculeInsert }}
}

interface MoleculeInsert {
  name: string
  versions : VersionItp[]
  dir: string; 
}

interface InsertedVersion {
  inserted : VersionItp[]
  not_inserted : {reason:string, version : VersionItp}[]
  parent?: string; 
}


/**
 * Read a bacth of molecule objects contained in a Json object and send their infos to the database to insert them sequentially
 * @param batch - a list of molecules informations
 */
export const CreateMoleculeFromJson = async (batch : InfosJson[]) : Promise< InsertionRecap > => {
  const recap : InsertionRecap = {'inserted': {}, 'not_inserted' : {'other': {}}}
  for (const info of batch) {
    logger.info(`Try to insert ${info.name}`)
    try {
      const insertedInfo = await CreateMoleculeFromJsonAux(info)
      if(insertedInfo.inserted.length > 0) recap.inserted[info.alias] = {name: info.name, versions : insertedInfo.inserted, dir : info.directory}
      if(insertedInfo.not_inserted.length > 0) {
        for(const notInserted of insertedInfo.not_inserted) {
          if(!(notInserted.reason in recap.not_inserted)){
            recap.not_inserted[notInserted.reason] = {}
          }
          if(!(info.alias in recap.not_inserted[notInserted.reason])) recap.not_inserted[notInserted.reason][info.alias] = {name: info.name, versions : [], dir : info.directory}
          recap.not_inserted[notInserted.reason][info.alias].versions.push(notInserted.version)
          
        }
        
      }
      

    } catch(e) {
      logger.warn("insertion failed")
      if(e.data && e.data.message){
        if(!(e.data.message in recap.not_inserted)) recap.not_inserted[e.data.message] = {}
        recap.not_inserted[e.data.message][info.alias] = {name: info.name, versions : info.versions, dir: info.directory}
      }
      else {
        console.error(e)
        recap.not_inserted['other'][info.alias] = {name: info.name, versions : info.versions, dir: info.directory}
      }
      
    }
  }
  return new Promise((res, rej) => res(recap))
}
/**
 * Auxiliary function reading the infos from one molecule of the batch and inserting them in the databse
 * @param infos - The informations on the moelcule
 * @returns a Promise with resolve true :)
 */
const CreateMoleculeFromJsonAux = async (infos : InfosJson): Promise<InsertedVersion> => {
  return new Promise(async (resolve, reject) => {
    let inserted : VersionItp[] = []
    const not_inserted : {reason: string, version: VersionItp}[] = []
  let tree_id; 
  //Get versions already in database based on alias 
  const alreadyInDb = await Database.molecule.getVersionsFromAlias(infos.alias)
  if(alreadyInDb.length > 0) tree_id = alreadyInDb[0].tree_id
  //const tree = new MoleculeVersionTree(alreadyInDb); 
  const versionsToInsert = infos.versions.filter(v => {
    const inDb = alreadyInDb.find(mol => mol.force_field === v.force_field && mol.version === v.number)
    if (inDb) {
      logger.warn(`${infos.alias} ${v.force_field} ${v.number} already in database`)
      not_inserted.push({reason: "already in database", version : v})
      return false
    }
    return true
  })

  if(versionsToInsert.length === 0) {
    resolve({inserted, not_inserted})
  }

  const versionsByFf = separateVersionsByForcefield(versionsToInsert)
  
  for(const ff in versionsByFf){
    let to_insert = versionsByFf[ff]
    const versions_number = versionsByFf[ff].map(v => v.number)
    const duplicates = versions_number.filter((item,index) => versions_number.indexOf(item) != index)
    if(duplicates.length > 0){
      for (const v of duplicates){
        logger.warn(`${infos.alias} v${v} is duplicated`)
        to_insert = [versionsByFf[ff].find(ver => ver.number === v)!]
        const allOthers = versionsByFf[ff].filter(ver => ver.number === v).slice(1, versionsByFf[ff].length)
        for (const mol of allOthers){
          not_inserted.push({reason : "duplicate version", version: mol})
        }
      }
      
    }
    const {nodes, notTreated} = versionTree(to_insert)
    if(notTreated.length > 0){
      for (const v of notTreated){
        logger.warn(`${infos.alias} v${v} has no available parent version. Not inserted.`)
        not_inserted.push({reason : 'no previous sibling available', version : v})
      }
      
    }
    for (const root of nodes){
      if(!tree_id){
        const inDb = await Database.molecule.getVersionsFromAlias(infos.alias) //This is probably shit but my brain is tired and it works
        if(inDb.length > 0) tree_id = inDb[0].tree_id
      }
      let parentId; 
      const parentNumber = getParent(root.version.number)
      if(parentNumber){
        const parentDb = alreadyInDb.find(mol => mol.force_field === root.version.force_field && mol.version === parentNumber)
        if(!parentDb && root.version.number.split(".").length > 1){
          logger.warn(`${infos.alias} v${root.version.number} has no available parent version. Not inserted.`)
          not_inserted.push({reason : 'no parent available', version : root.version})
          continue
        }


        if(parentDb) parentId = parentDb.id
      }

      if(!getPreviousSiblingExistence(root.version.number, root.version.force_field, nodes, alreadyInDb)) {
        logger.warn(`${infos.alias} v${root.version.number} has no previous sibling version. Not inserted.`)
        not_inserted.push({version :root.version, reason : "no previous sibling available"})
        continue
      }

      try {
        await insert(root, infos, tree_id, parentId)
        inserted = inserted.concat(flatRoot(root))
      } catch(e) {
        reject(e)
      }


    }
  }

  resolve({inserted, not_inserted})

    
  })
    
}

const insert = async (node: MoleculeVersion, globalInfos: InfosJson, treeId?:string, parentId?: string)=> {
    const parent = await insertOne(globalInfos, node.version, treeId, parentId)
    const newParentId = parent.id
    let newTreeId = treeId; 
    if(!treeId) {
      const parentDoc = await Database.molecule.get(newParentId)
      newTreeId = parentDoc.tree_id
    }

    for(const child of node.children){
      await insert(child, globalInfos, newTreeId, newParentId).then(() => {})
    }
    
    return parent
  

}

const insertOne = async (globalInfos: InfosJson, version: VersionItp, treeId?: string, parentId?: string) => {
  
  const req : SimuRequest = {
    full_user: {
      id: CONNECTED_USER_CLI.id,
      role: CONNECTED_USER_CLI.role
    },
    body: {
      name: globalInfos.name,
      alias: globalInfos.alias,
      smiles: '',
      version: version.number,
      category: globalInfos.category,
      command_line: version.command_line,
      comments: version.comments,
      create_way: globalInfos.create_way,
      force_field: version.force_field,
      validation: '',
      citation: version.citation,
      parent: parentId ?? null,
      tree_id : treeId ?? undefined
    },
    files: {
      itp: [version.itp,],
      pdb: [version.gro,],
      top: [version.top,],
      map: version.map
    }

    
  }

  if (req.full_user.role != 'admin') {
    return Errors.throw(ErrorType.Unallowed);
  }

  // Error if the user connected is not an admin
  if (req.full_user.role != 'admin') {
    return Errors.throw(ErrorType.Unallowed);
  }

  const checker = new MoleculeChecker(req);
  const molecule = await checker.check()
  return await Database.molecule.save(molecule as Molecule);
}

const getParent = (version: string) => {
  const _versionSplit = version.split(".")
  const versionSplit = _versionSplit.slice(-1)[0] === '0' ? _versionSplit.slice(0,_versionSplit.length - 1): _versionSplit
  if(versionSplit.length === 1) return null
  const parentVersion = versionSplit.slice(0, versionSplit.length -1)
  if (parentVersion.length === 1) return parentVersion[0] + ".0"
  return parentVersion.join(".")
}

const getPreviousSibling = (version : string) => {
  const _versionSplit = version.split(".")
  const versionSplit = _versionSplit.slice(-1)[0] === '0' ? _versionSplit.slice(0,_versionSplit.length - 1): _versionSplit
  const last = versionSplit.slice(-1)[0]
  const newLast = parseInt(last) - 1
  if(newLast === 0) {
    return
  }
  const newVersion = versionSplit.length > 1 ? versionSplit.slice(0, versionSplit.length - 1).join(".") + "." + newLast.toString() : newLast.toString() + ".0"
  return newVersion
}

const getPreviousSiblingExistence = (version: string, force_field: string, treeVersions: MoleculeVersion[], alreadyInDb?: Molecule[]) => {
  const previousSibling = getPreviousSibling(version)
  if(previousSibling) {
    const siblingDb = alreadyInDb ? alreadyInDb.find(mol => mol.force_field === force_field && mol.version === previousSibling) : alreadyInDb
    if(!siblingDb) {
      const siblingRoot = treeVersions.find(root => root.version.number === previousSibling)
      if(!siblingRoot) {
        logger.warn(`v${version} has no previous sibling version. Not inserted.`)
        return false
      }
    }
  }
  return true
}

const getAllPreviousSiblings = (version: string) => {
  const _versionSplit = version.split(".")
  const versionSplit = _versionSplit.slice(-1)[0] === '0' ? _versionSplit.slice(0,_versionSplit.length - 1): _versionSplit
  const last = versionSplit.slice(-1)[0]
  const newLast = parseInt(last) - 1
  if(newLast === 0) {
    return
  }
  const _allLast = [...Array(newLast).keys()];
  const allLast = _allLast.slice(1, _allLast.length)
  const newVersion = versionSplit.length > 1 ? allLast.map(l => versionSplit.slice(0, versionSplit.length - 1).join(".") + "." + l.toString()) : allLast.map(l => newLast.toString() + ".0")
  return newVersion
}

const versionTree = (versions: VersionItp[]) => {
  const allNodes: MoleculeVersion[] = versions.map(v => ({version:v, children:[], root:true}))
  const not_treated: VersionItp[] = []
  
  for (const n of allNodes){
    const previousSiblings = getAllPreviousSiblings(n.version.number)
    if(previousSiblings && previousSiblings.length > 0) {
      const inTree = allNodes.filter(n => previousSiblings?.includes(n.version.number))
      if(inTree.length != previousSiblings.length)
        not_treated.push(n.version)
        n.root = false; 
        continue
    }
    const parentNumber = getParent(n.version.number)
    if(parentNumber){
      const parentNode = allNodes.find(node => node.version.number === parentNumber)
      if(parentNode){
        parentNode.children.push(n)
        n.root = false
      }
    }
  }

  return {nodes: allNodes.filter(n => n.root), notTreated: not_treated}

}

function separateVersionsByForcefield(versions: VersionItp[]) {
  const sorted : {[ff: string]: VersionItp[]} = {}
  for (const v of versions){
    if(!(v.force_field in sorted)) sorted[v.force_field] = []
    sorted[v.force_field].push(v)
  }
  return sorted
}

const flatRoot = (root: MoleculeVersion): VersionItp[] => {
  const _flatRootRec = (root: MoleculeVersion) => {
    all.push(root.version)
    for (const child of root.children) _flatRootRec(child)
  }

  const all: VersionItp[] = []
  _flatRootRec(root)
  return all
  
}

