import { Request } from 'express';
import { BaseMolecule, Molecule, StashedMolecule } from '../../Entities/entities';
import Errors, { ErrorType, ApiError } from '../../Errors';
import { MAX_ITP_FILE_SIZE, NAME_REGEX, ALIAS_REGEX, VERSION_REGEX } from '../Uploader';
import { generateSnowflake } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import { SETTINGS_FILE } from '../../constants';
import { promises as FsPromise } from 'fs';
import { SettingsJson, CategoryTree } from '../../types';
import MoleculeOrganizer, { MoleculeSave } from '../../MoleculeOrganizer';
import logger from '../../logger';


import { SimuFile } from '../../MoleculeOrganizer';
import { GoTerms } from '../../types';

export interface SimuRequest{
  full_user: {
    id: string,
    role: "admin"
  };
  body: {
    name: string,
    alias: string,
    smiles: string,
    version: string,
    category: keyof typeof GoTerms,
    command_line: string,
    comments: string,
    create_way: string,
    force_field: string,
    validation: string,
    citation: string,
    parent: string | null
  };
  files: {
    itp: SimuFile[],
    pdb: SimuFile[],
    top: SimuFile[] | [],
    map: SimuFile[] | []
  };
}



type File = Express.Multer.File;

export class MoleculeChecker {
  constructor(protected req: Request | SimuRequest) {}

  /**
   * Check a molecule about to be published to {Molecule} database
   */
  public async check() {
    const molecule = await this.checker(false, false) as Molecule;

    molecule.last_update = molecule.created_at;
    molecule.approved_by = this.req.full_user!.id;

    return molecule;
  }

  /**
   * Check a molecule about to be edited in {Molecule} database
   */
  public async checkEdition() {
    const molecule = await this.checker(false, true) as Molecule;

    molecule.last_update = new Date().toISOString();
    molecule.approved_by = this.req.full_user!.id;

    return molecule;
  }

  /**
   * Check a molecule about to be published to {Stashed} database
   */
  public async checkStashed() {
    const molecule = await this.checker(true, false) as StashedMolecule;
    
    return molecule;
  }

  /**
   * Check a molecule about to be edited in {Stashed} database
   */
  public async checkStashedEdition() {
    const molecule = await this.checker(true, true) as StashedMolecule;

    return molecule;
  }

  protected async checker(stashed: boolean, edition: boolean) {
    let actual_version: BaseMolecule | undefined = undefined;

    
    if (edition) {
      const id = this.req.body.id;
      if (!id) {
        return Errors.throw(ErrorType.MissingParameters);
      }
      
      try {
        if (stashed) {
          actual_version = await Database.stashed.get(id);
        }
        else {
          actual_version = await Database.molecule.get(id);
        }
      } catch (e) {
        return Errors.throw(ErrorType.MoleculeNotFound);
      }
    }
    
    const molecule = await this.constructBaseMoleculeFromRequest(actual_version);

    // Must set the following fields: files, created_at, hash

    // Test if files are attached to the request
    if (!this.areFilesPresent()) {
      if (!edition) {
        
        return Errors.throw(ErrorType.MissingParameters);
      }
      else {
        // Edition mode, we must assure that files is defined
        if (!molecule.files) {
          return Errors.throw(ErrorType.MissingParameters);
        }

        // Check if the file exists
        if (!(await MoleculeOrganizer.exists(molecule.files))) {
          return Errors.throw(ErrorType.MissingFiles);
        }

        // Refresh the hash
        molecule.hash = await MoleculeOrganizer.hash(molecule.files);
      }
    }
    // Must insert files into upload directory
    else {
      const files = await this.getFilesFromRequest();

      // Register the files
      // Save the molecule in ZIP format if files changed
      let save: MoleculeSave;
      try {
        save = await MoleculeOrganizer.save(
          files.itps, 
          files.molecule, 
          files.top, 
          files.maps,
          molecule.force_field!
        );
      } catch (e) {
        if (e instanceof ApiError) {
          throw e;
        }
        if (e instanceof Error) {
          console.log(files.itps);
          return Errors.throw(ErrorType.InvalidMoleculeFiles, { detail: e.message });
        }
        return Errors.throw(ErrorType.InvalidMoleculeFiles);
      }

      // Remove the old files, if they exists
      if (molecule.files) {
        try {
          await MoleculeOrganizer.remove(molecule.files);
        } catch (e) {
          logger.error("Unable to remove file", e);
        }
      }

      // Register hash + file id
      molecule.hash = save.infos.hash;
      molecule.files = save.id;
    }

    if (!molecule.created_at) {
      molecule.created_at = new Date().toISOString();
    }

    // Every field is set (except molecule specific fields) !
    return molecule as BaseMolecule;
  }

  protected areFilesPresent() {
    if (!this.req.files || Array.isArray(this.req.files)) {
      return false;
    }

    // Find the ITP files
    const itps_files: (File | SimuFile)[] = this.req.files.itp;
    const pdb_files: (File | SimuFile)[] = this.req.files.pdb;
    const top_files: (File | SimuFile)[] = this.req.files.top;
    

    // Requires one itp file at least
    if (!itps_files || !top_files || !pdb_files || !itps_files.length || !pdb_files.length || !top_files.length) {
      return false;
    }
    return true;
  }

  protected async getFilesFromRequest() {
    if (!this.req.files || Array.isArray(this.req.files)) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    // Find the ITP files
    const itps_files: (File | SimuFile)[] = this.req.files.itp;


    // Check the weight of each ITP file
    if (itps_files.some(f => f.size > MAX_ITP_FILE_SIZE)) {
      return Errors.throw(ErrorType.FileTooLarge);
    }

    // Get the PDB files
    const pdb_files: (File | SimuFile)[] = this.req.files.pdb;
    const is_pdb = !!this.req.files.pdb;

    if (pdb_files.length !== 1) {
      return Errors.throw(ErrorType.TooManyFiles);
    }

    const map_files: (File | SimuFile)[] = this.req.files.map || [];

    // Find if top files are present
    let top_file: File | SimuFile = this.req.files.top[0];

    return {
      molecule: pdb_files[0],
      is_pdb,
      itps: itps_files,
      top: top_file,
      maps: map_files,
    };
  }

  /**
   * Construct a base molecule that contain 
   * id, name, alias, formula, category, version, comments, command_line, martinize_version, force_field, parent, tree_id, owner.
   * 
   * MISSING : files, created_at, hash, <last_update>, <approved_by>
   */
  protected async constructBaseMoleculeFromRequest(actual_version?: BaseMolecule) : Promise<Partial<BaseMolecule>> {
    const nullOrString = (str?: string | null) => {
      if (!str || str === undefined || str === null || str === 'null') {
        return null;
      }
      return typeof str === 'string' ? str : undefined;
    };

    const mol: Partial<BaseMolecule> = actual_version ? { ...actual_version } : {};
    let parent: Molecule | undefined = undefined;
    const settings: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));

    if (!mol.owner) {
      mol.owner = this.req.full_user!.id;
    }

    const body = this.req.body;


    // If the molecule doesn't have any ID, set it
    if (!mol.id) {
      mol.id = generateSnowflake();
    }

    // Download the parent if the molecule had one
    const parent_id = nullOrString(body.parent) || mol.parent;

    if (parent_id) {
      try {
        parent = await Database.molecule.get(parent_id);
        this.swallowCopyFromParent(mol, parent);
      } catch (e) {
        return Errors.throw(ErrorType.UnknownParent);
      }
    }

    // Set the right parent
    if (!mol.parent) {
      mol.parent = null;
      if (!mol.tree_id) {
        mol.tree_id = generateSnowflake();
      }
    }

    // Check the parent-linked fields (copy from them only if molecule is not parented)
    if (!parent) {
      if (!body.name || !body.alias || !body.category) {
        return Errors.throw(ErrorType.MissingParameters);
      }

      await this.checkName(body.name, mol.tree_id!);
      mol.name = body.name;

      await this.checkAlias(body.alias, mol.tree_id!);
      mol.alias = body.alias;

      // TODO introduce check
      mol.smiles = body.smiles || "";

      this.checkCategory(body.category, settings);
      mol.category = body.category;
    }

    // OK now: id, name, tree_id, alias, formula, category, parent, owner

    // Copy the version-specific fields
    if (!body.version || !body.force_field || !body.create_way) {
      return Errors.throw(ErrorType.MissingParameters);
    }

    body.version = body.version.trim();
    this.checkVersion(body.version);
    mol.version = body.version;

    // Check if version exists in tree_id
    if (await this.versionExistsInTreeId(mol.tree_id!, mol.version!, mol.id!)) {
      return Errors.throw(ErrorType.VersionAlreadyExists);
    }

    // Optional fields TODO limit length
    mol.comments = body.comments || "";
    mol.command_line = body.command_line || "";
    mol.citation = body.citation || "";
    mol.validation = body.validation || "";

    // Check force field and martinize version
    this.checkCreateWay(body.create_way, settings);
    this.checkForceField(body.force_field, settings);

    mol.force_field = body.force_field;
    mol.create_way = body.create_way;

    return mol;
  }

  protected swallowCopyFromParent(molecule: Partial<BaseMolecule>, parent: Molecule) {
    molecule.name = parent.name;
    molecule.alias = parent.alias;
    molecule.smiles = parent.smiles;
    molecule.category = parent.category;
    molecule.tree_id = parent.tree_id;
    molecule.parent = parent.id;
  }

  protected async checkName(name: string, tree_id: string) {
    if (!name.match(NAME_REGEX)) {
      return Errors.throw(ErrorType.InvalidName);
    }

    const mols = await Database.molecule.find({ limit: 99999, selector: { name } });

    for (const mol of mols) {
      if (mol.tree_id !== tree_id) {
        return Errors.throw(ErrorType.NameAlreadyExists);
      }
    }
  }

  protected async checkAlias(name: string, tree_id: string) {
    if (!name.match(ALIAS_REGEX)) {
      return Errors.throw(ErrorType.InvalidAlias);
    }

    const mols = await Database.molecule.find({ limit: 99999, selector: { alias: name } });

    for (const mol of mols) {
      if (mol.tree_id !== tree_id) {
        return Errors.throw(ErrorType.AliasAlreadyExists);
      }
    }
  }

  protected checkVersion(v: string) {
    if (!v.match(VERSION_REGEX)) {
      return Errors.throw(ErrorType.InvalidVersion);
    }
  }

  protected checkCategory(cat: string, settings: SettingsJson) {
    // todo
    const findInCategoryTree = (val: string, node: CategoryTree) : boolean => {
      for (const go_id in node) {
        if (go_id === val) {
          return true;
        }
        if (node[go_id].children && findInCategoryTree(val, node[go_id].children)) {
          return true;
        }
      }

      return false;
    };
    
    if (!findInCategoryTree(cat, settings.category_tree)) {
      return Errors.throw(ErrorType.InvalidCategory);
    }
  }

  protected checkCreateWay(v: string, settings: SettingsJson) {
    if (!(v in settings.create_way)) {
      return Errors.throw(ErrorType.InvalidMartinizeVersion);
    }
  }

  protected checkForceField(v: string, settings: SettingsJson) {
    if (!settings.force_fields.includes(v)) {
      return Errors.throw(ErrorType.InvalidForceField);
    }
  }

  protected async versionExistsInTreeId(tree_id: string, version: string, current_id: string) {
    const versions = await Database.molecule.find({ limit: 20, selector: { tree_id, version } });

    for (const v of versions) {
      if (v.id !== current_id) {
        return true;
      }
    }

    return false;
  }
}
