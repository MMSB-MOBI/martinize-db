import { GoTerms, UserRole } from '../types';

export interface BaseCouchDocument {
  _id?: string;
  _rev?: string;
}

export interface BaseMolecule extends BaseCouchDocument {
  /** Molecule snowflake ID */
  id: string;
  /** Molecule name (free text) */
  name: string;
  /** Molecule version (free text) */
  version: string;
  /** Category, should be a GO Term */
  category: keyof typeof GoTerms;
  /** Molecule parent version. If string, ref to <Molecule.id> */
  parent: null | string;
  /** Tree snowflake ID. Shared between parent and children */
  tree_id: string;
  /** Hash of generated zip file attached to this module */
  hash: string;
  /** Reference to <User.id> owner/curator of this mol */
  owner: string;
  /** Filepath to ZIP file containing `.itp` and `.gro`/`.pdb` files */
  files: string;
  /** Free comment text. */
  comments: string;
  /** Stringified ISO date of creation date */
  created_at: string;
  /** String version of the used command line (other parameters) */
  command_line: string;
  /** Martinize version used to generate files */
  martinize_version: string;
  /** Force field version */
  force_field: string;
}

export interface Molecule extends BaseMolecule {
  /** <User.id> that have approved this molecule */
  approved_by: string;
  /** Last time as ISO date the user/admin edited this molecule */
  last_update: string;
}

export interface StashedMolecule extends BaseMolecule {}

export interface User extends BaseCouchDocument {
  /** User snowflake ID */
  id: string;
  /** User unique e-mail address */
  email: string;
  /** Display name */
  name: string;
  /** Stringified ISO Date of the user creation */
  created_at: string;
  /** bcrypt-hashed password of the user */
  password: string;
  /** User role */
  role: UserRole;
}

export interface Token extends BaseCouchDocument {
  /** JTI UUID snowflake */
  id: string;
  /** <User.id> who own this token */
  user_id: string;
  /** Stringified ISO date of the token creation */
  created_at: string;
}

