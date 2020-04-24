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
  /** Molecule short alias */
  alias: string;
  /** Mol formula */
  formula: string;
  /** Category, should be a GO Term */
  category: keyof typeof GoTerms;
  /** Molecule version (free text) */
  version: string;
  /** Free comment text. */
  comments: string;
  /** String version of the used command line (other parameters) */
  command_line: string;
  /** Martinize version used to generate files */
  martinize_version: string;
  /** Force field version */
  force_field: string;
  /** Stringified ISO date of creation date */
  created_at: string;
  /** Molecule parent version. If string, ref to <Molecule.id> */
  parent: null | string;
  /** Tree snowflake ID. Shared between parent and children */
  tree_id: string;
  /** Hash of generated zip file attached to this module */
  hash: string;
  /** Reference to <User.id> owner/curator of this mol */
  owner: string;
  /** ID of related file containing `.itp` and `.gro`/`.pdb` files */
  files: string;
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
  /** Is approved or not */
  approved: boolean;
  /** Lost token ID */
  lost_token?: string;
}

export interface Token extends BaseCouchDocument {
  /** JTI UUID snowflake */
  id: string;
  /** <User.id> who own this token */
  user_id: string;
  /** Stringified ISO date of the token creation */
  created_at: string;
}

// extract all lennar johns > grep -P "^\s*(?'g1'\S+)\s+(\k{g1})\b" kwalp/martini_v.3.0.4.26/martini_v3.0.4.itp | sort | cut -d' ' -f-2,4- | uniq
export interface VanDerWaalsRadius extends BaseCouchDocument {
  /** Martinize force field related to the defined radius. It is a reference to _id, when it's present. */
  id: string;
  /** Map between atom name > van der waals radii. */
  atoms: {
    /** Van der Waals Radii */
    [name: string]: number
  };
}

export interface Lipid extends BaseCouchDocument {
  /** Lipid id. DO NOT USE. */
  id: string;
  /** Lipid short name. */
  name: string;
  /** Content of the ITP file for this lipid. */
  itp: string;
}

