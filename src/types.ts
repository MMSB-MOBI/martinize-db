
export const GoTerms = {
  "GO:0001": "<Go Term Regular Name>"
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
  martinize_versions: string[];
  category_tree: CategoryTree;
}

export interface CategoryTree {
  [go_id: string]: {
    children: CategoryTree,
    name: string;
  };
}
