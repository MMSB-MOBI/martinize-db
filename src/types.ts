
export const GoTerms = {
  "MC:0001": "<Go Term Regular Name>",
  "MC:0005": "",
  "MC:0002":"",
  "MC:0003":"",
  "MC:0004":"",
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
  create_way: { [wayId: string]: string };
  category_tree: CategoryTree;
}

export interface CategoryTree {
  [go_id: string]: {
    children: CategoryTree,
    name: string;
  };
}
