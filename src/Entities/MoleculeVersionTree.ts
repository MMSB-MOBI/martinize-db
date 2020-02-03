import { BaseMolecule } from "./entities";

export default class MoleculeVersionTree {
  protected _root: MoleculeNode<BaseMolecule>;

  constructor(molecules: BaseMolecule[]) {
    if (!molecules.length) {
      throw new Error("Tree can't be empty");
    }

    let root: BaseMolecule | undefined;
    let tid: string = "";
    
    for (const mol of molecules) {
      if (!tid) {
        tid = mol.tree_id;
      }

      if (tid !== mol.tree_id) {
        throw new Error("Molecules should have the same tree id in order to construct tree.");
      }
      if (mol.parent === null) {
        if (typeof root === 'undefined') {
          root = mol;
        }
        else {
          throw new Error("Tree have multiple roots. Have you inserted the molecules properly ?");
        }
      }
    }

    if (!root) {
      throw new Error("Tree doesn't have root: The original parent must be present in molecule list. Maybe first element of tree has been removed ?");
    }

    this._root = new MoleculeNode(root);

    const inserted = new Set([root.id]);
    let to_insert = molecules.filter(m => !inserted.has(m.id));
    let last_to_insert = molecules;

    while (to_insert.length !== last_to_insert.length) {
      for (const mol of to_insert) {
        const appended = this._root.appendAt(mol, el => el.id === mol.parent);
        
        if (appended) {
          inserted.add(mol.id);
        }
      }

      last_to_insert = to_insert;
      to_insert = molecules.filter(m => !inserted.has(m.id));
    }
  }

  get root() {
    return this._root;
  }
}

export class MoleculeNode<T> {
  protected _content: T;
  protected _children: MoleculeNode<T>[] = [];

  constructor(molecule: T) {
    this._content = molecule;
  }

  append(molecule: T) {
    const mol = new MoleculeNode(molecule);
    this._children.push(mol);
    return mol;
  }

  appendAt(molecule: T, predicate: (el: T) => boolean) : MoleculeNode<T> | undefined {
    if (predicate(this._content)) {
      return this.append(molecule);
    }

    for (const child of this._children) {
      const inserted = child.appendAt(molecule, predicate);
      if (inserted) {
        return inserted;
      }
    }
  }

  find(predicate: (el: T) => boolean) : MoleculeNode<T> | undefined {
    if (predicate(this._content)) {
      return this;
    }

    for (const child of this._children) {
      const validated = child.find(predicate);
      if (validated) {
        return validated;
      }
    }
  }

  get content() {
    return this._content;
  }

  get children() {
    return this._children;
  }
}

export { MoleculeVersionTree };
