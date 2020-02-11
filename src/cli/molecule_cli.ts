import { CliListener } from "interactive-cli-helper";
import CouchHelper, { Database } from "../Entities/CouchHelper";
import MoleculeOrganizer from "../MoleculeOrganizer";
import { Molecule } from "../Entities/entities";

const MOLECULE_CLI = new CliListener("Available commands are list, get, wipe");

MOLECULE_CLI.addSubListener('list', async () => {
  const mols = await Database.molecule.all();

  if (!mols.length) {
    return `This server does not contain any molecule.`;
  }

  return `Available molecules are \n- ${mols.map(m => m._id).join('\n- ')}`;
});

MOLECULE_CLI.addSubListener('get', rest => {
  rest = rest.trim();

  if (!rest) {
    return `Please specify a molecule id.`;
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  return Database.molecule.get(rest);
});

MOLECULE_CLI.addSubListener('wipe', async rest => {
  rest = rest.trim();
  
  if (!rest) {
    return `Please specify a molecule id or "all".`;
  }

  if (rest === "all") {
    await Database.delete(CouchHelper.MOLECULE_COLLECTION);
    await MoleculeOrganizer.removeAll();
    await Database.create(CouchHelper.MOLECULE_COLLECTION);
    return `Molecule database is wiped`;
  }

  try {
    BigInt(rest);
  } catch (e) {
    return `ID ${rest} is not valid, please enter a valid number.`;
  }

  let mol: Molecule;
  try {
    mol = await Database.molecule.get(rest);
  } catch (e) {
    return `Unable to get molecule (${rest})`;
  }

  if (mol) {
    await MoleculeOrganizer.remove(mol.files);
    return Database.molecule.delete(mol);
  }
  return `Unable to find molecule.`
});

export default MOLECULE_CLI;
