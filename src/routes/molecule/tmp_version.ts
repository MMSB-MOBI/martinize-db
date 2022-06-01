import { Database } from '../../Entities/CouchHelper';

export const correctVersions = async () => {
    const mols = await Database.molecule.all()
    for(const mol of mols){
        if(mol.version === '01'){
            mol.version = "1.0"
            mol.last_update = new Date().toISOString();
            await Database.molecule.save(mol)
        }
        else if(mol.version === '02'){
            mol.version = "2.0"
            mol.last_update = new Date().toISOString();
            await Database.molecule.save(mol)
        }
    }
}