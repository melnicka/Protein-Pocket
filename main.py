from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif

# testing
if __name__ == '__main__':
    cif_path = fetch_cif("9QFX")
    entry = Entry(cif_path, "9QFX")
    entry.find_pockets()
    entry.save_pocket_cif_files()
    metadata = entry.extract_metadata()

    for ligand in entry.ligands:
        print(ligand.smiles)
        print(ligand.pocket.atom_array.shape[0])

