from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif

# testing
if __name__ == '__main__':
    cif_path = fetch_cif("4HHB")
    entry = Entry(cif_path, "4HHB")
    entry.find_pockets()
    metadata = entry.extract_metadata()

    for ligand in entry.ligands:
        print(ligand.pocket.cif_file_path)

