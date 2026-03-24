from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif

if __name__ == '__main__':
    cif_path = fetch_cif("4HHB")
    entry = Entry(cif_path, "4HHB")
    entry.find_pockets()
    # entry.save_pocket_cif_files()
    metadata = entry.extract_metadata()


