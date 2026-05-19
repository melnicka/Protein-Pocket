from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif, extract_protein_name, fetch_protein_papers

if __name__ == '__main__':
    cif_path = fetch_cif("4HHB")
    entry = Entry(cif_path, "4HHB")
    entry.find_pockets()
    # entry.save_pocket_cif_files()
    metadata = entry.extract_metadata()
    # TEST PUBMED
    #protein_name = extract_protein_name(cif_path, "4HHB")
    papers = fetch_protein_papers("Human Deoxyhaemoglobin")
    #print(" PAPERS ")
    #for link, abstract in papers:
        #print(link)
        #print(abstract)
       # print("-" * 80)



