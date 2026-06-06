import os
import requests

def fetch_and_clean_alphafold(uniprot_id: str, output_dir: str = "data") -> str:
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"AF_{uniprot_id}.cif")
    
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(api_url)
    
    if response.status_code != 200:
        raise ValueError(f"Structure not found in AlphaFold Database.")
        
    cif_url = response.json()[0]['cifUrl'] 
    cif_response = requests.get(cif_url)

    cif_response = requests.get(cif_url)
    if not cif_response.ok:
        raise ConnectionError("Error downloading CIF file from AlphaFold servers.")
        
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(cif_response.text)
        
    return output_path