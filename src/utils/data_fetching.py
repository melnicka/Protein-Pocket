import requests
import os
from Bio import Entrez
from Bio import Medline
from biotite.structure.io.pdbx import CIFFile


# przydało by dać email ( lepsze rezultaty czy cus ) 
#Entrez.email = "test@example.com"
def fetch_cif(pdb_id: str, root_dir="data") -> str:
    """Downloads protein's 3D structure file in CIF format from RCSB PDB.
    Creates a unique directory for protein's data.
    If the file already exists, just returns the file path.

    Args:
        pdb_id: Unique PDB identifier.
        root_dir: Root data directory path.

    Returns:
        Path to the CIF file.

    Raises:
        RuntimeError: If failed to get the response from the server.
    """
    pdb_id = pdb_id.upper()
    dir_path = f"{root_dir}/{pdb_id}"
    file_path = f"{dir_path}/{pdb_id}.cif"

    if os.path.exists(file_path):
        return file_path

    os.makedirs(dir_path, exist_ok=True)

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try: 
        resp = requests.get(url)
        resp.raise_for_status()
    except requests.RequestException:
        raise RuntimeError(f"Failed to download {pdb_id}.")

    with open(file_path, "wb") as f:
        f.write(resp.content)

    return file_path
    
def extract_protein_name(cif_path: str, pdb_id: str) -> str:
    """
    #to niedziała narazie ( potrzebuje jakoś wyextractować nazwe białka )
    """
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        return data["struct"]["title"]
    except Exception:
        return pdb_id

#ta część działa 
def fetch_protein_papers(protein_name: str, max_results: int = 5) -> list[list[str]]:
    """
    Improved PubMed fetch with better query + safety.
    """
    clean_name = protein_name.replace("STRUCTURE OF", "").strip()
    search_term = f"{clean_name} protein structure"

    try:
        search_handle = Entrez.esearch(
            db="pubmed",
            term=search_term,
            retmax=20,
            sort="relevance"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        paper_ids = search_results.get("IdList", [])
        if not paper_ids:
            return []
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=",".join(paper_ids),
            rettype="medline",
            retmode="text"
        )
        records = list(Medline.parse(fetch_handle))
        fetch_handle.close()
        papers = []
        for record in records:
            pmid = record.get("PMID")
            title = record.get("TI") or "No title"
            abstract = record.get("AB") or "No abstract available."
            if not pmid:
                continue
            papers.append([
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                f"{title}\n\n{abstract[:250]}..."
            ])
        return papers[:max_results]
    except Exception as e:
        print(f"PubMed error: {e}")
        return []
