import requests
import os
from Bio import Entrez
from Bio import Medline
from biotite.structure.io.pdbx import CIFFile

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
    
def extract_protein_name(cif_path: str, pdb_id: str = None) -> str:
    """
    Robust protein name extraction from CIF with fallback.
    """
    try:
        cif = CIFFile.read(cif_path)
        block = list(cif.values())[0]
        title = block.get("_struct.title")
        if title is not None:
            return str(title.as_item())
        entity = block.get("_entity_poly.pdbx_strand_id")
        if entity is not None:
            return f"Protein {entity.as_item()}"
    except Exception:
        pass
    return pdb_id if pdb_id else "unknown protein"
def fetch_protein_papers(protein_name: str, max_results: int = 5) -> list[list[str]]:
    """
    Fetch top PubMed papers (more relevant + widely used results).
    Returns:
        [[link, abstract], ...]
    """
    search_term = protein_name.split("STRUCTURE OF")[-1].strip()
    search_handle = Entrez.esearch(
        db="pubmed",
        term=search_term,
        retmax=20,         
        sort="relevance")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    paper_ids = search_results["IdList"]
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
        pmid = record.get("PMID", "")
        title = record.get("TI", "")
        abstract = record.get("AB", "No abstract available.")
        if not pmid:
            continue
        papers.append([
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            f"{title}\n\n{abstract[:250]}..."
        ])
    return papers[:max_results]

