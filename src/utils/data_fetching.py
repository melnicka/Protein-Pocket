import requests
import os

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
