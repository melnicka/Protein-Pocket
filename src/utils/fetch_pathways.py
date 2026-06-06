import requests
from src.utils.fetch_uniprot import get_uniprot_accession


def get_pathways(id, is_uniprot=False):
    if is_uniprot:
        uniprot_id = id
    else:
        uniprot_id = get_uniprot_accession(id)
    if not uniprot_id:
        return []
    uni_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(uni_url)
    if response.status_code != 200:
        return []
    uni_data = response.json()
    pathways = []
    for ref in uni_data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "KEGG":
            kegg_gene_id = ref.get("id")
            path_url = f"https://rest.kegg.jp/link/pathway/{kegg_gene_id}"
            path_response = requests.get(path_url)
            if path_response.status_code != 200:
                continue
            for line in path_response.text.splitlines():
                parts = line.split()
                if len(parts) < 2:
                    continue
                pathway_id = (
                    parts[1]
                    .replace("pathway:", "")
                    .replace("path:", "")
                )
                name_url = f"https://rest.kegg.jp/get/{pathway_id}"
                name_response = requests.get(name_url)
                pathway_name = pathway_id
                pathway_class = ""
                if name_response.status_code == 200:
                    for name_line in name_response.text.splitlines():
                        if name_line.startswith("NAME"):
                            pathway_name = (
                                name_line
                                .replace("NAME", "")
                                .strip()
                                .split(" - ")[0]
                            )
                        if name_line.startswith("CLASS"):
                            pathway_class = (
                                name_line
                                .replace("CLASS", "")
                                .strip()
                            )
                # tylko metaboliczne pathways
                if "Metabolism" not in pathway_class:
                    continue
                pathways.append({
                    "id": pathway_id,
                    "name": pathway_name,
                    "url": (
                        f"https://www.kegg.jp/kegg-bin/"
                        f"show_pathway?{pathway_id}+{kegg_gene_id}"
                    )
                })
    return pathways
