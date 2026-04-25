import requests
import streamlit as st


def get_uniprot_accession(pdb_id):
    """
    Retrieves the UniProt accession number for a given PDB ID.

    Parameters
    ----------
    pdb_id : str
        The PDB identifier.

    Returns
    -------
    str or None
        The UniProt accession number if found, otherwise None.
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        refs = data.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "reference_sequence_identifiers", []
        )
        for ref in refs:
            if ref.get("database_name") == "UniProt":
                return ref.get("database_accession")
    except:
        pass
    return None


def get_protein_name_from_uniprot(accession, translate=False):
    """
    Retrieves the recommended protein name from UniProt using an accession number.

    Parameters
    ----------
    accession : str
        The UniProt accession number.
    Returns
    -------
    str or None
        The full recommended name of the protein, or None if not found.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        name = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", None)
        )
        return name
    except:
        pass
    return None


def get_similar_proteins(accession, level=90):
    """
    level = 100, 90, 50
    """
    params = {
        "query": f"(uniref_cluster_{level}:UniRef{level}_{accession}) NOT (accession:{accession}) AND (database:pdb)",
        "fields": ["accession", "protein_name", "organism_name"],
        "sort": "accession desc",
        "size": "10",
    }
    headers = {"accept": "application/json"}
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    try:
        response = requests.get(base_url, headers=headers, params=params)
        response.raise_for_status()
        data = response.json()
        return data["results"]
    except:
        pass

    return None

def get_pdb_from_uniprod_acc(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        xrefs = data["uniProtKBCrossReferences"]
        xrefs = filter(lambda x: x["database"] == "PDB", xrefs)
        pdb_xref = next(xrefs)

        return pdb_xref["id"]
    except:
        pass
    return None

def display_aa_seq(sequences: dict, chunk_size=60):
    expander_label = "📄 View FASTA sequence"

    fasta_output = ""

    for chain, seq in sequences.items():
        fasta_output += f">Chain_{chain}\n"
        for i in range(0, len(seq), chunk_size):
            fasta_output += seq[i : i + chunk_size] + "\n"

    with st.expander(expander_label):
        st.code(fasta_output.strip(), language="fasta")
