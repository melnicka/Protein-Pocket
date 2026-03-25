import requests
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
        refs = data.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
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
        name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", None)
        if translate and name:
            try:
                name = GoogleTranslator(source='en', target='pl').translate(name)
            except Exception as e:
                print(f"Tłumaczenie nazwy białka nie powiodło się: {e}")
        return name
    except:
        pass
    return None

def get_function_from_uniprot(accession, translate=False):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        descriptions = []
        has_function = False
        cofactor_names = set()

        for comment in data.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    raw_text = texts[0].get("value", "")
                    cleaned_text = re.sub(r'\(PubMed:[^)]+\)', '', raw_text)
                    descriptions.append(cleaned_text.strip())
                    has_function = True
            if comment.get("commentType") == "COFACTOR":
                cofactors = comment.get("cofactors", [])
                for cofactor in cofactors:
                    name = cofactor.get("name")
                    if name:
                        cofactor_names.add(name)

        # Dodajemy kofaktory jako osobne linie z myślnikiem
        if cofactor_names:
            descriptions.append(f"Kofaktory: {', '.join(sorted(cofactor_names))}")

        # Dodajemy cechy jako osobne linie z myślnikiem, ale pomijamy 'Disordered'
        for feature in data.get("features", []):
            feature_type = feature.get("type", "").lower()
            if feature_type in {"domain", "region", "motif"}:
                desc = feature.get("description")
                if desc != "Disordered":
                    descriptions.append(f"Cecha: {desc}")

        if not has_function:
            descriptions.insert(0, "Ogólna funkcja nie została opisana w UniProt.")

        full_description = "\n".join(descriptions) if descriptions else None
        if translate and full_description:
            try:
                full_description = GoogleTranslator(source='en', target='pl').translate(full_description)
            except Exception as e:
                print(f"Tłumaczenie opisu funkcji nie powiodło się: {e}")

        return full_description

    except:
        return None
