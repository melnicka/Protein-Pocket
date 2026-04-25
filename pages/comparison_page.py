import streamlit as st

from src.utils.fetch_uniprot import get_pdb_from_uniprod_acc
from util_components import pocket_radius, protein_details, structure_viewer
from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif

acc1 = st.query_params.get("acc1")
acc2 = st.query_params.get("acc2")

pdb1, pdb2 = None, None
with st.spinner("Fetching PDB ids"):
    pdb1 = get_pdb_from_uniprod_acc(acc1)
    pdb2 = get_pdb_from_uniprod_acc(acc2)

if not pdb1 or not pdb2:
    st.warning("Please specify proteins to compare")
else:
    st.set_page_config(page_title=f"{pdb1} vs {pdb2}", layout="wide")
    st.title(f"{pdb1} vs {pdb2}")

    col1, col2 = st.columns(2)

    with st.sidebar:
        bg_color = st.selectbox("Background", options=["white", "black", "lightgray"], index=0)

    with col1:
        cif_path = fetch_cif(pdb1)
        entry = Entry(cif_path, pdb1)
        entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
        entry.save_pocket_cif_files()
        metadata = entry.extract_metadata()
        available_chains = list(metadata.get('chains', []))
        selected_chains = st.sidebar.multiselect("Visible chain(s)", available_chains, default=available_chains, key="chains1")
        unique_ligand_ids = list(set([lig.comp_id for lig in entry.ligands]))

        st.sidebar.markdown("### 🎨 Chain colors")
        chain_colors = {}
        for chain in available_chains:
            picked = st.sidebar.color_picker(f"Chain {chain}", value="#FFFFFF", key=f"colapicka1_{chain}")
            chain_colors[chain] = "spectrum" if picked.upper() == "#FFFFFF" else picked

        structure_viewer(entry, cif_path, bg_color, selected_chains, chain_colors)
        protein_details(pdb1, entry)

    with col2:
        cif_path = fetch_cif(pdb2)
        entry = Entry(cif_path, pdb2)
        entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
        entry.save_pocket_cif_files()
        metadata = entry.extract_metadata()
        available_chains = list(metadata.get('chains', []))
        selected_chains = st.sidebar.multiselect("Visible chain(s)", available_chains, default=available_chains, key="chains2")
        unique_ligand_ids = list(set([lig.comp_id for lig in entry.ligands]))

        st.sidebar.markdown("### 🎨 Chain colors")
        chain_colors = {}
        for chain in available_chains:
            picked = st.sidebar.color_picker(f"Chain {chain}", value="#FFFFFF", key=f"colapicka2_{chain}")
            chain_colors[chain] = "spectrum" if picked.upper() == "#FFFFFF" else picked

        structure_viewer(entry, cif_path, bg_color, selected_chains, chain_colors)
        protein_details(pdb2, entry)
