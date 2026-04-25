import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import requests

from util_components import pocket_radius, protein_details, structure_viewer
from src.utils.data_fetching import fetch_cif
from src.utils.cif_parsing import extract_metadata
from src.utils.fetch_uniprot import get_similar_proteins, get_uniprot_accession, get_protein_name_from_uniprot, display_aa_seq
from src.engine.entry import Entry
from src.engine.descriptors import (
    format_descriptor, calc_ligand_buried_surface, calc_sasa_protein, calc_gyration_radius,
    calc_amino_acid_composition, calc_instability_index, calc_aromaticity,
    calc_helix_fraction, calc_isoelectric_point, calc_pocket_centroid,
    calc_pocket_hydrophobicity, calc_dipole_moment, calc_hydrogen_bond_features, calc_charged_surface_fraction
)
from src.engine.protein_visualization import visualize_structure

if __name__ == "__main__":
    # main page
    st.set_page_config(page_title="Protein Viewer", layout="wide")
    st.title("🧬 3D Protein Viewer & Pocket Analysis")

    st.sidebar.header("Search & Settings")
    pdb_id = st.text_input(
        "Enter PDB ID:",
        "1J91"
    ).strip().upper()

    if pdb_id:

        cif_path = fetch_cif(pdb_id)
        entry = Entry(cif_path, pdb_id)
        entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
        entry.save_pocket_cif_files()
        metadata = entry.extract_metadata()

        col1, col2 = st.columns([1, 1])
        with col1:
            protein_details(pdb_id, entry)

        # Sidebar Controls
        available_chains = list(metadata.get('chains', []))
        selected_chains = st.sidebar.multiselect("Visible chain(s)", available_chains, default=available_chains)
        unique_ligand_ids = list(set([lig.comp_id for lig in entry.ligands]))

        st.sidebar.markdown("### 🎨 Chain colors")
        chain_colors = {}
        for chain in available_chains:
            picked = st.sidebar.color_picker(f"Chain {chain}", value="#FFFFFF")
            chain_colors[chain] = "spectrum" if picked.upper() == "#FFFFFF" else picked

        with st.sidebar:
            bg_color = st.selectbox("Background", options=["white", "black", "lightgray"], index=0)

        with col2:
            st.subheader("3D Structure Viewer")
            structure_viewer(entry, cif_path, bg_color, selected_chains, chain_colors)
    else:
        st.info("Please enter a PDB ID to get started.")
