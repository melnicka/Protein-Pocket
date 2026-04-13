import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import requests

from src.utils.data_fetching import fetch_cif
from src.utils.cif_parsing import extract_metadata
from src.utils.fetch_uniprot import get_uniprot_accession, get_protein_name_from_uniprot, get_function_from_uniprot
from src.engine.entry import Entry
from src.engine.descriptors import (
    calc_ligand_buried_surface, calc_sasa_protein, calc_gyration_radius,
    calc_amino_acid_composition, calc_instability_index, calc_aromaticity,
    calc_helix_fraction, calc_isoelectric_point, calc_pocket_centroid,
    calc_pocket_hydrophobicity, calc_dipole_moment
)
from src.engine.protein_visualization import visualize_structure
from src.utils.cif_parsing import extract_metadata

st.set_page_config(page_title="Protein Viewer", layout="wide")
st.title("🧬 3D Protein Viewer & Pocket Analysis")

st.sidebar.header("Search & Settings")
pdb_id = st.text_input(
        "Enter PDB ID:",
        "1J91"
    ).strip().upper()

pocket_radius = 5.0 

if pdb_id:
    try:
        # 1. Fetch PDB Data 
        with st.spinner(f"Fetching {pdb_id} from RCSB PDB..."):
            cif_path = fetch_cif(pdb_id)

            pdb_response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
            raw_pdb_string = pdb_response.text if pdb_response.status_code == 200 else None
            custom_metadata = {}
            if raw_pdb_string:
                custom_metadata = extract_metadata(raw_pdb_string) 


        #  2. Parse PDB Data and find pockets 
        with st.spinner("Analyzing structure and finding pockets..."):
            entry = Entry(cif_path, pdb_id)
            entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
            entry.save_pocket_cif_files()
            metadata = entry.extract_metadata()
 
            try: global_sasa = calc_sasa_protein(entry.atom_array)
            except: global_sasa = None
            
            try: rg = calc_gyration_radius(entry.atom_array)
            except: rg = None
            
            try: pi = calc_isoelectric_point(entry.atom_array)
            except: pi = None
            
            try: instability = calc_instability_index(entry.atom_array)
            except: instability = None
            
            try: helix_frac = calc_helix_fraction(entry.atom_array)
            except: helix_frac = None
            
            try: aa_comp = calc_amino_acid_composition(entry.atom_array)
            except: aa_comp = None

            try: aromaticity = calc_aromaticity(entry.atom_array)
            except: aromaticity = None

            try: 
                dipole_vec = calc_dipole_moment(entry.atom_array)
                dipole_mag = np.linalg.norm(dipole_vec)
            except: 
                dipole_mag = None

        # 3. Fetch UniProt Data 
        with st.spinner("Fetching UniProt annotations..."):
            accession = get_uniprot_accession(pdb_id)
            prot_name, prot_func = None, None

        # Sidebar Controls 
        available_chains = list(metadata['chains'])
        selected_chains = st.sidebar.multiselect("Visible chain(s)", available_chains, default=available_chains)
        unique_ligand_ids = list(set([lig.comp_id for lig in entry.ligands]))

        st.sidebar.markdown("### 🎨 Chain colors")
        chain_colors = {}
        for chain in available_chains:
            picked = st.sidebar.color_picker(f"Chain {chain}", value="#FFFFFF")
            chain_colors[chain] = "spectrum" if picked.upper() == "#FFFFFF" else picked

        with st.sidebar:
            bg_color = st.selectbox("Background", options=["white", "black", "lightgray"], index=0)
            #pocket_radius = st.sidebar.slider("Pocket Search Radius (Å)", 3.0, 10.0, 5.0)

        col1, col2 = st.columns([1, 1])

        # Left Column: Metadata & Ligands
        with col1:
            st.subheader(f"Structure information")

            if custom_metadata:
            
                st.markdown(f"""
                **Type:** {custom_metadata.get('header', 'N/A')}  \n
                - **Resolution:** {custom_metadata.get('resolution', 'N/A')} Å \n
                - **Experiment:** {custom_metadata.get('experiment', 'N/A')}  \n
                - **Organism:** {custom_metadata.get('organism', 'N/A')}  \n
                """)
            
            if accession:
                st.success(f"**UniProt Accession:** [{accession}](https://www.uniprot.org/uniprotkb/{accession}/entry)")
            else:
                st.warning("No UniProt linkage found for this structure.")
            
            st.divider() 

            st.write("**Whole Protein Properties:**")
            prop_col1, prop_col2, prop_col3, prop_col4 = st.columns(4)
            if global_sasa is not None:
                prop_col1.metric("Total SASA", f"{float(global_sasa):.0f} Å²")
            if rg is not None: 
                prop_col2.metric("Gyration Rad.", f"{float(rg):.2f} Å")
            if pi is not None: 
                prop_col3.metric("pI", f"{float(pi):.2f}")
            if aromaticity is not None: 
                prop_col4.metric("Aromaticity", f"{(float(aromaticity)*100):.1f}%")
            
            prop_col5, prop_col6, prop_col7, prop_col8 = st.columns(4)
            if instability is not None: 
                stability_status = "🔴" if instability > 40 else "🟢"
                prop_col5.metric("Instability", f"{float(instability):.1f}", delta=stability_status, delta_color="off")
            if helix_frac is not None:
                prop_col6.metric("Helix Frac.", f"{(float(helix_frac)*100):.1f}%")
            prop_col7.metric("Total Atoms", metadata['full_atom_count'])
            if dipole_mag is not None:
                prop_col8.metric("Dipole Moment", f"{float(dipole_mag):.1f}")
            
            if aa_comp:
                with st.expander("Amino Acid Composition (%)"):
                    aa_df = pd.DataFrame(list(aa_comp.items()), columns=['Amino Acid', 'Percentage']).set_index('Amino Acid')
                    st.bar_chart(aa_df)

            st.divider()

            st.subheader("Ligands & Binding Pockets")
            if len(entry.ligands) > 0:
                for lig in entry.ligands:
                    with st.expander(f"Ligand {lig.comp_id} (Chain {lig.auth_asym_id}, Seq {lig.auth_seq_id})"):
                        st.write(f"**SMILES:** `{lig.smiles}`")
                        if not lig.pocket.is_empty:
                            try:
                                buried_surf = calc_ligand_buried_surface(lig.pocket.atom_array, lig.atom_array)
                                hydrophobicity = calc_pocket_hydrophobicity(lig.pocket.atom_array)
                                st.write(f"- **Buried Surface:** {buried_surf:.2f} Å²")
                                st.write(f"- **KD Hydrophobicity:** {hydrophobicity:.2f}")
                            except Exception as e:
                                st.warning(f"Could not calculate pocket descriptors: {e}")
                        else:
                            st.warning("No pocket residues found.")
            else:
                st.write("No ligands found.")

        # Right Column: Native 3D Viewer
        with col2:
            st.subheader("3D Structure Viewer")
            
            view_options = ["Whole Protein"]
            pocket_ligands = [lig for lig in entry.ligands if not lig.pocket.is_empty]
            for lig in pocket_ligands:
                view_options.append(f"Pocket: {lig.comp_id} ({lig.auth_asym_id}:{lig.auth_seq_id})")
                
            selected_view = st.selectbox("Select structure to render:", view_options)
            
            active_sites = None
            is_pocket = False
            
            if selected_view == "Whole Protein":
                file_to_load = cif_path
            else:
                is_pocket = True
                target_lig = next(lig for lig in pocket_ligands if f"Pocket: {lig.comp_id} ({lig.auth_asym_id}:{lig.auth_seq_id})" == selected_view)
                file_to_load = target_lig.pocket.cif_file_path
                active_sites = [{'chain': target_lig.auth_asym_id, 'resnum': target_lig.auth_seq_id}]

            with open(file_to_load, "r") as f:
                cif_data = f.read()
            
            view = visualize_structure(
                cif_data,
                bg_color=bg_color,
                selected_chains=selected_chains,
                chain_colors=chain_colors,
                active_site_residues=active_sites,
                is_pocket_view=is_pocket
            )
            
            components.html(view._make_html(), height=600, width=600)

    except Exception as e:
        st.error(f"An error occurred: {e}")
else:
    st.info("Please enter a PDB ID to get started.")
