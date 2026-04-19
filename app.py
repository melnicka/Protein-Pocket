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
    format_descriptor, calc_ligand_buried_surface, calc_sasa_protein, calc_gyration_radius,
    calc_amino_acid_composition, calc_instability_index, calc_aromaticity,
    calc_helix_fraction, calc_isoelectric_point, calc_pocket_centroid,
    calc_pocket_hydrophobicity, calc_dipole_moment, calc_hydrogen_bond_features, calc_charged_surface_fraction
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
            def safe_calc(func, *args):
                try: 
                    return func(*args)
                except Exception as e:
                    print(f"[ERROR] {func.__name__}: {e}")
                    return None

            global_sasa = safe_calc(calc_sasa_protein, entry.atom_array)
            rg = safe_calc(calc_gyration_radius, entry.atom_array)
            pi = safe_calc(calc_isoelectric_point, entry.atom_array)
            instability = safe_calc(calc_instability_index, entry.atom_array)
            helix_frac = safe_calc(calc_helix_fraction, entry.atom_array)
            aa_comp = safe_calc(calc_amino_acid_composition, entry.atom_array)
            aromaticity = safe_calc(calc_aromaticity, entry.atom_array)

            try: 
                dipole_vec = calc_dipole_moment(entry.atom_array)
                dipole_mag = np.linalg.norm(dipole_vec) if dipole_vec is not None else None
            except Exception as e:
                print(f"[ERROR] Dipole moment: {e}") 
                dipole_mag = None

            try: global_sasa = calc_sasa_protein(entry.atom_array)
            except Exception as e:
                print(f"[ERROR] SASA: {e}")
                global_sasa = None

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
            prop_col1.metric("Total SASA", format_descriptor(global_sasa, unit="Å²", decimals=0))
            prop_col2.metric("Gyration Rad.", format_descriptor(rg, unit="Å", decimals=2))
            prop_col3.metric("pI", format_descriptor(pi, decimals=2))
            prop_col4.metric("Aromaticity", format_descriptor(aromaticity, decimals=1, is_percent=True))
            
            prop_col5, prop_col6, prop_col7, prop_col8 = st.columns(4)
            stability_status = "🔴" if (instability and instability > 40) else ("🟢" if instability else None)
            prop_col5.metric("Instability", format_descriptor(instability, decimals=1), delta=stability_status, delta_color="off")
            prop_col6.metric("Helix Frac.", format_descriptor(helix_frac, decimals=1, is_percent=True))
            prop_col7.metric("Total Atoms", metadata.get('full_atom_count', 'N/A'))
            prop_col8.metric("Dipole Moment", format_descriptor(dipole_mag, decimals=1))
            
            if aa_comp:
                with st.expander("Amino Acid Composition (%)"):
                    aa_df = pd.DataFrame(list(aa_comp.items()), columns=['Amino Acid', 'Percentage']).set_index('Amino Acid')
                    st.bar_chart(aa_df)

            st.divider()

            if hasattr(entry, 'ligands') and len(entry.ligands) > 0:
                for lig in entry.ligands:
                    with st.expander(f"Ligand {lig.comp_id} (Chain {lig.auth_asym_id}, Seq {lig.auth_seq_id})"):
                        st.write(f"**SMILES:** `{getattr(lig, 'smiles', 'N/A')}`") 
            
                        if hasattr(lig, 'pocket') and not lig.pocket.is_empty:
                            try:
                                buried_surf = calc_ligand_buried_surface(lig.pocket.atom_array, lig.atom_array)
                                hydrophobicity = calc_pocket_hydrophobicity(lig.pocket.atom_array)
                                hbond_features = calc_hydrogen_bond_features(lig.pocket.atom_array)
                                charged_frac = calc_charged_surface_fraction(lig.pocket.atom_array)
                
                                st.write(f"- **Buried Surface:** {buried_surf:.2f} Å²")
                                st.write(f"- **KD Hydrophobicity:** {hydrophobicity:.2f}")
                                st.write(f"- **H-Bond Donors:** {hbond_features.get('HBD', 'N/A')}")
                                st.write(f"- **H-Bond Acceptors:** {hbond_features.get('HBA', 'N/A')}")
                                st.write(f"- **Charged Surface:** {(charged_frac * 100):.1f}%")
                    
                            except Exception as e:
                                st.warning(f"Could not calculate pocket descriptors: {e}")
                        else:
                            st.warning("No pocket residues found.")
            else:
                st.info("No ligands or pockets detected in this structure.")

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
