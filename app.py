import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd

from src.utils.data_fetching import fetch_cif
from src.utils.cif_parsing import load_cif, extract_entry_metadata, extract_ligands
from src.utils.fetch_uniprot import get_uniprot_accession, get_protein_name_from_uniprot, get_function_from_uniprot
from src.engine.entry import Entry
from src.engine.descriptors import (
    calc_ligand_buried_surface, calc_sasa_protein, calc_gyration_radius,
    calc_amino_acid_composition, calc_instability_index, calc_aromaticity,
    calc_helix_fraction, calc_isoelectric_point, calc_pocket_centroid,
    calc_pocket_hydrophobicity, calc_dipole_moment
)

st.set_page_config(page_title="Protein Viewer", layout="wide")
st.title("🧬 3D Protein Viewer & Pocket Analysis Dashboard")
st.sidebar.header("Search & Settings")
pdb_id = st.sidebar.text_input("Enter PDB ID:" "1HDA").strip().upper()
st.sidebar.subheader("Analysis Settings")

pocket_radius = 5.0

if pdb_id:
  try:
    with st.spinner(f"Fetching {pdb_id} from RCSB PDB..."):
      cif_path = fetch_cif(pdb_id)

    with st.spinner("Analyzing structure..."):
      entry = Entry(cif_path, pdb_id)
      entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
      entry.save_pocket_cif_files()
      metadata = entry.extract_metadata()

      # Calculate global descriptors
      try: global_sasa = calx_sasa_protein(entry.atom_array)
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

      try: aromacity = calc_aromacity(entry.atom_array)
      except: aromacity = None

      try: dipole_moment = calc_dipole_moment(entry.atom_array)
      except: dipole_moment = None

    with st.spinner("Fetching UniProt annotations..."):
      accession = get_uniprot_accession(pdb_id)
      prot_name = None
      prot_func = NOne
      if accession:
        prot_name - get_protein_name_from_uniprot(accession)
        prot_func = get_function_from_uniprot(accession)

    col1, col2 = st.columns([1,1])

    '''Left column for metadata and ligands''''
    with col1:
      st.subheader(f"Data for {pdb_id}")

      if accession:
        st.success(f"**UNiProt Accession:** [{accession}](https://www.uniprot.org/uniprotkb/{accession}/entry)")
        if prot_name:
          st.write(f"**Name:** {prot_name}")
        if prot_func:
          with st.expander("View Protein Function and Other Features", expanded=False):
            st.write(prot_func)

        else:
          st.warning("No UniProt linkage found for this structure.")

      st.divider

      st.write("**Whole Protein Properties:**")
      prop_col1, prop_col2, prop_col3, prop_col4 = st.columns(4)
      if global_sasa: prop_col1.metric("Total SASA", f"{global_sasa:.0f} Å²")
      if rg: prop_col2.metric("Gyration Rad.", f"{rg:.2f} Å")
      if pi: prop_col3.metric("pI", f"{pi:.2f}")
      if aromaticity: prop_col4.metric("Aromaticity", f"{(aromaticity*100):.1f}%")
            
      prop_col5, prop_col6, prop_col7, prop_col8 = st.columns(4)
      if instability: 
        stability_status = "🔴" if instability > 40 else "🟢"
        prop_col5.metric("Instability", f"{instability:.1f}", delta=stability_status, delta_color="off")
      if helix_frac: prop_col6.metric("Helix Frac.", f"{(helix_frac*100):.1f}%")
      prop_col7.metric("Total Atoms", metadata['full_atom_count'])
      if dipole_moment: prop_col8.metric("Dipole Moment", f"{dipole_moment:.1f}")

      if aa_comp:
        with st.expander("Amino Acid Composition (%)"):
          aa+df = pd.DataFrame(list(aa_comp.items()), columns=['Amino Acid', 'Percentage']).set_index('Amino Acid')
          st.bar_chart(aa_df)
          
      export_dict = {
        "PDB_ID": [pdb_id],
        "Total_Atoms": [metadata['full_atom_count']],
        "SASA": [global_sasa],
        "Radius_of_Gyration": [rg],
        "Isoelectric_Point": [pi],
        "Aromaticity_Pct": [aromaticity * 100 if aromaticity else None],
        "Instability_Index": [instability],
        "Helix_Fraction_Pct": [helix_frac * 100 if helix_frac else None],
        "Dipole_Moment": [dipole_moment]
      }

      export_df = pd.DataFrame(export_dict)
      csv_data = export_df.to_csv(index=False).encode('utf-8')

      st.download_button(
        label="📥 Download Metrics as CSV",
        data = csv_data,
        file_name = f"{pdb_id}_metrics.csv",
        mime = "text/csv",
        use_container_width = True
      )

      st.divider()
      
      st.subheader("Ligands & Binding Pockets")
      if len(entry.ligands) >0:
        for lig in entry.ligands:
          with st.expander(f"Ligand {lig.comp_id} (Chain , Seq)"):
            st.write(f'**SMILES** `{lig.smiles}`')
            
            if not lig.pocket.is_empty:
              st.success(f"Pocket found within {pocket_radius}Å")

              #Calculate Protein Descriptors
              try:
                buried_surf = calc_ligand_buried_surface(lig.pocket.atom_array, lig.atom_array)
                hydrophobicity = calc_pocket_hydrophobicity(lig.pocket.atom_array)
                centroid = calc_pocket_centroid(lig.pocket.atom_array)

                st.write(f"- **Buried Surface:** {buried_surf:.2f} Å²")
                st.write(f"- **KD Hydrophobicity:** {hydrophobicity:.2f}")
                st.write(f"- **Centroid (x,y,z):** ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
              except Exception as e:
                st.warning(f"Could not calculate all pocket descriptors: {e}")

            else: 
              st.warning("No pocket residues found.")

        else:
          st.write("No ligands found in this structure.")
          
      '''Right column: 3D Viewer'''
      with col2: 
        st.subheader("3D Structure Viewer")

        view_options = ["Whole Protein"]
        pocket_ligands = [lig for lig in entry.ligands if not lig.pocket.is_empty]
        
        for lig in pocket_ligands:
          view_options.append(f"Pocket: {lig.comp_id}")

        selected_view = st.selectbox("Select structure to render:", view_options)
            
            if selected_view == "Whole Protein":
                file_to_load = cif_path
                style_func = lambda v: (v.setStyle({'cartoon': {'color': 'spectrum'}}), v.addStyle({'hetflag': True}, {'stick': {}}))
            else:
                target_lig = next(lig for lig in pocket_ligands if f"Pocket: {lig.comp_id} ({lig.auth_asym_id}:{lig.auth_seq_id})" == selected_view)
                file_to_load = target_lig.pocket.cif_file_path
                style_func = lambda v: (
                    v.setStyle({'stick': {'colorscheme': 'grayCarbon'}}),
                    v.addStyle({'resn': target_lig.comp_id}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
                )

            with open(file_to_load, "r") as f:
              cif_data = f.read()

            view = py3Dmol.view(width=500, height=500)
            view.addModel(cif_data, "cif")

            style_func(view)
            view.zoomTo()
            viewer_html = view._make_html()
            components.html(viewer_html, height=500, width=500)
        
  except Exception as e:
    st.error(f"An error occurred: {e}")
    
else:
  st.info("Please enter a PDB ID in the sidebar to get started.")



      
