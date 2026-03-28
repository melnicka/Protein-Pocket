import streamlit as st
import streamlit.components.v1 as components
import py3Dmol

from src.utils.data_fetching import fetch_cif
from src.utils.cif_parsing import load_cif, extract_entry_metadata, extract_ligands
from src.utils.fetch_uniprot import get_uniprot_accession, get_protein_name_from_uniprot, get_function_from_uniprot
from src.engine.entry import Entry

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

      st.subheader("Ligands & Binding Pockets")
      if len(entry.ligands) >0:
        for lig in entry.ligands:
          with st.expander(f"Ligand {lig.comp_id} (Chain , Seq)"):
            st.write(f'**SMILES** `{lig.smiles}`')
            
            if not lig.pocket.is_empty:
              st.success(f"Pocket found within {pocket_radius}Å")

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



      
