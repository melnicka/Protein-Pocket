import streamlit as st
import streamlit.components.v1 as components
import py3Dmol

from src.utils.data_fetching import fetch_cif
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




      
