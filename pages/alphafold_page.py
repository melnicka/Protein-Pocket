import streamlit as st
import pandas as pd
import numpy as np
import requests
import streamlit.components.v1 as components
from util_components import pocket_radius, show_plddt_legend
from src.engine.entry import Entry
from src.utils.fetch_alphafold import fetch_and_clean_alphafold
from src.utils.fetch_uniprot import get_protein_name_from_uniprot, get_papers_from_accession
from src.engine.protein_visualization import visualize_structure
from src.engine.descriptors import (
    format_descriptor, calc_sasa_protein, calc_gyration_radius,
    calc_instability_index, calc_aromaticity, calc_helix_fraction,
    calc_isoelectric_point, calc_amino_acid_composition, calc_dipole_moment
)
from src.utils.fetch_pathways import get_pathways

st.set_page_config(page_title="AlphaFold Explorer", layout="wide")
st.title("AlphaFold Explorer")
st.markdown("Analyze AI-predicted protein structures from the AlphaFold Database.")

uniprot_id = st.text_input("Enter UniProt ID:", "Q9H6R3").strip().upper()

if uniprot_id:
    with st.spinner(f"Processing {uniprot_id}..."):
        try:
            cif_path = fetch_and_clean_alphafold(uniprot_id)

            entry = Entry(cif_path, uniprot_id)
            entry.find_pockets(search_radius=pocket_radius, filter_out_solvent=True)
            metadata = entry.extract_metadata()

            col1, col2 = st.columns([1, 1])

            with col1:
                st.markdown("### Structure information")
                protein_name = get_protein_name_from_uniprot(uniprot_id)

                organism_name = "UNKNOWN"
                try:
                    res = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
                    if res.ok:
                        organism_name = res.json().get("organism", {}).get("scientificName", "UNKNOWN").upper()
                except Exception:
                    pass
                
                st.markdown(f"**Protein Name:** {protein_name}")
                st.markdown("**Type:** Predicted Model (AlphaFold)")
                st.markdown(f"- **Organism:** {organism_name}")
                
                pathways = []
                try:
                    pathways = get_pathways(uniprot_id, is_uniprot=True)
                except Exception:
                    pass
                
                if pathways:
                    st.markdown("**Metabolic Pathways:**")
                    for p in pathways:
                        st.markdown(f"[{p['name']}]({p['url']})")

                st.write("")

                with st.expander("📄 View FASTA sequence"):
                    try:
                        fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
                        fasta_res = requests.get(fasta_url)
                        if fasta_res.status_code == 200:
                            st.code(fasta_res.text, language="fasta")
                        else:
                            st.warning("FASTA sequence not found in UniProt.")
                    except Exception:
                        st.error("Failed to connect to UniProt API.")
                
                st.markdown("### Model Confidence Legend")
                show_plddt_legend()
                st.markdown("---")

                st.markdown("### Whole Protein Properties:")
                arr = entry.atom_array
                
                sasa = calc_sasa_protein(arr)
                rg = calc_gyration_radius(arr)
                pi = calc_isoelectric_point(arr)
                arom = calc_aromaticity(arr)
                instab = calc_instability_index(arr)
                helix = calc_helix_fraction(arr)
                total_atoms = len(arr)
                
                dipole_vec = calc_dipole_moment(arr)
                dipole_mag = np.linalg.norm(dipole_vec) 
                
                m_col1, m_col2, m_col3, m_col4 = st.columns(4)
                m_col1.metric("Total SASA", format_descriptor(sasa, "Å²", 0))
                m_col2.metric("Gyration Rad.", format_descriptor(rg, "Å", 2))
                m_col3.metric("pI", format_descriptor(pi, "", 2))
                m_col4.metric("Aromaticity", format_descriptor(arom, "", 1, is_percent=True))

                m_col5, m_col6, m_col7, m_col8 = st.columns(4)
                m_col5.metric("Instability", format_descriptor(instab, "", 1))
                m_col6.metric("Helix Frac.", format_descriptor(helix, "", 1, is_percent=True))
                m_col7.metric("Total Atoms", str(total_atoms))
                m_col8.metric("Dipole Moment", format_descriptor(dipole_mag, "D", 1))

                st.write("") 

                with st.expander("› Amino Acid Composition (%)"):
                    aa_comp = calc_amino_acid_composition(arr)
                    aa_df = pd.DataFrame.from_dict(aa_comp, orient='index', columns=['Percentage'])
                    st.bar_chart(aa_df)

            st.sidebar.header("Visualization Settings")
            available_chains = list(metadata.get('chains', []))
            selected_chains = st.sidebar.multiselect("Visible chain(s)", available_chains, default=available_chains, key="af_chains")
            
            with st.sidebar:
                style_option = st.radio("Visualization Style", ["cartoon", "stick", "surface", "sphere", "line"], key="af_style")
                bg_color = st.selectbox("Background", options=["white", "black", "lightgray"], index=0, key="af_bg")

            with col2:
                st.subheader("3D Structure Viewer")
                with open(cif_path, "r", encoding="utf-8") as f:
                    cif_data = f.read()
                
                view = visualize_structure(
                    cif_data=cif_data,
                    style=style_option,
                    color_scheme="spectrum",
                    bg_color=bg_color,
                    selected_chains=selected_chains,
                    chain_colors=None,
                    is_pocket_view=False 
                )
                
                plddt_style = {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 100}}
                
                view.setStyle({}, {})
                if style_option == "surface":
                    view.addSurface("VDW", plddt_style, {"chain": selected_chains} if selected_chains else {})
                else:
                    if selected_chains:
                        for chain in selected_chains:
                            view.setStyle({"chain": chain}, {style_option: plddt_style})
                    else:
                        view.setStyle({}, {style_option: plddt_style})
                
                view.zoomTo()

                components.html(view._make_html(), height=600)

                papers = []
                try:
                    papers = get_papers_from_accession(uniprot_id)
                except Exception as e:
                    st.warning(f"Could not fetch PubMed articles: {e}")

                if papers:
                    st.subheader("PubMed Papers")
                    with st.container(height=450):
                        for p in papers:
                            st.markdown(f"### {p['title']}")
                            st.markdown(f"**Summary:** {p['summary']}")
                            st.link_button("Open PubMed", p["url"])
                            st.divider()
                else:
                    st.info("No PubMed articles found for this protein.")

        except Exception as e:
            st.error(f"An error occurred: {e}")
else:
    st.info("Please enter a UniProt ID to start analysis.")