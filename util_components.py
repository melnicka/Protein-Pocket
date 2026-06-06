import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import requests
from huggingface_hub import InferenceClient

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
from src.utils.fetch_pathways import get_pathways

pocket_radius = 5.0

hf_client = InferenceClient(api_key=st.secrets["HF_TOKEN"])
HF_MODEL_ID = "meta-llama/Meta-Llama-3-8B-Instruct"

def protein_details(pdb_id, entry):
    metadata = entry.extract_metadata()
    try:
        # 1. Fetch PDB Data
        with st.spinner(f"Fetching {pdb_id} from RCSB PDB..."):
            cif_path = fetch_cif(pdb_id)

            pdb_response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
            raw_pdb_string = pdb_response.text if pdb_response.status_code == 200 else None
            custom_metadata = {}
            if raw_pdb_string:
                custom_metadata = extract_metadata(raw_pdb_string)

            fasta_response = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb_id}")
            fasta_sequence = fasta_response.text if fasta_response.status_code == 200 else None

        #  2. Parse PDB Data and find pockets
        with st.spinner("Analyzing structure and finding pockets..."):
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

        # 3. Fetch UniProt Data
        with st.spinner("Fetching UniProt annotations..."):
            accession = get_uniprot_accession(pdb_id)
            prot_name, prot_func = None, None

            if accession:
                prot_name = get_protein_name_from_uniprot(accession)

        # fetch similar proteins
        with st.spinner("Fetching similar proteins..."):
            accession = get_uniprot_accession(pdb_id)
            similar_proteins = []
            if accession:
                similar_proteins = get_similar_proteins(accession, level=50)

        # Left Column: Metadata & Ligands
        st.subheader("Structure information")

        if accession:
            if prot_name:
                st.write(f"**Protein Name:** {prot_name}")

            if custom_metadata:
                st.markdown(f"""
                **Type:** {custom_metadata.get('header', 'N/A')}  \n
                - **Resolution:** {custom_metadata.get('resolution', 'N/A')} Å \n
                - **Experiment:** {custom_metadata.get('experiment', 'N/A')}  \n
                - **Organism:** {custom_metadata.get('organism', 'N/A')}  \n
                """)

            st.success(f"**UniProt Accession:** [{accession}](https://www.uniprot.org/uniprotkb/{accession}/entry)")
        else:
            st.warning("No UniProt linkage found for this structure.")

        if fasta_sequence:
            sequences_dict = {}
            current_header = ""

            for line in fasta_sequence.splitlines():
                if line.startswith(">"):
                    parts = line.split("|")
                    if len(parts) > 1:
                        current_header = parts[1].replace("Chains", "").replace("Chain", "").strip()
                    else:
                        current_header = "Unknown"
                    sequences_dict[current_header] = ""
                else:
                    if current_header:
                        sequences_dict[current_header] += line.strip()

            display_aa_seq(sequences_dict)

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
        prop_col8.metric("Dipole Moment", format_descriptor(dipole_mag, decimals=1, unit="D"))

        if aa_comp:
            with st.expander("Amino Acid Composition (%)"):
                aa_df = pd.DataFrame(list(aa_comp.items()), columns=['Amino Acid', 'Percentage']).set_index('Amino Acid')
                st.bar_chart(aa_df)

        st.divider()

        st.subheader("Ligands & Binding Pockets")
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

        st.divider()

        st.subheader("Metabolic Pathways")
        pathways = get_pathways(pdb_id)
        if pathways:
            for pathway in pathways:
                st.link_button(
                    pathway["name"],
                    pathway["url"]
                )

        else:
            st.info("No KEGG pathway data found.")

        st.divider()

        st.subheader("Similar proteins (UniRef90)")
        if similar_proteins:
            for entry in similar_proteins:
                accession_id = entry["primaryAccession"]
                name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                organism = entry["organism"]["scientificName"]

                id_col, name_col, org_col, btn_col = st.columns(4)
                id_col.write(f"{accession_id}")
                name_col.write(f"{name}")
                org_col.write(f"{organism}")
                with btn_col:
                    st.page_link("pages/comparison_page.py", label="Compare", query_params={"acc1": f"{accession}", "acc2": f"{accession_id}"})
        else:
            st.info("No similar proteins with PDB entries found")

        # ✅ FIXED: Explicitly executing the LLM interface block inside the dashboard workflow
        st.divider()
        st.subheader("🤖 AI Structure Assistant")
        
        safe_metadata = custom_metadata if custom_metadata else {}
        safe_name = prot_name if prot_name else "Unknown Protein"
        
        try:
            llm_implementation(
                pdb_id=pdb_id,
                entry=entry,
                custom_metadata=safe_metadata,
                prot_name=safe_name,
                global_sasa=global_sasa,
                rg=rg,
                pi=pi,
                instability=instability,
                helix_frac=helix_frac
            )
        except Exception as llm_err:
            st.error(f"Could not load protein details to feed to the AI: {llm_err}")

        return {
            "prot_name": safe_name,
            "custom_metadata": safe_metadata, 
            "global_sasa": global_sasa,
            "rg": rg,
            "pi": pi,
            "instability": instability,
            "helix_frac": helix_frac
        }

    except Exception as e:
        st.error(f"An error occurred: {e}")


def structure_viewer(entry, cif_path, bg_color, selected_chains, chain_colors, style_option, color_option):
    view_options = ["Whole Protein"]
    pocket_ligands = [lig for lig in getattr(entry, 'ligands', []) if hasattr(lig, 'pocket') and not lig.pocket.is_empty]
    for lig in pocket_ligands:
        view_options.append(f"Pocket: {lig.comp_id} ({lig.auth_asym_id}:{lig.auth_seq_id})")

    selected_view = st.selectbox("Select structure to render:", view_options, key=entry.pdb_id)

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
        is_pocket_view=is_pocket,
        style=style_option,
        color_scheme=color_option
    )

    components.html(view._make_html(), height=600, width=600)

def llm_implementation(pdb_id, entry, custom_metadata, prot_name, global_sasa, rg, pi, instability, helix_frac):

    try:
        avg_b_factor = np.mean(entry.atom_array.b_factor)
    except Exception:
        avg_b_factor = None

    protein_summary = f"""
    Protein: {pdb_id}
    Name: {prot_name}
    Resolution: {custom_metadata.get('resolution', 'N/A')} Å
    SASA: {f"{global_sasa:.2f} Å²" if global_sasa is not None else 'N/A'}
    Gyration Radius: {f"{rg:.2f} Å" if rg is not None else 'N/A'}
    Instability Index: {f"{instability:.1f} ({'Stable' if instability < 40 else 'Unstable'})" if instability is not None else 'N/A'}
    Helix Content: {f"{helix_frac:.1f}%" if helix_frac is not None else 'N/A'}
    Isoelectric Point: {f"{pi:.2f}" if pi is not None else 'N/A'}
    B-factor avg: {f"{avg_b_factor:.2f}" if avg_b_factor is not None else 'N/A'}
    Ligands found: {len(entry.ligands) if hasattr(entry, 'ligands') else 0}
    """
    
    chat_key = f"hf_chat_{pdb_id}"
    
    if chat_key not in st.session_state:
        system_instruction = f"""
        You are an expert computational biologist. You are analyzing the following protein data:
        {protein_summary}
        
        Answer all user questions based strictly on this data and standard biochemical principles. 
        Keep your answers concise, scientifically accurate, and format them clearly with markdown.
        """
        st.session_state[chat_key] = [
            {"role": "system", "content": system_instruction}
        ]

    chat_history = st.session_state[chat_key]
    
    with st.container():
        for message in chat_history:
            if message["role"] == "system":
                continue

            if message["role"] == "user" and "Analyze this protein structure data" in message["content"]:
                continue 
                
            with st.chat_message(message["role"]):
                st.markdown(message["content"])

    if len(chat_history) == 1:

        if st.button("✨ Generate Initial Structure Insights", key=f"btn_generate_{pdb_id}"):
            with st.spinner("AI is analyzing your protein..."):
                try:
                    prompt = """Analyze this protein structure data and provide:
                    1. 2-3 key structural features
                    2. Implications for protein function/stability
                    3. Notable characteristics compared to typical proteins
                    """
                    chat_history.append({"role": "user", "content": prompt})
                    
                    response = hf_client.chat.completions.create(
                        model=HF_MODEL_ID,
                        messages=chat_history,
                        max_tokens=500
                    )
                    reply = response.choices[0].message.content
                    
                    chat_history.append({"role": "assistant", "content": reply})
                    st.rerun()
                except Exception as e:
                    st.error(f"Failed to generate insights: {e}")

    if prompt := st.chat_input(f"Ask a follow-up question about {pdb_id}...", key=f"chat_input_{pdb_id}"):
        with st.chat_message("user"):
            st.markdown(prompt)
            
        chat_history.append({"role": "user", "content": prompt})
            
        with st.chat_message("assistant"):
            with st.spinner("Thinking..."):
                try:
                    response = hf_client.chat.completions.create(
                        model=HF_MODEL_ID,
                        messages=chat_history,
                        max_tokens=500
                    )
                    reply = response.choices[0].message.content
                    
                    st.markdown(reply)
                    chat_history.append({"role": "assistant", "content": reply})
                except Exception as e:
                    st.error(f"Error communicating with Hugging Face API: {e}")

    if len(chat_history) > 1:
        st.divider()
        chat_text = f"--- Structure Analysis for {pdb_id} ---\n\n"
        for message in chat_history:
            if message["role"] == "system":
                continue
                
            role = "Assistant" if message["role"] == "assistant" else "User"
            if role == "User" and "Analyze this protein structure data" in message["content"]:
                continue
                
            chat_text += f"**{role}:**\n{message['content']}\n\n"
        
        st.download_button(
            label="📥 Download Chat History",
            data=chat_text,
            file_name=f"{pdb_id}_analysis.txt",
            mime="text/plain",
            key=f"dl_button_{pdb_id}" 
        )

def show_plddt_legend():
    legend_items = [
        ("#0053D6", "Very high (pLDDT > 90)"),
        ("#65CBF3", "Confident (90 > pLDDT > 70)"),
        ("#FFDB13", "Low (70 > pLDDT > 50)"),
        ("#FF7045", "Very low (pLDDT < 50)")
    ]

    rows_html = ""
    for color, text in legend_items:
        rows_html += f'<div style="display: flex; align-items: center; gap: 12px;"><div style="width: 22px; height: 22px; background-color: {color}; border-radius: 4px;"></div><span style="color: #d8dee9; font-weight: 500; font-size: 14px;">{text}</span></div>'
        
    final_html = f'<div style="display: flex; flex-direction: column; gap: 8px; background-color: #11151c; padding: 14px; border-radius: 8px; border: 1px solid #2e3440;">{rows_html}</div>'
    
    st.markdown(final_html, unsafe_allow_html=True)
