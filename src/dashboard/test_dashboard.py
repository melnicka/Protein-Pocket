import streamlit as st
from deep_translator import GoogleTranslator

st.set_page_config(page_title="3D Protein Viewer", layout="wide")

col_a, col_b = st.columns([8, 1])
with col_b:
    translate = st.checkbox("🇵🇱", value=True, help="Tłumacz na polski")
st.title("🧬 Wizualizacja struktury białka 3D" if translate else "🧬 3D Protein Structure Visualization")

def ligand_color_controls(ligand_counts):
    st.subheader("🧪 Kolory Ligandów:" if translate else "🧪 Ligand Coloring:")

    color_options = ["green", "red", "blue", "yellow", "cyan", "magenta", "orange", "white"]
    color_names_pl = {
        "green": "zielony",
        "red": "czerwony",
        "blue": "niebieski",
        "yellow": "żółty",
        "cyan": "błękitny",
        "magenta": "magenta",
        "orange": "pomarańczowy",
        "white": "biały"
    }

    ligand_colors = {}
    ligands = list(ligand_counts.items())

    num_cols = 2 if len(ligands) > 1 else 1
    cols = st.columns(num_cols)

    for idx, (ligand, count) in enumerate(ligands):
        col = cols[idx % num_cols]
        with col:
            ligand_colors[ligand] = st.selectbox(
                f"{ligand} ({count})",
                options=color_options,
                index=idx % len(color_options),
                key=f"ligand_color_{ligand}",
                format_func=lambda x: color_names_pl.get(x, x) if translate else x
            )
    return ligand_colors

def display_amino_acid_sequences(sequences, translate=True):
    if translate:
        st.subheader("🧬 Sekwencje aminokwasowe")
        expander_label = "📄 Zobacz sekwencje w formacie FASTA"
        download_label = "📥 Pobierz FASTA"
        chain_label = "Łańcuch"
        length_label = "Długość"
    else:
        st.subheader("🧬 Amino Acid Sequences")
        expander_label = "📄 View FASTA sequences"
        download_label = "📥 Download FASTA"
        chain_label = "Chain"
        length_label = "Length"

    fasta_output = ""
    for chain, seq in sequences.items():
        st.markdown(f"**{chain_label} {chain}** ({length_label}: {len(seq)})")
        chunk_size = 60
        for i in range(0, len(seq), chunk_size):
            color = "#ffffff" if st.session_state.get("dark_mode", False) else "#000000"
            st.markdown(f"<pre style='color:{color}'>{i + 1:>4} {seq[i:i + chunk_size]}</pre>", unsafe_allow_html=True)     
        st.markdown("---")

        # Tworzenie FASTA
        fasta_output += f">Chain_{chain}\n"
        for i in range(0, len(seq), chunk_size):
            fasta_output += seq[i:i + chunk_size] + "\n"

    # Wyświetlanie FASTA
    with st.expander(expander_label):
        st.code(fasta_output.strip(), language="fasta")

    # Przycisk do pobrania
    st.download_button(
        label=download_label,
        data=fasta_output,
        file_name="protein_sequences.fasta",
        mime="text/plain"
    )

def main():
    pdb_id = st.text_input(
        "Wpisz PDB ID (np. 3PAV, 6LU7)" if translate else "Enter PDB ID (e.g., 3PAV, 6LU7):",
        "3PAV"
    ).strip().upper()

    if not pdb_id:
        return

    with st.spinner("Ładuję..." if translate else "Loading..."):
        pdb_data = fetch_pdb(pdb_id)

    if not pdb_data:
        return

    metadata = extract_metadata(pdb_data)
    if translate:
        try:
            for key in ['resolution', 'experiment', 'organism']:
                if metadata.get(key):
                    metadata[key] = GoogleTranslator(source='en', target='pl').translate(metadata[key])
        except Exception as e:
            st.warning(f"Błąd tłumaczenia metadanych: {e}")

accession = get_uniprot_accession(pdb_id)
    name = get_protein_name_from_uniprot(accession, translate=translate)
    function = get_function_from_uniprot(accession, translate=translate)

with st.sidebar:
    
        st.header("⚙️ Ustawienia" if translate else "⚙️ Settings")
        style = st.radio("Styl wizualizacji" if translate else "Visualization style", ["cartoon", "stick", "surface", "sphere", "line"])
        color = st.radio("Schemat kolorów" if translate else "Color scheme", ["spectrum", "chain", "residue"])
        ligands = st.checkbox("Pokaż ligandy" if translate else "Show ligands", value=True)

        chains = sorted(metadata['chains'])
        selected_chains = st.multiselect(
            "Wybierz łańcuch(y)" if translate else "Select chain(s)",
            chains,
            default=chains
        )

        st.markdown("### 🎨 Kolory dla łańcuchów" if translate else "### 🎨 Chain colors")
        chain_colors = {}
        for chain in chains:
            picked = st.color_picker(f"{chain}", value="#FFFFFF")
            chain_colors[chain] = None if picked.upper() == "#FFFFFF" else picked

         bg_options = ["white", "black", "lightgray"]
        bg_names_pl = {
            "white": "biały",
            "black": "czarny",
            "lightgray": "jasnoszary"
        }

        bg_color = st.selectbox(
            "Tło" if translate else "Background",
            options=bg_options,
            format_func=lambda x: bg_names_pl.get(x, x) if translate else x
        )


    ligand_colors = ligand_color_controls(metadata['ligand_counts']) if metadata['ligand_counts'] else None

    # Wyświetlanie
    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("🔍 Informacje o strukturze" if translate else "🔍 Structure Information")
        if translate:
            st.markdown(f"""
            - **Tytuł:** {metadata['title']}
            - **Rozdzielczość:** {metadata['resolution']}
            - **Organizm:** {metadata['organism']}
            - **Eksperyment:** {metadata['experiment']}
            - **Łańcuchy:** {', '.join(sorted(metadata['chains']))}
            - **Reszty:** {metadata['residues']}
            - **Atomy:** {metadata['atoms']}
            - **Mostki disulfidowe:** {len(metadata['disulfide_bonds'])}
            - **Miejsca aktywne:** {len(metadata['active_sites'])}
            """)
        else:
            st.markdown(f"""
            - **Title:** {metadata['title']}
            - **Resolution:** {metadata['resolution']}
            - **Organism:** {metadata['organism']}
            - **Experiment:** {metadata['experiment']}
            - **Chains:** {', '.join(sorted(metadata['chains']))}
            - **Residues:** {metadata['residues']}
            - **Atoms:** {metadata['atoms']}
            - **Disulfide Bonds:** {len(metadata['disulfide_bonds'])}
            - **Active Sites:** {len(metadata['active_sites'])}
            """)

        st.subheader("🔬 Informacje o białku" if translate else "🔬 Protein Information")
        if name:
            st.markdown(f"**Nazwa:** {name}" if translate else f"**Name:** {name}")

  with col2:
        if translate:
                    st.markdown("""
                    **Przewodnik interakcji:**
                    - Najedź kursorem na atomy, aby zobaczyć informacje o resztach
                    - Kliknij na atomy, aby przypiąć informacje
                    - Kliknij ponownie, aby usunąć przypiętą etykietę
                    """)

        else:
                    st.markdown("""
                    **Interaction Guide:**
                    - Hover over atoms to see residue information  
                    - Click on atoms to pin information  
                    - Click again to remove pinned label
                    """)

        viewer = visualize_protein(
            pdb_data,
            style=style,
            chain_colors=chain_colors,
            show_ss=True,
            show_ligands=ligands,
            ligand_colors=ligand_colors,
            show_active_sites=True,
            active_site_residues=metadata['active_sites'],
            selected_chains=selected_chains,
            bg_color=bg_color
        )


  display_amino_acid_sequences(metadata['sequences'], translate=translate)

    if translate:
        st.markdown(f"[🔗 Zobacz strukturę w bazie RCSB PDB](https://www.rcsb.org/structure/{pdb_id})")
        st.markdown("[🔗 Odwiedź bazę RCSB PDB](https://www.rcsb.org)")
        st.markdown(f"[🔗 Zobacz na UniProt](https://www.uniprot.org/uniprotkb/{accession})")
    else:
        st.markdown(f"[🔗 View this structure on RCSB PDB](https://www.rcsb.org/structure/{pdb_id})")
        st.markdown("[🔗 Go to RCSB PDB](https://www.rcsb.org)")
        st.markdown(f"[🔗 View on UniProt](https://www.uniprot.org/uniprotkb/{accession})")

if __name__ == "__main__":
    main()




