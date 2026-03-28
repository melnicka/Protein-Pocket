from src.engine.entry import Entry
from src.utils.data_fetching import fetch_cif
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

