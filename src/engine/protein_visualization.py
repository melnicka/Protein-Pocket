import py3Dmol
import streamlit as st

def visualize_structure(cif_data, style="cartoon", bg_color="white", 
                        show_ss=True, chain_colors=None, selected_chains=None,
                        show_ligands=True, ligand_colors=None,
                        active_site_residues=None, active_site_color="magenta",
                        is_pocket_view=False):
    
    view = py3Dmol.view(width=600, height=600)
    view.addModel(cif_data, "cif")
    view.setBackgroundColor(bg_color)

    # 1. Handle Chain Coloring / Selection
    if not is_pocket_view:
        if selected_chains is not None and chain_colors is not None:
            # Color selected chains; fade out the rest
            all_chains = list(chain_colors.keys()) if chain_colors else []
            for chain in all_chains:
                if chain in selected_chains:
                    color = chain_colors.get(chain, "spectrum")
                    view.setStyle({'chain': chain}, {style: {'color': color}})
                else:
                    view.setStyle({'chain': chain}, {'cartoon': {'color': 'gray', 'opacity': 0.2}})
        else:
            # Default spectrum view
            view.setStyle({style: {'color': 'spectrum'}})
    else:
        # POCKET VIEW: Stick representation for the local environment
        view.setStyle({'stick': {'colorscheme': 'grayCarbon'}})

    # 2. Secondary Structure Highlights (Helices = Red, Sheets = Yellow)
    if show_ss is True and not is_pocket_view:
        view.addStyle({'helix': True}, {'cartoon': {'color': 'red'}})
        view.addStyle({'sheet': True}, {'cartoon': {'color': 'yellow'}})

    # 3. Ligand Styling
    if show_ligands:
        if ligand_colors is not None:
            for ligand_resn, color in ligand_colors.items():
                radius = 0.2 if ligand_resn == "HOH" else 0.4
                view.addStyle({'resn': ligand_resn}, {'stick': {'color': color, 'radius': radius}})
        else:
            # Default highlight for ligands
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})

    # 4. Active Site / Pocket Highlight
    if active_site_residues is not None:
        for site in active_site_residues:
            # site expected as {'chain': 'A', 'resnum': 123}
            view.addStyle({'chain': site['chain'], 'resi': str(site['resnum'])},
                          {'stick': {'color': active_site_color, 'radius': 0.3},
                           'sphere': {'color': active_site_color, 'radius': 0.8}})

    view.zoomTo()
    return view

def ligand_color_controls(unique_ligands, translate=False):
    """Generates sidebar selectboxes for choosing ligand colors."""
    
    color_options = ["green", "red", "blue", "yellow", "cyan", "magenta", "orange", "white", "gray"]
    ligand_colors = {}
    
    if len(unique_ligands) > 0:
        st.sidebar.subheader("Ligand Colors" if not translate else "Kolory Ligandów")
        
        for i, lig_id in enumerate(unique_ligands):
            # Give each ligand a different default color to start
            default_color = color_options[i % len(color_options)]
        

            label = f"Color for {lig_id}" if not translate else f"Kolor dla {lig_id}"
            
            selected_color = st.sidebar.selectbox(
                label,
                options=color_options,
                index=color_options.index(default_color),
                key=f"color_{lig_id}" # Unique key required by Streamlit
            )
            ligand_colors[lig_id] = selected_color

    return ligand_colors
