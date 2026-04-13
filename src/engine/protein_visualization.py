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

    if show_ligands:
        if ligand_colors is not None:
            for ligand_resn, color in ligand_colors.items():
                radius = 0.2 if ligand_resn == "HOH" else 0.4
                view.addStyle({'resn': ligand_resn}, {'stick': {'color': color, 'radius': radius}})
        else:
            # Default highlight for ligands
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})

    if active_site_residues is not None:
        for site in active_site_residues:
            # site expected as {'chain': 'A', 'resnum': 123}
            view.addStyle({'chain': site['chain'], 'resi': str(site['resnum'])},
                          {'stick': {'color': active_site_color, 'radius': 0.3},
                           'sphere': {'color': active_site_color, 'radius': 0.8}})

    view.zoomTo()
    return view
