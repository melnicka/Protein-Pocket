import biotite.structure as struct 
import numpy as np

def calc_ligand_buried_surface(
        pocket_array: struct.AtomArray,
        ligand_array: struct.AtomArray
) -> float:
    """Calculates ligand buried surface.

    Args:
        pocket_array: Structure of the pocket, ligand included.
        ligand_array: Structure of the ligand alone.

    Returns:
        float: Ligand buried surface.
    """
    ligand_sasa = struct.sasa(ligand_array).sum()
    ligand_in_pocket_mask = pocket_array.res_name == ligand_array.res_name[0] 
    pocket_sasa = struct.sasa(pocket_array)
    ligand_in_pocket_sasa = pocket_sasa[ligand_in_pocket_mask].sum()

    return ligand_sasa - ligand_in_pocket_sasa
        
"""
Solvent Accessible Surface Area (SASA):
Measures how much of the protein's surface area is physically 
exposed to the surrounding solvent (usually water).
"""
def calc_sasa_protein(structure: struct.AtomArray) -> float:
        protein_mask = struc.filter_amino_acids(structure) #removing ligands, water molecules, etc.
        protein_array = structure[protein_mask]
        
        sasa_protein = struc.sasa(protein_array)
        total_sasa = np.nansum(sasa_protein) 
        
        return total_sasa
        
"""
Radius of Gyration (Rg):
Measures the overall size and compactness of the protein structure.
It calculates the mass-weighted average distance of all atoms from 
the protein's center of mass.
- Low Rg: The protein is tightly folded and compact
- High Rg: The protein is extended, unfolded, or highly flexible
"""
def calc_gyration_radius(protein_array: struct.AtomArray) -> float:
        gyration_radius = struct.gyration_radius(protein_array, masses=none)
        
        return gyration_radius


