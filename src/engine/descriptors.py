import biotite.structure as struct 
import biotite.sequence as seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
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
        protein_mask = struct.filter_amino_acids(structure) #removing ligands, water molecules, etc.
        protein_array = structure[protein_mask]
        
        sasa_protein = struct.sasa(protein_array)
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
        gyration_radius = struct.gyration_radius(protein_array, masses=None)
        
        return gyration_radius

def calc_amino_acid_composition(
        atom_array: struct.AtomArray
    ) -> dict:
    '''
    Calculate the percentage occurrence of each amino acid in the protein sequence based on symbol frequency
    '''
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    symbols, counts=np.unique(list(protein_seq), return_counts=True)
    freq = dict(zip(symbols,counts))
    percent = {str(aa): float(count/len(protein_seq)) * 100 for aa, count in freq.items()}
    
    return percent

def calc_instability_index(
        atom_array: struct.AtomArray
    ) -> float:
    '''
    Calculate the instability index according to Guruprasad et al 1990.
    Any value above 40 means the protein is unstable (has a short half life).
    '''
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)

    analysis = ProteinAnalysis("".join(protein_seq))

    return analysis.instability_index()


'''opisuje jaki ułamek aminokwasów w białku jest aromatyczny.'''
def calc_aromaticity(atom_array: struct.AtomArray) -> float:
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    analysis = ProteinAnalysis("".join(protein_seq))
    return analysis.aromaticity()


'''mierzy udział helis w strukturze drugorzędowej białka.'''
def calc_helix_fraction(atom_array: struct.AtomArray) -> float:
    protein_atoms = atom_array[struct.filter_amino_acids(atom_array)]
    sse_annotation = struct.annotate_sse(protein_atoms)
    helix_fraction = np.sum(sse_annotation == "H") / len(sse_annotation)
    return helix_fraction


'''oblicza pH przy którym białko ma ładunek 0'''        
def calc_isoelectric_point(atom_array: struct.AtomArray) -> float:
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    analysis = ProteinAnalysis("".join(protein_seq))
    return analysis.isoelectric_point()

def calc_pocket_centroid(pocket_array: struct.AtomArray) -> tuple:
    """Calculates the geometric center (centroid) of the pocket."""
    if len(pocket_array) == 0:
        return 0.0, 0.0, 0.0
    centroid_coords = struct.centroid(pocket_array)
    
    return float(centroid_coords[0]), float(centroid_coords[1]), float(centroid_coords[2])

def calc_pocket_hydrophobicity(pocket_array: struct.AtomArray) -> float:
    """Calculates the sum of hydrophobicity of the binding pocket 
    using the Kyte-Doolittle scale."""
    if len(pocket_array) == 0:
        return 0.0

    kd_scale = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
        'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
        'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9,
        'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
        'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }

    protein_mask = struct.filter_amino_acids(pocket_array)
    protein_pocket = pocket_array[protein_mask]

    unique_residues = set(zip(protein_pocket.chain_id, protein_pocket.res_id, protein_pocket.res_name))

    total_hydrophobicity = sum(kd_scale.get(res_name, 0.0) for _, _, res_name in unique_residues)

    return float(total_hydrophobicity)
