import biotite.structure as struct 
import biotite.sequence as seq
import biotite.sequence.seqstats as seqstats
import Bio.SeqUtils.ProtParam import ProteinAnalysis
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

'''
Calculate the percentage occurrence of each amino acid in the protein sequence based on symbol frequency
'''
def calc_amino_acid_composition(
        atom_array: struct.AtomArray
    ) -> dict:
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    freq = seqstats.get_symbol_frequency(residues)
    percent = freq / len(protein_seq)*100
    
    return percent

'''
Calculate the instability index according to Guruprasad et al 1990.
Any value above 40 means the protein is unstable (has a short half life).
'''
def calc_instability_index(
        atom_array: struct.AtomArray
    ) -> float:
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
def calc_aromaticity(atom_array: struct.AtomArray) -> float:
    residues = struct.get_residues(atom_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    analysis = ProteinAnalysis("".join(protein_seq))
    return analysis.aromaticity()
