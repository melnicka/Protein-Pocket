import io
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile, set_structure  # <-- Explicit imports
from rdkit import Chem
from rdkit.Chem import Lipinski
import biotite.sequence as seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

def format_descriptor(value, unit="", decimals=2, is_percent=False):
    """Formats values, adding units and handling missing data (None)."""
    if value is None:
        return "N/A"
    
    try:
        val = float(value)
        if is_percent:
            return f"{(val * 100):.{decimals}f}%"
        if unit:
            return f"{val:.{decimals}f} {unit}"
        
        return f"{val:.{decimals}f}"
    except (ValueError, TypeError):
        return str(value)

def _get_protein_analysis(atom_array: struc.AtomArray) -> ProteinAnalysis:
    """Helper function to extract sequence and initialize Biopython's ProteinAnalysis."""
    protein_array = atom_array[struc.filter_amino_acids(atom_array)]  #<- dodalam filtrowanie, by ignorowalo 'smieci'~NB
    residues = struc.get_residues(protein_array)[1]
    protein_seq = seq.ProteinSequence(residues)
    return ProteinAnalysis(str(protein_seq))

        
def calc_ligand_buried_surface(
        pocket_array: struc.AtomArray,
        ligand_array: struc.AtomArray
) -> float:
    """Calculates ligand buried surface.

    Args:
        pocket_array: Structure of the pocket, ligand included.
        ligand_array: Structure of the ligand alone.

    Returns:
        float: Ligand buried surface.
    """
    ligand_sasa = struc.sasa(ligand_array).sum()
    ligand_in_pocket_mask = pocket_array.res_name == ligand_array.res_name[0] 
    pocket_sasa = struc.sasa(pocket_array)
    ligand_in_pocket_sasa = pocket_sasa[ligand_in_pocket_mask].sum()

    return ligand_sasa - ligand_in_pocket_sasa
        

def calc_sasa_protein(structure: struc.AtomArray) -> float:
    """Solvent Accessible Surface Area (SASA):
Measures how much of the protein's surface area is physically 
exposed to the surrounding solvent (usually water).
"""
    protein_mask = struc.filter_amino_acids(structure) #removing ligands, water molecules, etc.
    protein_array = structure[protein_mask]
        
    sasa_protein = struc.sasa(protein_array, ignore_ions=True) #<- na wszelki wypadek dopisalam ignore~NB
    return float(np.nansum(sasa_protein))
        
def calc_gyration_radius(atom_array: struc.AtomArray) -> float:
    """
Radius of Gyration (Rg):
Measures the overall size and compactness of the protein structure.
It calculates the mass-weighted average distance of all atoms from 
the protein's center of mass.
- Low Rg: The protein is tightly folded and compact
- High Rg: The protein is extended, unfolded, or highly flexible
"""
    protein_array = atom_array[struc.filter_amino_acids(atom_array)]
    return struc.gyration_radius(protein_array, masses=None)


def calc_amino_acid_composition(atom_array: struc.AtomArray) -> dict:
    '''Calculate the percentage occurrence of each amino acid in the protein sequence based on symbol frequency'''
    analysis = _get_protein_analysis(atom_array)
    return {aa: percent * 100 for aa, percent in analysis.get_amino_acids_percent().items()}

def calc_instability_index(atom_array: struc.AtomArray) -> float:
    '''
    Calculate the instability index according to Guruprasad et al 1990.
    Any value above 40 means the protein is unstable (has a short half life).
    '''
    analysis = _get_protein_analysis(atom_array)
    return analysis.instability_index()

def calc_aromaticity(atom_array: struc.AtomArray) -> float:
    """Calculates the fraction of aromatic amino acids in the protein."""
    analysis = _get_protein_analysis(atom_array)
    return analysis.aromaticity()

def calc_helix_fraction(atom_array: struc.AtomArray) -> float:
    """Measures the fraction of helices in the secondary structure."""
    protein_atoms = atom_array[struc.filter_amino_acids(atom_array)]
    sse_annotation = struc.annotate_sse(protein_atoms)
    
    if len(sse_annotation) == 0:
        return 0.0
        
    helix_fraction = np.sum(sse_annotation == "a") / len(sse_annotation)   #tu chyba mialo byc a zamiast H?~NB
    return float(helix_fraction)
       
def calc_isoelectric_point(atom_array: struc.AtomArray) -> float:
    '''oblicza pH przy którym białko ma ładunek 0''' 
    analysis = _get_protein_analysis(atom_array)
    return analysis.isoelectric_point()

def calc_pocket_centroid(pocket_array: struc.AtomArray) -> tuple:
    """Calculates the geometric center (centroid) of the pocket."""
    if len(pocket_array) == 0:
        return 0.0, 0.0, 0.0
    centroid_coords = struc.centroid(pocket_array)
    
    return float(centroid_coords[0]), float(centroid_coords[1]), float(centroid_coords[2])

def calc_pocket_hydrophobicity(pocket_array: struc.AtomArray) -> float:
    """Calculates the sum of hydrophobicity of the binding pocket using the Kyte-Doolittle scale."""
    if len(pocket_array) == 0:
        return 0.0

    kd_scale = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
        'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
        'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9,
        'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
        'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }

    protein_mask = struc.filter_amino_acids(pocket_array)
    protein_pocket = pocket_array[protein_mask]

    unique_residues = set(zip(protein_pocket.chain_id, protein_pocket.res_id, protein_pocket.res_name))

    total_hydrophobicity = sum(kd_scale.get(res_name, 0.0) for _, _, res_name in unique_residues)

    return float(total_hydrophobicity)

def calc_dipole_moment(atom_array: struc.AtomArray):

    # ladunki czastkowe, niektore atomy nie maja nw jakies bledy ...
    charges = struc.partial_charges(atom_array)

    centroid = np.sum([atom.coord for atom in atom_array]) / len(atom_array)
    dipole_moment = centroid
    for atom, charge in zip(atom_array, charges):
        partial_moment = atom.coord * charge
        dipole_moment += np.nan_to_num(partial_moment)

    return dipole_moment

def calc_hydrogen_bond_features(pocket_array: struc.AtomArray) -> dict:
    if len(pocket_array) == 0:
        return {"HBA": 0, "HBD": 0}
    
    pdb_file = PDBFile()
    set_structure(pdb_file, pocket_array)

    with io.StringIO() as buffer:
        pdb_file.write(buffer)
        pdb_block = buffer.getvalue()

    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True)
    if mol is None:
        return {"HBA": 0, "HBD": 0}

    HBD = Lipinski.NumHDonors(mol)
    HBA = Lipinski.NumHAcceptors(mol)

    return {"HBD": int(HBD), "HBA": int(HBA)}


def calc_charged_surface_fraction(pocket_array: struc.AtomArray) -> float:
    if len(pocket_array) == 0:
        return 0.0
    
    charged_residues = {"ARG", "LYS", "HIS", "ASP", "GLU"}
    sasa = struc.sasa(pocket_array)

    amino_acid_mask = struc.filter_amino_acids(pocket_array)
    amino_acid_atoms = pocket_array[amino_acid_mask]
    amino_acid_sasa = sasa[amino_acid_mask]

    if len(amino_acid_atoms) == 0:
         return 0.0
    
    total_sasa = np.nansum(amino_acid_sasa)

    charged_mask = np.isin(amino_acid_atoms.res_name, list(charged_residues))
    charged_sasa = np.nansum(amino_acid_sasa[charged_mask])     
    
    if total_sasa == 0:
         return 0.0
    
    return float(charged_sasa / total_sasa)
