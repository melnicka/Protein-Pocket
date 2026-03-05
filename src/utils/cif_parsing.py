import numpy as np
import biotite.structure as struct
from biotite.structure import AtomArray
import biotite.structure.io.pdbx as pdbx

def load_cif(file_path: str) -> AtomArray:
    """Loads mmCIF file to an AtomArray object. 

    Args:
        file_path: Path to the mmCIF file.

    Returns:
        atom_arr: AtomArray ready for further processing.
    """
    cif_file = pdbx.CIFFile.read(file_path)
    atom_arr = pdbx.get_structure(cif_file, model=1, extra_fields=['atom_id'])
    atom_arr.bonds = struct.connect_via_residue_names(atom_arr)
    return atom_arr

def extract_ligands(entry_arr: AtomArray) -> list[tuple[AtomArray, dict]]:
    """Extracts all ligands present in the entry, with their PDB identifiers.

    Args:
        entry_arr: A target entry.

    Returns:
        ligand_list: each tuple within the list represents one ligand:
            - ligand_arr: AtomArray of the ligand.
            - ligand_ids: PDB identifiers of the ligand.
    """
    ligands = []
    ligand_filter = (
    struct.filter_polymer(entry_arr, pol_type="peptide") |
    struct.filter_polymer(entry_arr, pol_type="nucleotide") |
    struct.filter_polymer(entry_arr, pol_type="carbohydrate") |
    struct.filter_solvent(entry_arr)
    )

    ligand_atoms = entry_arr[~ligand_filter]
    identifiers = set(
        zip(
            ligand_atoms.res_name,
            ligand_atoms.res_id,
            ligand_atoms.chain_id
        )
    )

    for ids in identifiers:
        mask = (
            (ligand_atoms.res_name == ids[0]) &
                (ligand_atoms.res_id == ids[1]) &
                (ligand_atoms.chain_id == ids[2])
        )
        ligand_arr = ligand_atoms[mask],
        ligand_ids = {
            'comp_id': ids[0],
            'auth_seq_id': ids[1],
            'auth_asym_id': ids[2]
        }
        ligands.append((ligand_arr, ligand_ids))
        

    return ligands

def extract_entry_metadata(entry_arr: AtomArray) -> dict:
    """Extracts basic metadata from PBD entry's AtomArray. 

    Args:
        atom_arr: A target entry.

    Returns:
        A dictionary with metadata about the entry.
    """
    full_atom_count = entry_arr.array_length()
    protein_atom_count = struct.filter_polymer(entry_arr, pol_type='peptide').sum()
    nuc_acid_atom_count = struct.filter_polymer(entry_arr, pol_type='nucleotide').sum()
    carb_atom_count = struct.filter_polymer(entry_arr, pol_type='carbohydrate').sum()
    ion_atom_count = struct.filter_monoatomic_ions(entry_arr).sum()
    solvent_atom_count = struct.filter_solvent(entry_arr).sum()
    chains = np.unique(struct.get_chains(entry_arr))
    residue_count = struct.get_residue_count(
            entry_arr[struct.filter_polymer(entry_arr)]
        )

    ligand_filter = (
    struct.filter_polymer(entry_arr, pol_type="peptide") |
    struct.filter_polymer(entry_arr, pol_type="nucleotide") |
    struct.filter_polymer(entry_arr, pol_type="carbohydrate") |
    struct.filter_solvent(entry_arr)
    )

    ligand_comp_ids = np.unique(entry_arr[~ligand_filter].res_name)

    return {
        'full_atom_count': full_atom_count,
        'protein_atom_count': protein_atom_count, 
        'nuc_acid_atom_count': nuc_acid_atom_count,
        'carb_atom_count': carb_atom_count,
        'ion_atom_count': ion_atom_count,
        'solvent_atom_count': solvent_atom_count,
        'chains': chains,
        'residue_count': residue_count,
        'ligand_comp_ids': ligand_comp_ids
    }






