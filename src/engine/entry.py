import numpy as np
import os
from ..utils.cif_parsing import load_cif, extract_ligands, extract_entry_metadata
from biotite.structure.io.pdbx import CIFFile, set_structure
from biotite.structure import AtomArray, CellList, get_residue_masks, filter_solvent, connect_via_residue_names
from biotite.interface.rdkit import to_mol
from rdkit.Chem.rdmolfiles import MolToSmiles


class Entry:
    """Representation of the whole PDB entry. 

    Attributes:
        pdb_id (str): A unique PDB identifier.
        atom_array (AtomArray): whole structure. 
        cif_file_path (str): Path to the mmCIF file.
        entry_dir (str): Path to the entry directory.
        ligands (list[Ligand]): All ligands found in the entry.
    """
    def __init__(self, cif_file_path: str, pdb_id: str):
        self.pdb_id = pdb_id.upper()
        self.atom_array = load_cif(cif_file_path)
        self.cif_file_path = cif_file_path
        self.entry_dir = os.path.dirname(cif_file_path)
        self.ligands = []
        
        ligand_list = extract_ligands(self.atom_array)
        for lig_arr, lig_ids in ligand_list:
            self.ligands.append(Ligand(lig_arr, lig_ids))

    def find_pockets(self, search_radius: float = 5.0, filter_out_solvent: bool = True):
        """Extracts pockets for all ligands. 

        Here binding pocket is defined as all residues within cutoff distance 
        from every ligand atom.

        Args:
            search_radius: The cutoff distance in Angstroms.
            filter_out_solvent: Whether to filter out solvent molecules.
        """
        cell_list = CellList(self.atom_array, search_radius)
        pocket_idx = 0
        for ligand in self.ligands:
            ligand_coord = ligand.atom_array.coord
            raw_indices = cell_list.get_atoms(ligand_coord, search_radius)
            indices = np.unique(raw_indices[raw_indices != -1])

            residue_masks = get_residue_masks(self.atom_array, indices)
            combined_residue_mask = np.logical_or.reduce(residue_masks) 
            pocket_array = self.atom_array[combined_residue_mask]
            if filter_out_solvent:
                pocket_array = pocket_array[~filter_solvent(pocket_array)]
            
            if len(pocket_array) == 0:
                ligand.pocket = Pocket(None, None)

            else: 
                pocket_array.bonds = connect_via_residue_names(pocket_array)

                pocket_file_name = f"{ligand.comp_id}_{ligand.auth_seq_id}_{ligand.auth_asym_id}_{pocket_idx}.cif"
                pocket_cif_file_path = f"{self.entry_dir}/{pocket_file_name}"

                ligand.pocket = Pocket(pocket_array, pocket_cif_file_path)

            pocket_idx += 1

    def save_pocket_cif_files(self):
        """Saves all pocket structures to mmCIF files."""
        for ligand in self.ligands:
            ligand.pocket.save()
            
    def extract_metadata(self) -> dict:
        """Extracts metadata from the entry.

        Returns:
            entry_metadata: Basic structural metadata and ligand identifiers.
        """
        entry_metadata = extract_entry_metadata(self.atom_array)
        ligand_ids = []
        for ligand in self.ligands:
            ligand_ids.append([ligand.comp_id, ligand.auth_seq_id, ligand.auth_asym_id])
        entry_metadata['ligand_identifiers'] = ligand_ids
        return entry_metadata

class Ligand:
    """Representation of the ligand molecule.

    Attributes:
        atom_array (AtomArray): Structure of the ligand. 
        comp_id (str): A chemical compound ID.
        auth_seq_id (str): The residue ID.
        auth_asym_id (str): The chain ID.
        pocket (Pocket): The associated binding pocket.
        smiles (str): The canonical SMILES for the ligand.
    """
    def __init__(
        self,
        atom_array: AtomArray,
        identifiers: dict
    ):
        self.atom_array = atom_array
        self.comp_id = identifiers['comp_id']
        self.auth_seq_id = identifiers['auth_seq_id']
        self.auth_asym_id = identifiers['auth_asym_id']

        self.pocket = None

        lig_mol = to_mol(atom_array)
        self.smiles = MolToSmiles(lig_mol)


class Pocket:
    """Representation of the ligand binding pocket.

    Attributes:
        atom_array (AtomArray | None): Structure of the pocket, including the ligand. 
        descriptors: (dict | None): Descriptors of the pocket.
        is_empty (bool): True if atom array could not be extracted.
        cif_file_path: Path to the saved pocket structure, ligand included.
    """
    def __init__(
            self,
            atom_array: AtomArray | None,
            cif_file_path: str | None,
    ):
        self.atom_array = atom_array
        if cif_file_path is None:
            self.is_empty = True
            self.cif_file_path = None
        else:
            self.is_empty = False
            self.cif_file_path = cif_file_path 

        self.descriptors = None

    def save(self):
        """Saves pocket's structure to a mmCIF file."""
        if self.cif_file_path is not None:
            cif_file = CIFFile()
            set_structure(cif_file, self.atom_array)
            cif_file.write(self.cif_file_path)


        
