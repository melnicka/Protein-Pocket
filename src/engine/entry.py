from ..utils.cif_parsing import load_cif, extract_ligands, extract_entry_metadata
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotite.structure import AtomArray

class Entry:
    """Representation of the whole PDB entry. 

    Attributes:
        pdb_id (str): A unique PDB identifier.
        atom_array (AtomArray): whole structure. 
        ligands (list[Ligand]): all ligands found in the entry.
    """
    def __init__(self, cif_file_path: str, pdb_id: str):
        self.pdb_id = pdb_id.upper()
        self.atom_array, self.sequences = load_cif(cif_file_path)
        self.ligands = []
        
        ligand_list = extract_ligands(self.atom_array)
        for lig_arr, lig_ids in ligand_list:
            self.ligands.append(Ligand(lig_arr, lig_ids))

    # TODO: will include pocket metadata in the future
    def extract_metadata(self) -> dict:
        """Extracts metadata from the entry.

        Returns:
            entry_metadata: Basic structural metadata and ligand identifiers.
        """
        metadata = {}
        entry_metadata = extract_entry_metadata(self.atom_array)
        ligand_ids = [lig.identifiers for lig in self.ligands]
        metadata['entry'] = entry_metadata
        metadata['ligands'] = ligand_ids
        return metadata


class Ligand:
    """Representation of the ligand molecule.

    Attributes:
        atom_array (AtomArray): structure of the ligand. 
        identifiers (dict): PDB identifiers.
        calculated_pockets (bool): True if pockets have been calculated.
        pocket (Pocket): The associated binding pocket.
    """
    def __init__(
        self,
        atom_array: AtomArray,
        identifiers: dict
    ):
        self.atom_array = atom_array
        self.identifiers = identifiers
        self.calculated_pockets = False
        self.pocket = None

#TODO: Implement Pocket class
class Pocket:
    pass
