from ..utils.cif_parsing import load_cif, extract_ligands, extract_entry_metadata
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotite.structure import AtomArray

class Entry:
    def __init__(self, cif_file_path: str, pdb_id: str):
        self.pdb_id = pdb_id
        self.atom_array = load_cif(cif_file_path)
        self.ligands = []
        
        ligand_list = extract_ligands(self.atom_array)
        for lig_arr, lig_ids in ligand_list:
            self.ligands.append(Ligand(lig_arr, lig_ids))

    def extract_metadata(self) -> dict:
        """Extracts metadata from the entry.

        Returns:
            entry_metadata: Basic structural metadata and ligand identifiers.
        """
        entry_metadata = extract_entry_metadata(self.atom_array)
        ligand_ids = [lig.identifiers for lig in self.ligands]
        entry_metadata['ligand_identifiers'] = ligand_ids
        return entry_metadata


class Ligand:
    def __init__(
        self,
        atom_array: AtomArray,
        identifiers: dict
    ):
        self.atom_array = atom_array
        self.identifiers = identifiers
        self.calculated_pockets = False
        self.pockets = None
