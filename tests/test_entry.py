import pytest
import biotite.structure as struct
from engine.entry import Entry
from utils.data_fetching import fetch_cif

@pytest.fixture(scope='module')
def testing_entry():
    cif_path = fetch_cif('4HHB')
    entry = Entry(cif_path, '4HHB')
    entry.find_pockets()

    return entry

def test_entry_ligand_atom_array(testing_entry):
    if testing_entry.ligands:
        for ligand in testing_entry.ligands:
            assert isinstance(ligand.atom_array, struct.AtomArray)
            

def test_entry_ligand_pocket_atom_array(testing_entry):
    if testing_entry.ligands:
        for ligand in testing_entry.ligands:
            pocket_array = ligand.pocket.atom_array
            assert isinstance(
                    pocket_array,
                    (struct.AtomArray, type(None))
            )
        

