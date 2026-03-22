import pytest
import biotite.structure as struct
from engine.descriptors import calc_ligand_buried_surface
from engine.entry import Entry
from utils.data_fetching import fetch_cif

@pytest.fixture(scope='module')
def testing_setup():
    cif_path = fetch_cif('4HHB')
    entry = Entry(cif_path, '4HHB')
    entry.find_pockets()
    pocket_arrays = []
    ligand_arrays = []
    for ligand in entry.ligands:
        pocket_arrays.append(ligand.pocket.atom_array)
        ligand_arrays.append(ligand.atom_array)

    return (pocket_arrays, ligand_arrays)

def test_ligand_buried_surface_far_away(testing_setup):
    pocket_arrays, ligand_arrays = testing_setup
    for pocket, ligand in zip(pocket_arrays, ligand_arrays):
        test_pocket = ligand.copy()
        ligand_in_pocket_mask = pocket.res_name == ligand.res_name[0]
        pocket[ligand_in_pocket_mask].coord += 100

    near = calc_ligand_buried_surface(pocket, ligand)
    far = calc_ligand_buried_surface(test_pocket, ligand)

    assert near > far
    assert far == pytest.approx(0.0, abs=1e-4)

def test_ligand_buried_surface_valid_float(testing_setup):
    pocket_arrays, ligand_arrays = testing_setup
    for pocket, ligand in zip(pocket_arrays, ligand_arrays):
        buried = calc_ligand_buried_surface(pocket, ligand)
        assert buried >= 0 




