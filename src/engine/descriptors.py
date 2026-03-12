import biotite.structure as struct 

def calc_ligand_buried_surface(
        pocket_array: struct.AtomArray,
        ligand_array: struct.AtomArray
) -> float:
    ligand_sasa = struct.sasa(ligand_array).sum()
    ligand_in_pocket_mask = pocket_array.res_name == ligand_array.res_name[0] 
    pocket_sasa = struct.sasa(pocket_array)
    ligand_in_pocket_sasa = pocket_sasa[ligand_in_pocket_mask].sum()

    return ligand_sasa - ligand_in_pocket_sasa



