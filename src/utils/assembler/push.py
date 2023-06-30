import numpy as np
from tqdm import tqdm
from ..geom.pbc import distance_pbc2, delta_pbc

def push_atoms_apart(structure, min_distance2, max_distance2, iteration_lim=10):
    min_distance = np.sqrt(min_distance2)
    max_distance = np.sqrt(max_distance2)

    new_structure = structure.copy()
    ins = 0
    overlap = True
    while overlap:
        overlap = False
        overlap_counter = 0
        ins += 1
        if (ins > iteration_lim):
            break
        print(f"Iteration {ins}")
        for i, atom in enumerate(tqdm(new_structure.atoms[:-1], desc='Atoms')):
            for j, other_atom in enumerate(new_structure.atoms[i+1:]):
                if atom.mol_idx != other_atom.mol_idx:
                    delta = delta_pbc(structure.box, atom.xyz, [other_atom.xyz])[0]

                    if np.sum(abs(delta) < min_distance):
                        overlap = distance_pbc2(structure.box, atom.xyz, [other_atom.xyz])[0] < min_distance2
                        if overlap:
                            overlap_counter += 1
                            atom.xyz += delta / np.linalg.norm(delta) * \
                                np.random.uniform(min_distance, max_distance)
        print(f"{overlap_counter} overlaps detected")
    return new_structure.apply_pbc()
