import numpy as np
from tqdm import tqdm
import sys
from src.utils_py.io.gro import read_gro
from src.utils_py.assembler.insert import insert_point_into_shape, find_position
from src.utils_py.assembler.push import push_atoms_apart
from .grid import add_mol_grid

def build_system(
    dir, structure, names=None, numbers=None, shapes=None, points = list(),
    insertion_limit = int(1e5),
    rotation_limit = 10,
    package = 0.4,
    min_dist2 = 0.08**2
):
    if names == None or numbers == None or shapes == None:
        sys.exit('Please, give all the parametrs')

    if not(len(names) == len(numbers) == len(shapes)):
        sys.exit('Parameters have different lengths')


    # structure = Structure(box=system_size, atoms=[])
    system_size = structure.box
    print('Number of molecules:')
    for i, name in enumerate(names):
        print(name, '\t', numbers[i])

    # N = np.ceil(structure.box / np.sqrt(min_dist2)).astype(int)
    # dr = structure.box / N
    # grid = np.array(np.empty(N, dtype=np.object_))
    # for i in range(grid.shape[0]):
    #     for j in range(grid.shape[1]):
    #         for k in range(grid.shape[2]):
    #             grid[i, j, k] = []

    print('\nFilling system:')
    for i, name in enumerate(names):
        mol = read_gro(f'{dir}/gro/{name}.gro').center_atoms_to_center().center_atoms_to_zero()
        mol_size = np.max(np.linalg.norm(mol.get_XYZ(), axis=1))

        for mol_id in tqdm(range(numbers[i])):
            new_point = insert_point_into_shape(shapes[i], points, system_size, mol_size, mol_id, insertion_limit, package)
            points.append(new_point)

            new_mol = find_position(structure, new_point, mol, mol_id, rotation_limit, min_dist2)
            # new_mol = find_position(structure, new_point, mol, mol_id, grid, N, dr, rotation_limit, min_dist2)
            for atom in new_mol.atoms:
                new_atom = atom.copy()
                new_atom.id = len(structure.atoms) + 1
                new_atom.mol_id = mol_id + 1
                structure.atoms.append(new_atom)

                # grid = add_mol_grid(grid, structure, new_atom.xyz, new_atom.id, N, dr)

    return structure
