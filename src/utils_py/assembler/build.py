import numpy as np
from tqdm import tqdm
import sys
from src.utils_py.io.gro import read_gro
from src.utils_py.assembler.insert import insert_point_into_shape, find_position
from src.utils_py.assembler.push import push_atoms_apart

def build_system(
    dir_gro, structure, names=None, density=None, shapes=None, points = list(),
    insertion_limit = int(1e5),
    rotation_limit = 10,
    package = 0.4,
    min_dist2 = 0.08**2
):
    if names == None or density == None or shapes == None:
        sys.exit('Please, give all the parametrs')

    if not(len(names) == len(density) == len(shapes)):
        sys.exit('Parameters have different lengths')

    # structure = Structure(box=system_size, atoms=[])
    system_size = structure.box
    numbers = dict(zip(names, np.round([shapes[name].get_volume() * density[name] for name in names]).astype(int)))
    print('Number of molecules:')
    for key, item in numbers.items():
        print(key, '\t', item)

    print('\nFilling system:')
    for name in names:
        mol = read_gro(f'{dir_gro}{name}.gro').center_atoms_to_center().center_atoms_to_zero()
        mol_size = np.max(np.linalg.norm(mol.get_XYZ(), axis=1))

        for mol_id in tqdm(range(numbers[name])):
            new_point = insert_point_into_shape(shapes[name], points, system_size, mol_size, mol_id, insertion_limit, package)
            points.append(new_point)

            new_mol = find_position(structure, new_point, mol, mol_id, rotation_limit, min_dist2)
            for atom in new_mol.atoms:
                new_atom = atom.copy()
                new_atom.id = len(structure.atoms) + 1
                new_atom.mol_id = mol_id + 1
                structure.atoms.append(new_atom)

    return structure
