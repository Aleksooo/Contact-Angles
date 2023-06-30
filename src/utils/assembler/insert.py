import numpy as np
from scipy.spatial.transform import Rotation
import sys
from ..geom.pbc import distance_pbc2
from tqdm import tqdm


# def insert_com_into_shape(
#     insertion_limit, system_size,
#     package, old_coms, shape, mol_size, mol_idx):
#     overlap = True
#     insertion_counter = 0
#     while overlap:
#         insertion_counter += 1
#         if 5*insertion_counter % insertion_limit == 0:
#             print(f'Please, wait ... mol #{mol_idx} is tried to be inserted\n')
#         if insertion_counter > insertion_limit:
#             sys.exit('Decrease package')

#         new_com = shape.generate_point()
#         if not len(old_coms):
#             return new_com
#         overlap = np.sum(
#             distance_pbc2(system_size, new_com, old_coms) < \
#             ((2 - insertion_counter/insertion_limit) * package * mol_size))
#     return new_com

def insert_point_into_shape(system_size, points, shape, mol_size, mol_idx, insertion_limit, package):
    overlap = True
    insertion_counter = 0

    if not mol_idx:
        return shape.generate_point()

    while overlap:
        insertion_counter += 1
        if 5*insertion_counter % insertion_limit == 0:
            print(f'Please, wait ... mol #{mol_idx} is tried to be inserted\n')
        if insertion_counter > insertion_limit:
            sys.exit('Decrease package')

        new_point = shape.generate_point()
        overlap = np.sum(
            distance_pbc2(system_size, new_point, points) < \
            ((2 - insertion_counter/insertion_limit) * package * mol_size)
        )
    return new_point

def insert_com(insertion_limit, mol_idx, mol_size, system_size, package, old_coms):
    overlap = True
    insertion_counter = 0
    while overlap:
        insertion_counter += 1
        if 5*insertion_counter % insertion_limit == 0:
            print(f'Please, wait ... mol #{mol_idx} is tried to be inserted\n')
        if insertion_counter > insertion_limit:
            sys.exit('Decrease package')
        new_com = np.random.random(3)*system_size
        if not len(old_coms):
            return new_com
        overlap = np.sum(
            distance_pbc2(system_size, new_com, old_coms) < \
            ((2 - insertion_counter/insertion_limit) * package * mol_size))
    return new_com


def find_position(structure, new_point, mol, mol_idx, rotation_limit, min_distance2):
    overlap = True
    rotation_counter = 0
    old_atoms = structure.get_XYZ()

    if not mol_idx:
        new_mol = mol.copy()
        new_mol = new_mol.set_XYZ(rotate_random(new_mol.get_XYZ()) + new_point)
        return new_mol

    while overlap and rotation_counter < rotation_limit:
        # overlap = False
        rotation_counter += 1

        new_mol = mol.copy()
        new_mol = new_mol.set_XYZ(rotate_random(new_mol.get_XYZ()) + new_point)

        for atom in new_mol.get_XYZ():
            overlap = np.sum(distance_pbc2(structure.box, atom, old_atoms) < min_distance2)
            if overlap:
                break
    return new_mol

def rotate_random(xyz):
    theta = np.random.random(3)*2*np.pi
    rot = Rotation.from_euler('xyz', theta)
    new_xyz = rot.apply(xyz.copy())
    return new_xyz
