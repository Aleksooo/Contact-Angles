import numpy as np 
from scipy.spatial.transform import Rotation
import sys
from ..geom.pbc import distance_pbc2
from tqdm import tqdm


def insert_com_into_shape(
    insertion_limit, mol_idx, mol_size, system_size, 
    package, old_coms, shape):
    overlap = True
    insertion_counter = 0
    while overlap:
        insertion_counter += 1
        if 5*insertion_counter % insertion_limit == 0:
            print(f'Please, wait ... mol #{mol_idx} is tried to be inserted\n')
        if insertion_counter > insertion_limit:
            sys.exit('Decrease package')
        new_com = shape.generate_point()
        if not len(old_coms):
            return new_com
        overlap = np.sum(
            distance_pbc2(system_size, new_com, old_coms) < \
            ((2 - insertion_counter/insertion_limit) * package * mol_size))
    return new_com


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


def find_position(rotation_limit, min_distance2, structure, new_com, mol):
    overlap = True
    rotation_counter = 0
    old_atoms = structure.get_XYZ()
    while overlap and rotation_counter < rotation_limit:
        rotation_counter += 1
        overlap = False
        new_mol = mol.copy()
        new_mol = new_mol.set_XYZ(rotate_random(new_mol.get_XYZ()) + new_com)
        if not len(old_atoms):
            return new_mol
        for atom in new_mol.get_XYZ():
            overlap = np.sum(distance_pbc2(structure.box, atom, old_atoms) < min_distance2)
            if overlap: 
                break
    return new_mol


def rotate_random(xyz):
    new_xyz = xyz.copy()
    theta = np.random.random(3)*2*np.pi
    for i in range(3):
        axis = np.zeros(3)
        axis[i] = 1
        rot = Rotation.from_rotvec(theta[i] * axis)
        new_xyz = rot.apply(new_xyz)  
    return new_xyz


