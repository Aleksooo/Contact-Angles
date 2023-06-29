import groio
from ..gro.Structure import Structure
from ..gro.Atom import Atom
import numpy as np

def read_gro(file_name):
    title, atoms, box = groio.parse_file(file_name)
    structure = Structure(title=title[:-1], box=np.fromstring(box, sep=' '))
    for a in atoms:
        atom = Atom(
            mol_idx=a['resid'],
            mol_name=a['resname'],
            name=a['atom_name'],
            idx=a['atomid'],
            xyz=[a['x'], a['y'], a['z']]
            )
        structure.atoms.append(atom)
    return structure


def write_gro(structure):
    gro = [structure.title, str(len(structure.atoms))]
    gro.extend(atom.make_gro_line() for atom in structure.atoms)
    gro.extend([f'{structure.box[0]:10.5f}{structure.box[1]:10.5f}{structure.box[2]:10.5f}\n'])
    return '\n'.join(gro)
