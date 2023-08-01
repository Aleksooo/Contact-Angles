import groio
from ..gro.Structure import Structure
from ..gro.Atom import Atom
import numpy as np

def read_gro(file_name):
    title, atoms, box = groio.parse_file(file_name)
    structure = Structure(
        title=title[:-1],
        box=np.fromstring(box, sep=' '),
        atoms=np.zeros(len(atoms), dtype=Atom),
        atoms_xyz=np.zeros((len(atoms), 3))
    )
    for i, a in enumerate(atoms):
        atom = Atom(
            mol_id=a['resid'],
            mol_name=a['resname'],
            name=a['atom_name'],
            id=a['atomid'],
            xyz=[a['x'], a['y'], a['z']]
        )

        structure.add_atom(atom, i)

    return structure


def write_gro(structure):
    gro = [structure.title, str(len(structure.atoms))]
    gro.extend(atom.make_gro_line() for atom in structure.atoms)
    gro.extend([f'{structure.box[0]:10.5f}{structure.box[1]:10.5f}{structure.box[2]:10.5f}\n'])

    return '\n'.join(gro)
