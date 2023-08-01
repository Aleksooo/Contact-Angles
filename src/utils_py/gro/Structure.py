from dataclasses import dataclass, field
import numpy as np
from ..assembler.pbc import apply_pbc_to_points

@dataclass
class Structure:
    title: str = 'SYSTEM'
    box: np.array = field(default_factory=np.array)
    atoms: np.array = field(default_factory=np.array)
    atoms_xyz: np.array = field(default_factory=np.array)

    def add_atom(self, atom, atom_id):
        self.atoms[atom_id] = atom
        self.atoms_xyz[atom_id, :] = atom.xyz

    def get_center(self) -> np.array:
        return self.box/2

    def get_XYZ(self, mol_names=None) -> np.array:
        if mol_names is None:
            return self.atoms_xyz

        return np.array([atom.xyz for atom in self.atoms if atom.mol_name in mol_names])

    def get_center_pbc(self, mol_names=None):
        theta = self.atoms_xyz / self.box * 2 * np.pi
        center_pbc = np.zeros(3)

        for i in range(3):
            cos_phi_mean = np.average(np.cos(theta[:, i]))
            sin_psi_mean = np.average(np.sin(theta[:, i]))

            theta_mean = np.arctan2(-sin_psi_mean, -cos_phi_mean) + np.pi
            center_pbc[i] = self.box[i] * theta_mean / 2 / np.pi

        return center_pbc

    def set_XYZ(self, new_coords):
        new_structure = self.copy()
        for i, atom in enumerate(self.atoms):
            new_structure.atoms[i].xyz = new_coords[i, :]

        new_structure.atoms_xyz = new_coords.copy()

        return new_structure

    def apply_pbc(self):
        return self.set_XYZ(apply_pbc_to_points(self.atoms_xyz, self.box))

    def center_atoms_to_zero(self, mol_names=None):
        return self.set_XYZ(self.atoms_xyz - self.get_center_pbc(mol_names))

    def center_atoms_to_center(self, mol_names=None):
        new_coords = self.atoms_xyz - self.get_center_pbc(mol_names) + self.box/2
        return self.set_XYZ(new_coords).apply_pbc()

    def copy(self):
        return Structure(
            title = self.title,
            box = self.box.copy(),
            atoms = self.atoms.copy(),
            atoms_xyz =  self.atoms_xyz.copy()
        )
