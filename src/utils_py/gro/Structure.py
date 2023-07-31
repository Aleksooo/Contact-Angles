from dataclasses import dataclass, field
import numpy as np
from ..assembler.pbc import apply_pbc_to_points

@dataclass
class Structure:
    title: str = 'SYSTEM'
    box: np.array = field(default_factory=np.array)
    atoms: list = field(default_factory=list)

    def get_center(self) -> np.array:
        return self.box/2

    def get_XYZ(self, mol_names=None) -> np.array:
        if mol_names is None:
            return np.array([atom.xyz for atom in self.atoms])

        return np.array([atom.xyz for atom in self.atoms if atom.mol_name in mol_names])

    # ???
    def get_center_pbc(self, mol_names=None):
        cos_theta = np.zeros(3)
        sin_theta = np.zeros(3)
        center_pbc = np.zeros(3)

        for i in range(3):
            theta = self.get_XYZ(mol_names)[:, i]/self.box[i] * 2.0 * np.pi
            cos_theta[i] = np.sum(np.cos(theta))
            sin_theta[i] = np.sum(np.sin(theta))
            center_pbc[i] = (
                np.pi + \
                np.arctan2(
                    -sin_theta[i] / len(self.atoms),
                    -cos_theta[i] / len(self.atoms))) * \
                self.box[i] / (2.0 * np.pi)

        return center_pbc

    def set_XYZ(self, new_coords):
        new_structure = self.copy()
        for i, atom in enumerate(self.atoms):
            new_structure.atoms[i].xyz = new_coords[i, :]

        return new_structure

    def apply_pbc(self):
        return self.set_XYZ(apply_pbc_to_points(self.get_XYZ(), self.box))

    def center_atoms_to_zero(self, mol_names=None):
        return self.set_XYZ(self.get_XYZ() - self.get_center_pbc(mol_names))

    def center_atoms_to_center(self, mol_names=None):
        new_coords = self.get_XYZ() - self.get_center_pbc(mol_names) + self.box/2
        return self.set_XYZ(new_coords).apply_pbc()

    def copy(self):
        return Structure(
            title = self.title,
            box = self.box.copy(),
            atoms = [atom.copy() for atom in self.atoms]
        )
