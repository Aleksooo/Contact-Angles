from dataclasses import dataclass, field
import numpy as np

@dataclass
class Atom:
    mol_idx: int = 1
    mol_name: str = 'MOL'
    name: str = 'ATOM'
    idx: int = 1
    xyz: np.array = field(default_factory=np.array)
    
    def make_gro_line(self) -> str:
        start = f'{self.mol_idx:5d}{self.mol_name:<5}{self.name:>5}{self.idx:5d}' 
        coords = ''.join([f'{x:8.3f}' for x in self.xyz])
        return f'{start}{coords}'

    def copy(self):
        return Atom(
            mol_idx = self.mol_idx,
            mol_name = self.mol_name,
            name = self.name,
            idx = self.idx,
            xyz = self.xyz.copy())