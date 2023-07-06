from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from scipy.spatial.transform import Rotation
from .Shape import Shape
from .Box import Box
from .Cylinder import Cylinder

@dataclass
class AntiCylinder(Shape):
    radius: float = 0
    length: float = 0
    axis: np.array = field(default_factory=np.array)
    borders_center: np.array = field(default_factory=np.array)
    borders:  np.array = field(default_factory=np.array)

    def get_volume(self) -> float:
        return self.get_box().get_volume() - self.get_cylinder().get_volume()

    def get_surface(self) -> float:
        return self.get_box().get_surface() + self.get_cylinder().get_surface()

    def check_point(self, point) -> bool:
        return self.get_box().check_point(point) and (not self.get_cylinder().check_point(point))

    def generate_point(self) -> np.array:
        inside_cylinder = True
        while inside_cylinder:
            point = self.get_box().generate_point()
            inside_cylinder = self.get_cylinder().check_point(point)
        return point

    def get_box(self):
        return Box(center=self.borders_center, borders=self.borders)

    def get_cylinder(self):
        return Cylinder(center=self.center, radius=self.radius, length=self.length, axis=self.axis)
