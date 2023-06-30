from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from scipy.spatial.transform import Rotation
from .Shape import Shape

@dataclass
class Cylinder(Shape):
    radius: float = 0
    length: float = 0
    axis: np.array = field(default_factory=np.array)

    def get_volume(self) -> float:
        return np.pi * self.radius**2 * self.length

    def get_surface(self) -> float:
        return 2 * np.pi * self.radius * (self.radius + self.length)

    def check_point(self, point) -> bool:
        d = np.linalg.norm(np.cross(point - self.center, self.axis)) / np.linalg.norm(self.axis)
        l = np.dot(self.axis, point - self.center) / np.linalg.norm(self.axis)

        return d < self.radius and l < self.length / 2

    def generate_point(self) -> np.array:
        phi = np.random.uniform(0, 2*np.pi)
        theta = np.random.uniform(0, np.pi)
        sphere_vec = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

        r_vec = sphere_vec - np.dot(sphere_vec, self.axis) * self.axis / np.linalg.norm(self.axis)**2
        r = np.random.uniform(0, self.radius)
        r_vec = r * r_vec / np.linalg.norm(r_vec)

        psi = np.random.uniform(0, 2*np.pi)
        z = np.random.uniform(-self.length / 2, self.length / 2)
        return self.center + z * self.axis / np.linalg.norm(self.axis) + r_vec
