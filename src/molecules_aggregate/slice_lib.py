import numpy as np

class Slice:
    def __init__(self, atoms):
        self.atoms = atoms
        self._z_mean = None

    @property
    def z_mean(self):
        if not self._z_mean:
            self._z_mean = np.mean([atom.position[2] for atom in self.atoms])
        return self._z_mean