import math
import numpy as np

class Membrane:
    def __init__(self, system, molecule, slice_size):
        self.frame = system.trajectory.frame
        self.time = system.trajectory.time

        self.slices = self.slice_membrane(system, molecule, slice_size)
        self.highest = max(self.slices)
        self.lowest = min(self.slices)

    def slice_membrane(self, system, molecule, slice_size):
        system.trajectory[0]
        membrane_atoms = system.select_atoms(f"resname {molecule}")
        max_x = system.dimensions[0]
        nb_slices = math.ceil(max_x / slice_size)
        
        slices = []
        prev = 0
        for i in range(1,nb_slices + 1):
            slice_atoms = membrane_atoms.select_atoms(f"prop x >= {prev * slice_size} and prop x < {i * slice_size}")
            prev = i
            z_mean = np.mean([atom.position[2] for atom in slice_atoms])
            slices.append(z_mean)

        return slices


