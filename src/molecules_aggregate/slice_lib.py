import numpy as np
import math 

class Slice:
    def __init__(self, atoms):
        self.atoms = atoms
        self._z_mean = None

    @property
    def z_mean(self):
        if not self._z_mean:
            self._z_mean = np.mean([atom.position[2] for atom in self.atoms])
        return self._z_mean

    @property
    def y_mean(self):
        return np.mean([atom.position[1] for atom in self.atoms])
    
    @property
    def x_mean(self):
        return np.mean([atom.position[0] for atom in self.atoms])

    @property
    def z_std(self):
        return np.std([atom.position[2] for atom in self.atoms])



    @property
    def nb_atoms(self):
        return len(self.atoms)

class Slices:
    def __init__(self, parent_system, slice_size, slice_axis):
        self.slice_size = slice_size
        self.slice_axis = slice_axis
        self.parent_system = parent_system
        self.axis_idx = {'x':0, 'y':1, 'z':2}
        if not self.slice_axis in self.axis_idx:
            raise InvalidAxis(f"{axis} is not a valid value (must be x, y or z)")

        self.slices_list = self._compute_slices()

    def __iter__(self):
        return iter(self.slices_list)

    def __getitem__(self, i):
        return self.slices_list[i]

    def __repr__(self):
        return str(self.slices_list)

    def _compute_slices(self):
        max_slice = self.parent_system.box_dimension[self.axis_idx[self.slice_axis]]

        intervals = slice_intervals(max_slice, self.slice_size)
        slices = []
        for i in intervals:
            slices.append(self._compute_single_slice(i))

        return slices

    def _compute_single_slice(self, interval):
        slice_atoms = self.parent_system.membrane_atom_group.select_atoms(f"prop x >= {interval[0]} and prop x < {interval[1]}")
        return Slice(slice_atoms)

    def highest_z_slice(self):
        #Maybe try shorter version
        max_mean = 0
        max_idx = -1
        for i,s in enumerate(self):
            if s.z_mean > max_mean:
                max_idx = i
                max_mean = s.z_mean

        if max_idx == -1 : 
            raise Exception("Highest slice doesn't found.")

        return self[max_idx]

    def lowest_z_slice(self):
        min_mean = float('inf')
        min_idx = -1
        for i,s in enumerate(self):
            if s.z_mean < min_mean:
                min_idx = i
                min_mean = s.z_mean

        if min_idx == -1 : 
            raise Exception("Highest slice doesn't found.")

        return self[min_idx]


def slice_intervals(total_size, slice_size):
    nb_slices = math.ceil(total_size / slice_size)
    intervals = []
    prev = 0
    for i in range(1, nb_slices + 1):
        intervals.append((prev * slice_size, i*slice_size))
        prev = i
    return intervals
