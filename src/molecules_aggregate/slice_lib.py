import numpy as np
import math 

class Slice:
    def __init__(self, atoms_position):
        self.atoms_position = atoms_position
        self._z_mean = None

    @property
    def z_mean(self):
        if not self._z_mean:
            self._z_mean = np.mean([pos[2] for pos in self.atoms_position])
        return self._z_mean

    @property
    def y_mean(self):
        return np.mean([pos[1] for pos in self.atoms_position])
    
    @property
    def x_mean(self):
        return np.mean([pos[0] for pos in self.atoms_position])

    @property
    def z_std(self):
        return np.std([pos[2] for pos in self.atoms_position])



    @property
    def nb_atoms(self):
        return len(self.atoms)

class Slices:
    def __init__(self, slice_intervals, slice_axis, atomgroup):
        #self.slice_intervals = slice_intervals
        self.slice_axis = slice_axis
        self.axis_idx = {'x':0, 'y':1, 'z':2}
        if not self.slice_axis in self.axis_idx:
            raise InvalidAxis(f"{axis} is not a valid value (must be x, y or z)")

        self.slices_list = self._compute_slices(slice_intervals, atomgroup)
        #print("A", self.slices_list[0].atoms[0].position)

    def __iter__(self):
        return iter(self.slices_list)

    def __getitem__(self, i):
        return self.slices_list[i]

    def __repr__(self):
        return str(self.slices_list)

    def _compute_slices(self, intervals, atomgroup):
        #print("AN", atomgroup.universe.anchor_name)
        #print("frame", atomgroup.universe.trajectory.frame)
        slices = []
        for i in intervals:
            slices.append(self._compute_single_slice(i, atomgroup))
        #print("frame", atomgroup.universe.trajectory.frame, slices[0].atoms[0].position)
        return slices

    def _compute_single_slice(self, interval, atomgroup):
        slice_atoms = atomgroup.select_atoms(f"prop x >= {interval[0]} and prop x < {interval[1]}")
        return Slice([a.position for a in slice_atoms])

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