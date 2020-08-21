import math, numpy as np
from . import config

class Membrane:
    def __init__(self, atomgroup, slice_intervals):
        slices_atom_positions = self._compute_slices(atomgroup, slice_intervals)
        means_std = self._get_means_std(slices_atom_positions)

        self.highest_point = self._get_highest_point(means_std)
        self.lowest_point = self._get_lowest_point(means_std)

    def _compute_slices(self, atomgroup, slice_intervals):
        """[summary]

        Args:
            atomgroup ([type]): [description]
            slice_intervals ([type]): [description]

        Returns:
            [type]: [description]
        """
        slices = []
        for i in slice_intervals:
            slices.append([atom.position for atom in atomgroup.select_atoms(f"prop {config.SLICE_AXIS} >= {i[0]} and prop {config.SLICE_AXIS} < {i[1]}")])
        return slices
    
    def _get_means_std(self, slices_atom_positions):
        """[summary]

        Args:
            slices_atom_positions ([type]): [description]

        Returns:
            [type]: [description]
        """
        means = {"x":{'mean':[], 'std':[]}, "y":{'mean':[], 'std':[]}, "z":{'mean':[], 'std':[]}}
        for sl in slices_atom_positions:
            means["x"]["mean"].append(np.mean([pos[0] for pos in sl]))
            means["y"]["mean"].append(np.mean([pos[1] for pos in sl]))
            means["z"]["mean"].append(np.mean([pos[2] for pos in sl]))
            means["x"]["std"].append(np.std([pos[0] for pos in sl]))
            means["y"]["std"].append(np.std([pos[1] for pos in sl]))
            means["z"]["std"].append(np.std([pos[2] for pos in sl]))
        return means

    def _get_highest_point(self, mean_std):
        """[summary]

        Args:
            mean_std ([type]): [description]

        Raises:
            Exception: [description]

        Returns:
            [type]: [description]
        """
        max_mean = 0
        max_idx = -1
        for i, mean in enumerate(mean_std[config.MEMBRANE_AXIS]["mean"]):
            if mean > max_mean:
                max_idx = i
                max_mean = mean

        if max_idx == -1:
            raise Exception("Highest point not found")

        return {"x" : {"mean" : mean_std["x"]["mean"][max_idx], "std" : mean_std["x"]["std"][max_idx]},\
            "y": {"mean" : mean_std["y"]["mean"][max_idx], "std" : mean_std["y"]["std"][max_idx]},\
            "z" : {"mean" : mean_std["z"]["mean"][max_idx], "std" : mean_std["z"]["std"][max_idx]},\
            "i" : max_idx}

    def _get_lowest_point(self, mean_std):
        """[summary]

        Args:
            mean_std ([type]): [description]

        Raises:
            Exception: [description]

        Returns:
            [type]: [description]
        """
        min_mean = float('inf')
        min_idx = -1
        for i, mean in enumerate(mean_std[config.MEMBRANE_AXIS]["mean"]):
            if mean < min_mean:
                min_idx = i
                min_mean = mean

        if min_idx == -1:
            raise Exception("Lowest point not found")

        return {"x" : {"mean" : mean_std["x"]["mean"][min_idx], "std" : mean_std["x"]["std"][min_idx]},\
            "y": {"mean" : mean_std["y"]["mean"][min_idx], "std" : mean_std["y"]["std"][min_idx]},\
            "z" : {"mean" : mean_std["z"]["mean"][min_idx], "std" : mean_std["z"]["std"][min_idx]},\
            "i" : min_idx}

    

def get_slice_intervals(total_size, slice_size):
    """[summary]

    Args:
        total_size ([type]): [description]
        slice_size ([type]): [description]

    Returns:
        [type]: [description]
    """
    nb_slices = math.ceil(total_size / slice_size)
    intervals = []
    prev = 0
    for i in range(1, nb_slices + 1):
        intervals.append((prev * slice_size, i*slice_size))
        prev = i
    return intervals
