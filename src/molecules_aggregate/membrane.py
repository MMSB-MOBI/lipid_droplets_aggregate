import math
import numpy as np
from . import slice_lib
import asyncio
from multiprocessing import Pool

class Membrane:
    def __init__(self, atom_group, slice_size, axis, box_dimension):
        self.slices = slice_atom_group(atom_group, slice_size, axis, box_dimension)
        self.highest_z_mean = max([sl.z_mean for sl in self.slices])
        self.lowest_z_mean = min([sl.z_mean for sl in self.slices])

def slice_atom_group(atom_group, slice_size, axis, box_dimension):
    axis_idx = {'x':0, 'y':1, 'z':2}

    if not axis in axis_idx:
        raise InvalidAxis(f"{axis} is not a valid value (must be x, y or z)")

    max_slice = box_dimension[axis_idx[axis]]

    intervals = slice_intervals(max_slice, slice_size)

    with Pool(6) as p:
        slices = p.starmap(run_select_single, [(atom_group, interval) for interval in intervals])

    return slices

def slice_intervals(total_size, slice_size):
    nb_slices = math.ceil(total_size / slice_size)
    intervals = []
    prev = 0
    for i in range(1, nb_slices + 1):
        intervals.append((prev * slice_size, i*slice_size))
        prev = i
    return intervals


def run_select_single(atom_group, interval):
    slice_atoms = atom_group.select_atoms(f"prop x >= {interval[0]} and prop x < {interval[1]}")
    return slice_lib.Slice(slice_atoms)


