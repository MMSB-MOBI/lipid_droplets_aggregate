import math
import numpy as np
from . import slice_lib
import asyncio
from multiprocessing import Pool

class Membrane:
    def __init__(self, parent_system, slice_intervals, axis, atomgroup):
        self.parent_system = parent_system
        self.slices = slice_lib.Slices(slice_intervals, axis, atomgroup)
        self.highest = self.slices.highest_z_slice()
        self.lowest = self.slices.lowest_z_slice()
        #self.highest_z_mean = max([sl.z_mean for sl in self.slices])
        #self.lowest_z_mean = min([sl.z_mean for sl in self.slices])
