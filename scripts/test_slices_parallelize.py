
import MDAnalysis as mda 
import math
import time
from pmda.parallel import ParallelAnalysisBase
from MDAnalysis.lib import distances
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

def slice_intervals(total_size, slice_size):
    nb_slices = math.ceil(total_size / slice_size)
    intervals = []
    prev = 0
    for i in range(1, nb_slices + 1):
        intervals.append((prev * slice_size, i*slice_size))
        prev = i
    return intervals

def get_slices(atom_group):
    #intervals = slice_intervals(atom_group.universe.dimensions[0], 10)
    some_group = atom_group.select_atoms(f"prop x >= {0} and prop x < {100}")
    return some_group

def serial_slice(atom_group, intervals):
    slices = []
    for i in intervals:
        slices.append(get_slice(atom_group, i))
    return slices

def center_of_mass(ag):
    ag[0].universe.trajectory[ag[1]]
    print(ag[0][0].position)


def multiprocessing(func, workers, args):
    with ProcessPoolExecutor(workers) as ex:
        res = ex.map(func, args)
    return list(res)


if __name__ == "__main__":
    system = mda.Universe("../data/gro456.PO4.tpr", "../data/md.PO4.xtc")
    membrane_atoms = system.select_atoms(f"resname DOPC")
    test_atoms = membrane_atoms.select_atoms("prop x >= 0 and prop x < 10", updating = True)
    intervals = slice_intervals(system.dimensions[0], 10)
    test_slices = []
    for i in intervals:
        test_slices.append(membrane_atoms.select_atoms(f"prop x >= {i[0]} and prop x < {i[1]}", updating = True))
    
    for t in system.trajectory[:2]:
        print(len(test_slices[0]))