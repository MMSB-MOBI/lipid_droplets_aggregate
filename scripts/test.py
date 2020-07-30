#!/home/chilpert/miniconda3/envs/mdanalysis/bin/python

import molecules_aggregate
import MDAnalysis as mda
import time 
import cProfile, pstats, io
from pstats import SortKey

'''pr = cProfile.Profile()
pr.enable()'''

start = time.time()

mda_system = mda.Universe("../data/gro456.PO4.tpr", "../data/md.PO4.xtc")
mol = "TO"
cluster_threshold = 13
to_keep = 2
membrane_mol = "DOPC"
slice_membrane_size = 10

system = molecules_aggregate.load_system(mda_system, mol, cluster_threshold, to_keep, membrane_mol, slice_membrane_size)

'''pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())'''

print("END", time.time() - start)