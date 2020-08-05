from . import clusters as clusters_lib
from . import membrane as membrane_lib
from . import error
import logging
import pmda.custom
import math, multiprocessing, time

global COPIED_UNIVERSES

class System:
    def __init__(self, parent, mda_system, membrane_atoms, previous = None):
        self.parent = parent
        self.frame = mda_system.trajectory.frame
        self.time = mda_system.trajectory.time
        self.previous = previous
        self.membrane_atoms = membrane_atoms
        self.clusters = None
        self.membrane = None
        #self.pouet = mda_system

    @property
    def box_dimension(self):
        return self.parent.box_dimension

    @property
    def cluster_atom_group(self):
        return self.parent.cluster_atom_group

    @property
    def cluster_threshold(self):
        return self.parent.cluster_threshold

    @property
    def to_keep_clusters(self):
        return self.parent.to_keep_clusters

    @property
    def nb_cluster(self):
        return len(self.clusters)

    @property
    def correspondance_method(self):
        return self.parent.correspondance_method

    @property
    def nb_correspondance(self):
        return self.parent.nb_correspondance

    def load_clusters(self):
        self.clusters = clusters_lib.OneFrameClusters(self)

    def load_membrane(self):
        self.membrane = membrane_lib.Membrane(self, self.parent.slice_intervals, self.parent.slice_axis, self.membrane_atoms)

    def add_previous(self, previous_system):
        self.previous = previous_system
    

class SystemIterator:
    def __init__(self, mda_system, cluster_mol_name, cluster_threshold, to_keep_clusters, membrane_mol_name, slice_membrane_size, nb_frames, correspondance_method, nb_correspondance, slice_axis = 'x'):
        self.cluster_atom_group = mda_system.select_atoms(f"resname {cluster_mol_name}")
        #self.membrane_atom_group = mda_system.select_atoms(f"resname {membrane_mol_name}")
        if not self.cluster_atom_group:
            raise error.MoleculeNotFound(f"Molecule {cluster_mol_name} doesn't exist in system")
        #if not self.membrane_atom_group:
        #    raise error.MoleculeNotFound(f"Molecule {membrane_mol_name} doesn't exist in system")

        self.box_dimension = mda_system.dimensions
        self.cluster_threshold = cluster_threshold
        self.to_keep_clusters = to_keep_clusters
        self.slice_membrane_size = slice_membrane_size
        self.slice_axis = slice_axis
        self.slice_intervals = slice_intervals(self.box_dimension[0], self.slice_membrane_size) # Get from axis ! 
        #self.mda_system = mda_system
        self.correspondance_method = correspondance_method
        self.nb_correspondance = nb_correspondance

        self.system_frames = self._load_frames(mda_system, nb_frames)

    @property 
    def nb_clusters(self):
        return self.system_frames[0].nb_cluster

    def __iter__(self):
        return iter(self.system_frames)

    def __getitem__(self,i):
        return self.system_frames[i]

    def _load_frames(self, mda_system, nb_frames, nb_process = 32):
        global COPIED_UNIVERSES
        #1. Load membrane
        logging.info("== Load membranes...")
        nb_frames = len(mda_system.trajectory) if nb_frames == -1 else nb_frames
        trajectory_chunks = slice_trajectories(nb_frames, nb_process)

        COPIED_UNIVERSES = copy_systems(trajectory_chunks, mda_system)
        
        p = multiprocessing.Pool(nb_process)

        packed_frames = p.map(self._some_frames_membrane, [(anchor) for anchor in COPIED_UNIVERSES])
        frames = [frame for process in packed_frames for frame in process ]
        p.close()
        p.join()

        #2. Register previous frames
        logging.info("== Register previous frames...")
        for i, frame in enumerate(frames[1:]): 
            frame.add_previous(frames[i])
        
        #3. Load clusters
        logging.info("Load clusters...")
        for trj in mda_system.trajectory[:nb_frames]:
            if trj.frame % 1000 == 0:
                logging.info(f"Frame {trj.frame}/{nb_frames}")
            frames[trj.frame].load_clusters()

        return frames

    def _some_frames_membrane(self, universe_idx):
        global COPIED_UNIVERSES
        frames = []
        universe = COPIED_UNIVERSES[universe_idx]["universe"]
        for trj in universe.trajectory[COPIED_UNIVERSES[universe_idx]["i"]:COPIED_UNIVERSES[universe_idx]["j"]]:
            frame = System(self, universe, COPIED_UNIVERSES[universe_idx]["membrane_atoms"])
            frame.load_membrane()
            frames.append(frame)
        return frames
        

def slice_trajectories(nb_frames, nb_process):
    nb_frames_by_pool = int(nb_frames / nb_process)
    nb_chunks = nb_frames if nb_frames_by_pool == 0 else nb_process
    #nb_frames_by_pool = 1 if nb_frames_by_pool == 0 else nb_frames_by_pool

    intervals = []
    rest = nb_frames - nb_frames_by_pool * nb_process

    i = 0
    for j in range(rest):
        intervals.append([i, i + nb_frames_by_pool + 1])
        i = i + nb_frames_by_pool + 1
    
    for j in range(nb_chunks - rest):
        intervals.append([i, i + nb_frames_by_pool])
        i = i + nb_frames_by_pool
    
    
    return intervals

    

def copy_systems(chunks_interval, system):
    dict_store = {i:{} for i in range(len(chunks_interval))}
    for i,frame_int in enumerate(chunks_interval):
        copied_universe = system.copy()
        membrane_atoms = copied_universe.select_atoms(f"resname DOPC")
        if not membrane_atoms:
            raise error.MoleculeNotFound(f"Molecule DOPC doesn't exist in system")

        dict_store[i]["universe"] = copied_universe
        dict_store[i]["membrane_atoms"] = membrane_atoms
        dict_store[i]["i"] = frame_int[0]
        dict_store[i]["j"] = frame_int[1]

    return dict_store

def slice_intervals(total_size, slice_size):
    nb_slices = math.ceil(total_size / slice_size)
    intervals = []
    prev = 0
    for i in range(1, nb_slices + 1):
        intervals.append((prev * slice_size, i*slice_size))
        prev = i
    return intervals