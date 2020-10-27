import logging
from molecules_aggregate.utils import parallelization
from multiprocessing import Pool 
from . import membrane
from . import config
import molecules_aggregate
from . import clusters
from typing import List, Set

class Frame:
    def __init__(self, trj):
        self.frame = trj.frame
        self.time = trj.time
        self.clusters = None
        self.membrane = None
        self.previous = None
    
    def load_clusters(self, atomgroup):
        self.clusters = clusters.load_clusters(atomgroup, self.previous)

    def load_membrane(self, atomgroup, slice_intervals):
        self.membrane = membrane.Membrane(atomgroup, slice_intervals)

    def add_previous(self, frame):
        self.previous = frame
        


class TrajectoryIterator:
    """Object to iterate over trajectory and compute membrane and clusters through each frame. 
    Will load each frame and compute membrane and clusters by using config parameters 
    """
    def __init__(self):
        self.frames = []

    def __iter__(self):
        return iter(self.frames)

    def __getitem__(self, i):
        return self.frames[i]

    def __len__(self):
        return len(self.frames)
    
    def _register_previous(self):
        for i, frame in enumerate(self.frames[1:]): 
            frame.add_previous(self.frames[i])

    def compute_clusters(self):
        """ Compute clusters for each frame. 
        All frames are serial iterated again, to compute clusters. This step is not parallelized because each frame needs clusters from the previous one to keep track of clusters through time. 
        """

        cluster_atomgroup = config.UNIVERSE.select_atoms(f"resname {config.CLUSTER_MOL}")
        self._register_previous()
        logging.info("== Load clusters...")
        for trj in config.UNIVERSE.trajectory[:config.NB_FRAMES]:
            if trj.frame % 1000 == 0:
                logging.info(f"Frame {trj.frame}/{config.NB_FRAMES}")
      
            frame = self.frames[trj.frame]

            frame.load_clusters(cluster_atomgroup)


    def compute_membranes(self, create_frames = False):
        """Compute membrane for each frame. 
        Step 1 : initiate parallelization. One copy of the universe per process is created, to avoid iterate concurrently on universe, which is impossible. Frames are separated into chunks, same number as process. 
        Step 2 : compute membrane with process. For each frames chunks, frames are serial iterated, and membrane is sliced to get highest and lowest point according to given axis. 
        """
        #Init parallelization
        logging.info(f"Copy universe for parallelization with {config.NB_PROCESS} process...")
        frames_chunks = parallelization.slice_frames(config.NB_FRAMES, config.NB_PROCESS)
        global COPIED_UNIVERSES
        COPIED_UNIVERSES = parallelization.copy_universe(frames_chunks, config.UNIVERSE, config.MEMBRANE_MOL)
        
        #Load membrane, with several process
        logging.info(f"Load membrane...") 
        slice_intervals = membrane.get_slice_intervals(config.UNIVERSE.dimensions[config.SLICE_AXIS_IDX], config.SLICE_SIZE)
        p = Pool(config.NB_PROCESS)
        packed_frames = p.starmap(self._some_frames_membrane, [(universe_idx, slice_intervals) for universe_idx in COPIED_UNIVERSES])
        p.close()
        p.join()
        frames = [frame for process in packed_frames for frame in process]


        for i, frame in enumerate(frames[1:]): 
            frame.add_previous(frames[i])

        self.frames = frames


    def _some_frames_membrane(self, universe_idx:int, slice_intervals:List[Set[int]]) -> List[Frame]:
        """Compute membrane for a chunk of frames
        Frame objects in self.frames are some kind of copy because this function is called with starmap, so we need to restore this objects into a new list that we return. 

        Args:
            universe_idx (int): index of the copied universe to work with
            slice_intervals (List[Set[int]]): interval for membrane slicing

        Returns:
            List[Frame]: List of computed frames
        """
        new_frames = [] #For some reason frames object are not the same as initial frame, maybe because it's called inside starmap
        universe = COPIED_UNIVERSES[universe_idx]["universe"]
        atomgroup = COPIED_UNIVERSES[universe_idx]["atomgroup"]
        i = COPIED_UNIVERSES[universe_idx]["i"]
        j = COPIED_UNIVERSES[universe_idx]["j"]
        for trj in universe.trajectory[i:j]: 
            frame = Frame(trj)
            frame.load_membrane(atomgroup, slice_intervals)
            new_frames.append(frame)

        return new_frames

        




