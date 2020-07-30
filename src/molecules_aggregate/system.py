from . import clusters as clusters_lib
from . import membrane as membrane_lib
from . import error
import logging

class System:
    def __init__(self, parent, mda_system, previous = None):
        self.parent = parent
        self.frame = mda_system.trajectory.frame
        self.time = mda_system.trajectory.time
        self.previous = previous
        self.clusters = self._load_clusters()
        #self.membrane = self._load_membrane()
        

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
    def membrane_atom_group(self):
        return self.parent.membrane_atom_group

    @property
    def to_keep_clusters(self):
        return self.parent.to_keep_clusters
        
    @property
    def mda_system(self):
        return self.parent.mda_system

    @property
    def nb_cluster(self):
        return len(self.clusters)

    def _load_clusters(self):
        return clusters_lib.OneFrameClusters(self)

    def _load_membrane(self):
        return membrane_lib.Membrane(self.membrane_atom_group, self.parent.slice_membrane_size, self.parent.slice_axis, self.parent.box_dimension)

    

class SystemIterator:
    def __init__(self, mda_system, cluster_mol_name, cluster_threshold, to_keep_clusters, membrane_mol_name, slice_membrane_size, nb_frames, slice_axis = 'x'):
        self.cluster_atom_group = mda_system.select_atoms(f"resname {cluster_mol_name}")
        self.membrane_atom_group = mda_system.select_atoms(f"resname {membrane_mol_name}")
        if not self.cluster_atom_group:
            raise error.MoleculeNotFound(f"Molecule {cluster_mol_name} doesn't exist in system")
        if not self.membrane_atom_group:
            raise error.MoleculeNotFound(f"Molecule {membrane_mol_name} doesn't exist in system")

        self.box_dimension = mda_system.dimensions
        self.cluster_threshold = cluster_threshold
        self.to_keep_clusters = to_keep_clusters
        self.slice_membrane_size = slice_membrane_size
        self.slice_axis = slice_axis
        self.mda_system = mda_system

        self.system_frames = self._load_frames(mda_system, nb_frames)

    @property 
    def nb_clusters(self):
        return self.system_frames[0].nb_cluster

    def __iter__(self):
        return iter(self.system_frames)

    def __getitem__(self,i):
        return self.system_frames[i]

    def _load_frames(self, mda_system, nb_frames):
        mda_system.trajectory[0]
        logging.info("Init first frame")
        frames = [System(self, mda_system)]

        if nb_frames == -1:
            traj = mda_system.trajectory[1:]
        else:
            traj = mda_system.trajectory[1:nb_frames]

        logging.info("Load other frames")
        nb_frames = len(traj)
        for frame in traj:
            if frame.frame % 1000 == 0:
                logging.info(f"Frame {frame.frame}/{nb_frames}")
            frames.append(System(self, mda_system, frames[frame.frame - 1]))

        return frames


