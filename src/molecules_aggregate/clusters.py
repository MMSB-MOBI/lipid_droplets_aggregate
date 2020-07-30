import logging
from MDAnalysis.lib import distances
from typing import List, Dict, Set
import numpy as np
import copy
import asyncio

class Cluster_bak:
    """A class that stores cluster for 1 frame"""

    def __init__(self, system, residues, atom_group, idx = None, com = True):
        self.system = system
        """MDAnalysis system"""
        self.idx :int = idx
        """Cluster index, 0 is for the first largest etc..."""
        self.residues : Set[int] = residues
        """Cluster residues indexes"""
        self.frame : int = system.trajectory.frame
        """Cluster frame through trajectory"""
        self.time : float = system.trajectory.time
        """Cluster frame time"""
        self.atom_group : MDAnalysis.core.groups.AtomGroup = atom_group
        """The atom group from which cluster has been computed"""
        self.center_of_mass : Optional[numpy.ndarray[float]] = self._center_of_mass() if com else None
        """Cluster center of mass"""

    def _center_of_mass(self):
        """Compute cluster center of mass with MDAnalysis functions

        Returns:
            np.array[float]: Cluster center of mass
        """
        atoms = self.atom_group.residues[list(self.residues)].atoms
        return atoms.center_of_mass()

    @property
    def size(self):
        return len(self.residues)

class Cluster:
    def __init__(self, parent_system, residues, idx = None):
        self.parent_system = parent_system
        self.residues = residues
        self.idx = idx
        self.center_of_mass = self._center_of_mass()

    def _center_of_mass(self):
        atoms = self.parent_system.cluster_atom_group[list(self.residues)].atoms
        return atoms.center_of_mass()

    @property
    def z_position_relative_to_membrane(self):  
        return self.parent_system.membrane.highest_z_mean - self.center_of_mass[2]

    @property
    def size(self):
        return len(self.residues)
        

class OneFrameClusters:
    def __init__(self, parent_system):
        self.parent_system = parent_system
        self.clusters_list = self._load_clusters()

    @property
    def frame(self):
        return self.parent_system.frame
        
    def _load_clusters(self):
        clusters = compute_clusters(self.parent_system.box_dimension, self.parent_system.cluster_atom_group, self.parent_system.cluster_threshold)

        if self.parent_system.previous: 
            clusters_obj = self.make_correspondance_by_position(clusters)
        else:
            clusters_obj = [Cluster(self.parent_system, residues, i) for i,residues in enumerate(clusters[:self.parent_system.to_keep_clusters])]
        
        return clusters_obj

    #def __repr__(self):
    #    return str(self.clusters_list)    

    #def __iter__(self):
    #    return iter(self.clusters_list)

    def __getitem__(self, i):
        return self.clusters_list[i]

    def __len__(self):
        return len(self.clusters_list)

    def make_correspondance_by_residues(self, new_clusters_residues):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with residue method. 
        For each cluster at time t-1, a common percentage is computed between current clusters. This common percentage is the number of common residues over length of the largest cluster. The cluster with highest common percentage is assigned to corresponding t-1 cluster. 

        Args:
            new_clusters (List[List[Cluster]]): Current clusters
            frame_nb (int): Current frame
            system (MDAnalysis.universe.Universe): MDAnalysis universe
        """

        previous_clusters = [c.residues for c in self.parent_system.previous.clusters]
        reordered_new_clusters = []
        for i,res_list in enumerate(previous_clusters):
            max_perc = 0
            closest_cluster_idx = -1
            for j, new_res_list in enumerate(new_clusters_residues):
                common_perc = len(res_list.intersection(new_res_list)) / max(len(res_list), len(new_res_list))
                if common_perc > max_perc : 
                    max_perc = common_perc 
                    closest_cluster_idx = j
            if closest_cluster_idx == -1:
                logging.error("No closest cluster found")
                exit()

            closest_cluster = Cluster(self.parent_system, new_clusters_residues[closest_cluster_idx], i) 
            reordered_new_clusters.append(closest_cluster)

        return reordered_new_clusters

    def make_correspondance_by_position(self, new_clusters_residues):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with position method. 
        For each cluster at time t-1, the distance (based on center of mass) with current clusters is computed with MDAnalysis functions. The cluster with closest distance is assigned to corresponding t-1 cluster.

        Args:
            new_clusters (List[List[Cluster]]): Current clusters
            frame_nb (int): Current frame
            system (MDAnalysis.universe.Universe): MDAnalysis universe
        """
        com_prev = np.array([c.center_of_mass for c in self.parent_system.previous.clusters])
        
        #Create tmp clusters to have center of mass calculation 
        tmp_clusters = [Cluster(self.parent_system, residues) for residues in new_clusters_residues]
        com_curr = np.array([c.center_of_mass for c in tmp_clusters])

        print(tmp_clusters[207].size)
        print(com_prev)
        print(com_curr[207])

        dist_matrix = distances.distance_array(com_prev, com_curr)

        #print(len(dist_matrix))
        
        new_clusters = []
        #assigned = set()
        for i, dist in enumerate(dist_matrix):
            #print(i, dist)
            closest_cluster_idx = np.where(dist == min(dist))[0]
            #print(i, closest_cluster_idx)
            #print(dist[207])
            if len(closest_cluster_idx) > 1:
                logging.warn(f"Frame {frame_nb}. More than one same center of mass distance for cluster {i} correspondance. Assign the first.")
    
            closest_cluster = tmp_clusters[closest_cluster_idx[0]]
            closest_cluster.idx = i

            new_clusters.append(closest_cluster)

        return new_clusters

            



class ClusterIterator:

    def __init__(self, mda_system, molecule, cluster_threshold, nb_to_keep, frames_to_process, cluster_correspondance, nb_cluster_for_corr):
        """Init function

        Args:
            mda_system (MDAnalysis.universe.Universe): MDAnalysis universe
            molecule (str): Molecule name we want to cluster
            cluster_threshold (float): Distance below which a residue is assigned to cluster
            nb_to_keep (int): Number of clusters to keep for results
            frames_to_process (int): Number of frames to process for results
            cluster_correspondance (residue|position): Cluster correspondance method to assign clusters at time t to clusters at time t - 1. "residue" is based on residue comparison and "position" is based on position comparison. 
            nb_cluster_for_corr (int): Number of largest clusters to check at time t for correspondance with clusters at time t - 1. 
        """
        self.molecule:str = molecule
        """Molecule name"""
        self.mol_group : MDAnalysis.core.groups.AtomGroup = mda_system.select_atoms(f"resname {molecule}")
        """MDAnalysis AtomGroup for given molecule"""
        self.cluster_threshold:int = cluster_threshold
        """Clustering threshold"""
        self.clusters : List[List[Cluster]] = self.init_clusters(mda_system, nb_to_keep)
        self.clusters_through_time(mda_system, frames_to_process, cluster_correspondance, nb_cluster_for_corr)

    @property
    def nb_frame(self):
        """Number of frames"""
        return len(self.clusters[0])

    @property
    def nb_cluster(self):
        """Number of clusters"""
        return len(self.clusters)

    
    def __iter__(self):
        for clusters in self.clusters:
            yield clusters

    def __getitem__(self,i):
        return self.clusters[i]

    def init_clusters(self, system, nb_to_keep):
        """Initialize clusters at frame 0

        Args:
            system (MDAnalysis.universe.Universe): MDAnalysis universe
            nb_to_keep (int): Number of clusters to keep for results

        Returns:
            List[List[Cluster]]: The clusters. The list is organized this way : each element is for one cluster and is a list that contains Cluster object. Each object is cluster for 1 frame. 
        """
        system.trajectory[0]
        if not self.mol_group:
            logging.warn(f"No molecules {molecule_name} in system")
            exit()
            
        clusters = compute_clusters(system, self.mol_group, self.cluster_threshold)
        return [[Cluster(system, res, self.mol_group, idx)] for idx,res in enumerate(clusters[:nb_to_keep])]

    def clusters_through_time(self, system, frames_to_process:int, cluster_correspondance:str, nb_corr:int):
        """Compute clusters after initialization, from frame 1 to the end. Complete self.clusters attribute. 

        Args:
            system (MDAnalysis.universe.Universe): MDAnalysis universe
            frames_to_process (int): Number of frames to process
            cluster_correspondance (str : residue|position): Method for cluster correspondance between time t and t-1
            nb_corr (int): Number of clusters to look at for cluster correspondance
        """
        if frames_to_process == -1:
            frames = system.trajectory[1:]
        else:
            frames = system.trajectory[1:frames_to_process + 1]
        
        com_computation = True if cluster_correspondance == "position" else False

        cluster_for_corr = '' if nb_corr == -1 else nb_corr

        nb_frames = len(frames)
        for ts in frames:
            if ts.frame % 1000 == 0:
                logging.info(f"Frame {ts.frame}/{nb_frames}")
            clusters = compute_clusters(system, self.mol_group, self.cluster_threshold)
            
            tmp_clusters = [Cluster(system, res, self.mol_group, idx = None, com = com_computation) for res in clusters] if nb_corr == -1 else [Cluster(system, res, self.mol_group, idx = None, com = com_computation) for res in clusters[:nb_corr]]
            
            if cluster_correspondance == "position":
                self.make_correspondance_by_position(tmp_clusters, ts.frame, system)
            elif cluster_correspondance == "residue":
                self.make_correspondance_by_residues(tmp_clusters, ts.frame, system)
        
    def make_correspondance_by_residues(self, new_clusters, frame_nb:int, system):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with residue method. 
        For each cluster at time t-1, a common percentage is computed between current clusters. This common percentage is the number of common residues over length of the largest cluster. The cluster with highest common percentage is assigned to corresponding t-1 cluster. 

        Args:
            new_clusters (List[List[Cluster]]): Current clusters
            frame_nb (int): Current frame
            system (MDAnalysis.universe.Universe): MDAnalysis universe
        """

        previous_clusters = [c[frame_nb - 1] for c in self.clusters]
        for i,res_list in enumerate([c.residues for c in previous_clusters]):
            max_perc = 0
            closest_cluster_idx = -1
            for j, new_res_list in enumerate([c.residues for c in new_clusters]):
                common_perc = len(res_list.intersection(new_res_list)) / max(len(res_list), len(new_res_list))
                if common_perc > max_perc : 
                    max_perc = common_perc 
                    closest_cluster_idx = j
            if closest_cluster_idx == -1:
                logging.error("No closest cluster found")
                exit()

            closest_cluster = new_clusters[closest_cluster_idx]
            self.clusters[i].append(Cluster(system, closest_cluster.residues, closest_cluster.atom_group, i))

    def make_correspondance_by_position(self, new_clusters, frame_nb, system):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with position method. 
        For each cluster at time t-1, the distance (based on center of mass) with current clusters is computed with MDAnalysis functions. The cluster with closest distance is assigned to corresponding t-1 cluster.

        Args:
            new_clusters (List[List[Cluster]]): Current clusters
            frame_nb (int): Current frame
            system (MDAnalysis.universe.Universe): MDAnalysis universe
        """
        com_prev = np.array([c[frame_nb - 1].center_of_mass for c in self.clusters])
        com_curr = np.array([c.center_of_mass for c in new_clusters])
        dist_matrix = distances.distance_array(com_prev, com_curr)
        
        #assigned = set()
        for i, dist in enumerate(dist_matrix):
            closest_cluster_idx = np.where(dist == min(dist))[0]
            if len(closest_cluster_idx) > 1:
                logging.warn(f"Frame {frame_nb}. More than one same center of mass distance for cluster {i} correspondance. Assign the first.")
    
            #assigned.add(closest_cluster_idx)
            closest_cluster = new_clusters[closest_cluster_idx[0]]

            closest_cluster.idx = i
            self.clusters[i].append(closest_cluster)

def compute_clusters(box_dimension, mol_group, threshold):
    """Compute clusters.
    Center of mass is computed for each TO (or other molecule if other is given) residue. Then neighborhood is computed with MDAnalysis function and this neighborhood is reformated to have direct access to the list of neighbors for each residue, and clustering is done from these lists with clustering() function.  

    Args:
        system (MDAnalysis.universe.Universe): MDAnalysis universe
        mol_group (MDAnalysis.core.groups.AtomGroup): MDAnalysis AtomGroup for residues to cluster
        threshold (int): Clustering threshold

    Returns:
        List[Set[int]]: The list of all clusters, each element of the list is a set containing the residues indexes that are in the same cluster. The list is ordered in descending order with largest cluster in first position. 
    """
    center_of_masses = mol_group.center_of_mass(compound="residues") #Center of mass by residues
    neighbors = distances.self_capped_distance(center_of_masses, threshold, box = box_dimension) #Compute neighborhood with TO center of masses
    formated_neighbors = format_neighborhood(neighbors[0])
    clusters = clustering(formated_neighbors, mol_group.n_residues)
    return clusters


def format_neighborhood(neighbors: List[List[int]]) -> Dict[int, List[int]]:
    """Take neighborhood computed by MDAnalysis and reorganize it to have direct access to the complete list of neighbors for each residue. 

    Args:
        neighbors (List[List[int]]): List of neighbors pairs. Each element of the list is a list containing 2 neighbors residues index. 

    Returns:
        Dict[int, List[int]]: Dictionnary that stores for each residue all its neighbors. Keys are residue index and values are the list of neighbors residues indexes. 
    """
    neighbors_list = {}
    for n in neighbors:
        if n[0] not in neighbors_list:
            neighbors_list[n[0]] = []
        if n[1] not in neighbors_list:
            neighbors_list[n[1]] = []
        neighbors_list[n[0]].append(n[1])
        neighbors_list[n[1]].append(n[0])
    return neighbors_list


def clustering(formated_neighbors:Dict[int, List[int]], nb_residues:int) -> List[Set[int]]:
    """Recursive clustering procedure inspired by the g_aggregate C code. It applies the same principle that consists of identify connected component on a graph. 

    Args:
        formated_neighbors (Dict[int, List[int]]): The dictionnary that represents neighborhood. Given by format_neighborhood() function.
        nb_residues (int): Number of residues to cluster

    Returns:
        List[Set[int]]: The list of all clusters, each element of the list is a set containing the residues indexes that are in the same cluster. The list is ordered in descending order with largest cluster in first position. 
    """

    def recursive_add_to_cluster(cluster, residue, assigned, neighbors_list):
        neighbors = neighbors_list.get(residue,[])
        cluster.add(residue)
        assigned.add(residue)
        for res in neighbors:
            if not res in assigned:
                cluster, assigned = recursive_add_to_cluster(cluster, res, assigned, neighbors_list)
            
        return cluster, assigned

    clusters_array = [] 
    assigned = set() #Clusters already assigned
    for res in range(nb_residues):
        if not res in assigned:
            cluster, assigned = recursive_add_to_cluster(set(), res, assigned, formated_neighbors)
            clusters_array.append(cluster)
    
    sorted_clusters = sorted(clusters_array, key = lambda cluster:len(cluster), reverse = True)
    return sorted_clusters