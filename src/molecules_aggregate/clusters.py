import logging
from MDAnalysis.lib import distances
from typing import List, Dict, Set
import numpy as np
import copy

class Cluster:
    def __init__(self, system, residues, atom_group, idx = None, com = True):
        #print("com", com)
        self.system = system
        self.idx = idx
        self.residues = residues
        self.frame = system.trajectory.frame
        self.time = system.trajectory.time
        self.atom_group = atom_group
        self.center_of_mass = self._center_of_mass(atom_group) if com else None

    def _center_of_mass(self, atom_group):
        #print("COM", self.system.trajectory.frame)
        #if not self._center_of_mass.any():
        atoms = atom_group.residues[list(self.residues)].atoms
        return atoms.center_of_mass()

    @property
    def size(self):
        return len(self.residues)


class ClusterIterator:
    def __init__(self, mda_system, molecule, cluster_threshold, nb_to_keep, frames_to_process, cluster_correspondance):
        self.molecule = molecule
        self.mol_group = mda_system.select_atoms(f"resname {molecule}")
        self.cluster_threshold = cluster_threshold
        self.clusters = self.init_clusters(mda_system, nb_to_keep)
        self.clusters_through_time(mda_system, frames_to_process, nb_to_keep, cluster_correspondance)

    @property
    def nb_frame(self):
        return len(self.clusters[0])

    @property
    def nb_cluster(self):
        return len(self.clusters)


    def __iter__(self):
        for clusters in self.clusters:
            yield clusters

    def __getitem__(self,i):
        return self.clusters[i]

    def init_clusters(self, system, nb_to_keep):
        system.trajectory[0]
        if not self.mol_group:
            logging.warn(f"No molecules {molecule_name} in system")
            exit()
            
        clusters = compute_clusters(system, self.mol_group, self.cluster_threshold)
        return [[Cluster(system, res, self.mol_group, idx)] for idx,res in enumerate(clusters[:nb_to_keep])]

    def clusters_through_time(self, system, frames_to_process, nb_to_keep, cluster_correspondance):
        if frames_to_process == 0:
            frames = system.trajectory[1:]
        else:
            frames = system.trajectory[1:frames_to_process + 1]
        
        com_computation = True if cluster_correspondance == "position" else False
        nb_frames = len(frames)
        for ts in frames:
            if ts.frame % 1000 == 0:
                logging.info(f"Frame {ts.frame}/{nb_frames}")
            clusters = compute_clusters(system, self.mol_group, self.cluster_threshold)
            tmp_clusters = [Cluster(system, res, self.mol_group, idx = None, com = com_computation) for res in clusters]
            if cluster_correspondance == "position":
                self.make_correspondance_by_position(tmp_clusters, ts.frame, system)
            elif cluster_correspondance == "residue":
                self.make_correspondance_by_residues(tmp_clusters, ts.frame, system)
        
    def make_correspondance_by_residues(self, new_clusters, frame_nb, system):
        for i,res_list in enumerate([c[frame_nb - 1].residues for c in self.clusters]):
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
        com_prev = np.array([c[frame_nb - 1].center_of_mass for c in self.clusters])
        com_curr = np.array([c.center_of_mass for c in new_clusters])

        dist_matrix = distances.distance_array(com_prev, com_curr)

        for i, dist in enumerate(dist_matrix):
            closest_cluster_idx = np.where(dist == min(dist))[0]
            if len(closest_cluster_idx) > 1:
                logging.warn(f"Frame {frame_nb}. More than one same center of mass distance for cluster {i} correspondance. Assign the first.")
    
            closest_cluster = new_clusters[closest_cluster_idx[0]]

            closest_cluster.idx = i
            self.clusters[i].append(closest_cluster)

def compute_clusters(system, mol_group, threshold):
    center_of_masses = mol_group.center_of_mass(compound="residues") #Center of mass by residues
    neighbors = distances.self_capped_distance(center_of_masses, threshold, box = system.dimensions) #Compute neighborhood with TO center of masses
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