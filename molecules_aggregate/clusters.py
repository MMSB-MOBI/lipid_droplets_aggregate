from . import config
from MDAnalysis.lib import distances
from typing import List, Dict, Set
import numpy as np
import logging

class Cluster:
    def __init__(self, residues, atomgroup, idx = None):
        self.residues = residues
        self.idx = idx
        self.center_of_mass = self._center_of_mass(atomgroup)

    def _center_of_mass(self, atomgroup):
        atoms = atomgroup.residues[list(self.residues)].atoms
        return atoms.center_of_mass()

    @property
    def size(self):
        return len(self.residues)


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

def compute_clusters(box_dimension, mol_group, threshold):
    """Compute clusters.
    Center of mass is computed for each TO (or other molecule if other is given) residue. Then neighborhood is computed with MDAnalysis function and this neighborhood is reformated to have direct access to the list of neighbors for each residue, and clustering is done from these lists with clustering() function.  

    Args:
        box_dimension (List[float]): MDAnalysis universe box dimension
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

def make_correspondance_by_residues(new_clusters_residues, previous_frame, atomgroup):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with residue method. 
        For each cluster at time t-1, a common percentage is computed between current clusters. This common percentage is the number of common residues over length of the largest cluster. The cluster with highest common percentage is assigned to corresponding t-1 cluster. 

        Args:
            new_clusters_residues (List[Set[int]]): Current clusters residue numbers 
            previous_frame (molecules_aggregate.trajectory.Frame): molecules_aggregate previous frame
            atomgroup (MDAnalysis.core.groups.AtomGroup) AtomGroup that have been clustered
        """
        previous_clusters = [c.residues for c in previous_frame.clusters]
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
                raise Exception("No closest cluster found")

            closest_cluster = Cluster(new_clusters_residues[closest_cluster_idx], atomgroup, i) 
            reordered_new_clusters.append(closest_cluster)

        return reordered_new_clusters

def make_correspondance_by_position(self, new_clusters_residues, previous_frame, atomgroup):
        """Compute correspondance between current clusters at time t and previous clusters at time t -1 with position method. 
        For each cluster at time t-1, the distance (based on center of mass) with current clusters is computed with MDAnalysis functions. The cluster with closest distance is assigned to corresponding t-1 cluster.

        Args:
            new_clusters_residues (List[Set[int]]): Current clusters residue numbers 
            previous_frame (molecules_aggregate.trajectory.Frame): molecules_aggregate previous frame
            atomgroup (MDAnalysis.core.groups.AtomGroup) AtomGroup that have been clustered
        """
        com_prev = np.array([c.center_of_mass for c in previous_frame.clusters])
        
        #Create tmp clusters to have center of mass calculation 
        tmp_clusters = [Cluster(residues, atomgroup) for residues in new_clusters_residues]
        com_curr = np.array([c.center_of_mass for c in tmp_clusters])

        dist_matrix = distances.distance_array(com_prev, com_curr)

        new_clusters = []
        for i, dist in enumerate(dist_matrix):
            closest_cluster_idx = np.where(dist == min(dist))[0]
            if len(closest_cluster_idx) > 1:
                logging.warn(f"Frame {previous_frame.frame + 1}. More than one same center of mass distance for cluster {i} correspondance. Assign the first.")
    
            closest_cluster = tmp_clusters[closest_cluster_idx[0]]
            closest_cluster.idx = i

            new_clusters.append(closest_cluster)

        return new_clusters


def load_clusters(atomgroup, previous_frame):
    """Compute clusters for current UNIVERSE frame

    Args:
        atomgroup (MDAnalysis.core.groups.AtomGroup): AtomGroup that needs to be clustered
        previous_frame (molecules_aggregate.trajectory.Frame): previous frame

    Returns:
        List[molecules_aggregate.clusters.Cluster]: List of clusters for the current frame
    """
    clusters = compute_clusters(config.UNIVERSE.dimensions, atomgroup, config.CLUSTER_THRESHOLD)
    clusters = clusters if config.NB_CORRESPONDANCE == -1 else clusters[:config.NB_CORRESPONDANCE]
    if previous_frame:
        if config.CORRESPONDANCE == "residue":
            clusters_obj = make_correspondance_by_residues(clusters, previous_frame, atomgroup)
        elif config.CORRESPONDANCE == "position":
            clusters_obj = make_correspondance_by_position(clusters, previous_frame, atomgroup)
        else:
            raise Exception(f"invalid correspondance method {config.CORRESPONDANCE}")

    else: 
        clusters_obj = [Cluster(residues, atomgroup, i) for i,residues in enumerate(clusters[:config.NB_CLUSTERS])]

    return clusters_obj