#!/usr/bin/env python3

import argparse
import logging
import time
import MDAnalysis as mda
from MDAnalysis.lib import distances
from typing import List, Dict, Set
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level = logging.INFO, format='%(levelname)s\t%(message)s')

def args_gestion():
    parser = argparse.ArgumentParser(description="Script to identify TO clusters size and position along the MD trajectory.")
    parser.add_argument("-f", "--traj", help = "Trajectory (all formats accepted by MDAnalysis)", required = True,  type=str)
    parser.add_argument("-s", "--topo", help = "Topology (all formats accepted by MDAnalysis)", required = True, metavar ="FILE", type=str)
    parser.add_argument("-o", "--outdir", help = "Output directory (default : .)", default=".", metavar = "DIR", type=str)
    parser.add_argument("-p", "--prefix", help = "Output prefix (default : Trajectory file name)", metavar = "STR", type=str)
    parser.add_argument("-t", "--threshold", help = "Threshold for clustering (default : 13)", metavar = "NUMBER", type=float, default=13)
    parser.add_argument("-n", "--to-keep", help = "Number of largest clusters to keep for the results report (default : 2)", metavar = "NUMBER", type=int, default = 2)
    parser.add_argument("--no-plot", help = "Disable plot creation", action="store_true")
    args = parser.parse_args()

    if not args.prefix :
        args.prefix = ".".join(args.traj.split("/")[-1].split(".")[:-1]) # Get trajectory file name without all directories path and without file extension.
    
    return args

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



def serialize_size(out_handle, clusters:List[Set[int]], nb_to_keep:int):
    pass

def add_size_results(time, clusters):
    """Complete RESULTS global variable with size results. Complete "time" list with given time and "size" list with given clusters size (first cluster size will be in first sublist, second in second sublist...)

    Args:
        time ([type]): frame time
        clusters ([type]): clusters to report in results
    """
    RESULTS["time"].append(time)
    for i,clust in enumerate(clusters):
        RESULTS["size"][i].append(len(clust))


def plot_size(out_prefix):
    """Plot clusters size results for the first X clusters (X is given by user), with one line per cluster. 

    Args:
        out_prefix ([type]): The prefix for output image name. 
    """
    for cluster_sizes in RESULTS["size"]:
        plt.plot(RESULTS["time"], cluster_sizes)
        plt.ylim(0,300)
        plt.savefig(f"{out_prefix}_size_check_clusters.png")

def make_correspondance_by_residues(clusters1, clusters2):
    if not clusters1:
        return clusters2[:ARGS.to_keep]
    kept_clusters = []
    for c1 in clusters1:
        print("c1", c1)
        max_percentage = 0
        to_add_cluster = set()
        for c2 in clusters2:
            common_residues_percentage = len(c1.intersection(c2)) / len(c1) 
            print("c2", c2, common_residues_percentage)
            if common_residues_percentage == 1:
                to_add_cluster = c2
                break

            if common_residues_percentage > max_percentage:
                to_add_cluster = c2
                max_percentage = common_residues_percentage

        if not to_add_cluster:
            logging.error(f"Can't assign cluster at frame {SYSTEM.trajectory.frame} ({SYSTEM.trajectory.time} ps) to previous clusters.")
            exit()
        
        kept_clusters.append(to_add_cluster)
    
    return kept_clusters

def make_correspondance_by_position2(clusters_prev, clusters_curr, to_group):
    if not clusters_prev:
        return clusters_curr[:ARGS.to_keep]
    
    kept_clusters = []

    com_prev = np.array([cluster_center_of_mass(c, to_group) for c in clusters_prev])
    com_curr = np.array([cluster_center_of_mass(c, to_group) for c in clusters_curr])
    
    dist_matrix = distances.distance_array(com_prev, com_curr)
    
    for i,dist in enumerate(dist_matrix):
        closest_cluster = np.where(dist == min(dist))[0]

        if len(closest_cluster) > 1:
            logging.warn(f"More than one same center of mass distance for cluster {i} correspondance. Assign the first.")

        if len(clusters_curr[closest_cluster[0]].residues) < 150:
            print("OOOO") 
        kept_clusters.append(clusters_curr[closest_cluster[0]])

    return kept_clusters


def cluster_center_of_mass(cluster, to_group):
    atoms = to_group.residues[list(cluster)].atoms
    return atoms.center_of_mass()

def init_results():
    init_frame = SYSTEM.trajectory[0]
    clusters = compute_clusters()
    results = [Cluster(c,i+1) for i,c in enumerate(clusters[:ARGS.to_keep])]
    return results

def compute_through_time():
    for ts in SYSTEM.trajectory[1:2]:
        to_group = SYSTEM.select_atoms("resname TO")
        clusters = compute_clusters()
        clusters_obj = [Cluster(c) for c in clusters] 
        reordered_clusters = make_correspondance_by_position2(clusters_obj, to_group)
    

class Cluster:
    def __init__(self, residues, idx = None):
        self.idx = idx
        self.residues = residues
        self.frame = SYSTEM.trajectory.frame
        self.frame_time = SYSTEM.trajectory.time

        @property
        def center_of_mass(self):
            if not self._center_of_mass.any():
                cluster_atoms = self.collection.selected_atoms.residues[list(self.residues)].atoms
                self._center_of_mass = cluster_atoms.center_of_mass()
            return self._center_of_mass


def compute_clusters():
    to_group = SYSTEM.select_atoms("resname TO")
    center_of_masses = to_group.center_of_mass(compound="residues") #Center of mass by residues
    neighbors = distances.self_capped_distance(center_of_masses, ARGS.threshold, box = SYSTEM.dimensions) #Compute neighborhood with TO center of masses
    formated_neighbors = format_neighborhood(neighbors[0])
    clusters = clustering(formated_neighbors, to_group.n_residues)
    return clusters


if __name__ == "__main__":
    ARGS = args_gestion()

    start = time.time()

    logging.info("Load system...")

    SYSTEM = mda.Universe(ARGS.topo, ARGS.traj)
    logging.info("Compute on all frames...")

    RESULTS = init_results()
    compute_through_time()
    
    #plot_size(ARGS.outdir + "/" + ARGS.prefix)
    logging.info(f"Analysis complete in {round(time.time() - start, 3)} s")

   


    

