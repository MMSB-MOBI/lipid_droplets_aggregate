#!/usr/bin/env python3

import argparse
import logging
import time
import MDAnalysis as mda
from MDAnalysis.lib import distances
from typing import List, Dict, Set
import matplotlib.pyplot as plt

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
    clusters_array = [] 
    assigned = set() #Clusters already assigned
    for res in range(nb_residues):
        if not res in assigned:
            cluster, assigned = recursive_add_to_cluster(set(), res, assigned, formated_neighbors)
            clusters_array.append(cluster)
    
    sorted_clusters = sorted(clusters_array, key = lambda cluster:len(cluster), reverse = True)
    return sorted_clusters

def recursive_add_to_cluster(cluster, residue, assigned, neighbors_list):
        neighbors = neighbors_list.get(residue,[])
        cluster.add(residue)
        assigned.add(residue)
        for res in neighbors:
            if not res in assigned:
                cluster, assigned = recursive_add_to_cluster(cluster, res, assigned, neighbors_list)
            
        return cluster, assigned

def serialize_size(out_handle, clusters:List[Set[int]], nb_to_keep:int):
    pass



def add_size_results(time, clusters):
    RESULTS["time"].append(time)
    for cluster_nb in range(ARGS.to_keep):
        RESULTS[f"cluster{cluster_nb + 1}_size"].append(len(clusters[cluster_nb]))

def plot_size():
    for nb_cluster in range(ARGS.to_keep):
        plt.plot(RESULTS["time"], RESULTS[f"cluster{nb_cluster + 1}_size"])
        plt.ylim(0,300)
        plt.savefig("test_size.png")


if __name__ == "__main__":
    ARGS = args_gestion()

    start = time.time()

    logging.info("Load system...")

    system = mda.Universe(ARGS.topo, ARGS.traj)

    logging.info("Compute on all frames...")

    nb_frames = len(system.trajectory)

    RESULTS = {"time": []}

    for cluster_nb in range(ARGS.to_keep):
        RESULTS[f"cluster{cluster_nb+1}_size"] = []

    for ts in system.trajectory[:1000]:
        if ts.frame % 1000 == 0:
            logging.info(f"Frame {ts.frame}/{nb_frames}")
        to_group = system.select_atoms("resname TO")
        center_of_masses = to_group.center_of_mass(compound="residues") #Center of mass by residues
        neighbors = distances.self_capped_distance(center_of_masses, ARGS.threshold, box = system.dimensions) #Compute neighborhood with TO center of masses
        formated_neighbors = format_neighborhood(neighbors[0])
        clusters = clustering(formated_neighbors, to_group.n_residues)
        #print(clusters[0])
        add_size_results(ts.time, clusters)
        #size_data[ts.time] = 
        #serialize_size(ARGS.outdir + "/" + ARGS.prefix, clusters, ARGS.to-keep)

    plot_size()
    logging.info(f"Analysis complete in {round(time.time() - start, 3)} s")

   


    

