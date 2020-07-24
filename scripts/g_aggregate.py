#!/usr/bin/env python3

import argparse
import MDAnalysis as mda
import molecules_aggregate
import time
import logging
logging.basicConfig(level = logging.INFO, format='%(levelname)s\t%(message)s')
import cProfile, pstats, io
from pstats import SortKey
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def args_gestion():
    parser = argparse.ArgumentParser(description="Script to identify TO clusters size and position along the MD trajectory.")
    parser.add_argument("-f", "--traj", help = "Trajectory (all formats accepted by MDAnalysis)", required = True,  type=str, metavar = "FILE")
    parser.add_argument("-s", "--topo", help = "Topology (all formats accepted by MDAnalysis)", required = True, metavar ="FILE", type=str)
    parser.add_argument("-o", "--outdir", help = "Output directory (default : .)", default=".", metavar = "DIR", type=str)
    parser.add_argument("-p", "--prefix", help = "Output prefix (default : Trajectory file name)", metavar = "STR", type=str)
    parser.add_argument("-t", "--threshold", help = "Threshold for clustering (default : 13)", metavar = "NUMBER", type=float, default=13)
    parser.add_argument("-n", "--to-keep", help = "Number of largest clusters to keep for the results report (default : 2)", metavar = "NUMBER", type=int, default = 2)
    parser.add_argument("-m", "--method", help = "Method for clusters correspondance through frames (residue or position)", type=str, default = "residue", choices=['residue', 'position'])
    parser.add_argument("-c", "--nb-corr", help = "Number of largest clusters to take into account for clusters correspondance at time t (default : number of clusters to keep, -1 for all)", type = int, metavar = "NUMBER")
    parser.add_argument("--frames", help = "Number of frames to process (default : all)", metavar = "NUMBER", type = int, default=0)

    args = parser.parse_args()

    if not args.prefix :
        args.prefix = ".".join(args.traj.split("/")[-1].split(".")[:-1]) # Get trajectory file name without all directories path and without file extension.

    if not args.nb_corr :
        args.nb_corr = args.to_keep

    if args.nb_corr != -1 and args.nb_corr < args.to_keep : 
        logging.error("Number of clusters for correspondance calculation (-c/--nb-corr) can't be lower than number of clusters to keep (-n/--to-keep).")
        exit()
    
    return args

def plot_size(out_prefix:str, clusters:molecules_aggregate.clusters.ClusterIterator):
    """Plot size of clusters provided in ClusterIterator through time

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """

    fig, ax = plt.subplots()
    ax.set_title(f"Size of the {clusters.nb_cluster} largest clusters through time ({ARGS.method} for cluster correspondance)")
    ax.set_ylim(0,300)
    ax.set_xlabel("Time")
    ax.set_ylabel("Cluster size")

    for cluster_frames in clusters:
        ax.plot([c.time for c in cluster_frames], [c.size for c in cluster_frames], label = f"Cluster {cluster_frames[0].idx}")
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_size_clusters.png")

def plot_com(out_prefix:str, clusters:molecules_aggregate.clusters.ClusterIterator):
    """Plot absolute center of masses of clusters for all frames along x,y and z axis on "3d" plot. 

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f"Absolute center of mass of clusters ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0, system.dimensions[0])
    ax.set_ylim(0, system.dimensions[1])
    ax.set_zlim(0, system.dimensions[2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[0] for c in cluster_frames], [c.center_of_mass[1] for c in cluster_frames], [c.center_of_mass[2] for c in cluster_frames], label =  f"Cluster {cluster_frames[0].idx}", alpha = 0.1, s=1, depthshade = False)

    leg = ax.legend(loc="center left")
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([10])
    fig.savefig(f"{out_prefix}_absolute_com.png")  

def plot_2d(out_prefix:str, clusters:molecules_aggregate.clusters.ClusterIterator):
    """Plot absolute center of masses of clusters for all frames with 2d projection along xz, yx and yz axis. 

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """
    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through x/z axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[0])
    ax.set_ylim(0,system.dimensions[2])
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[0] for c in cluster_frames], [c.center_of_mass[2] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}", alpha = 0.1)
    
    leg = ax.legend()
    [lh.set_alpha(1) for lh in leg.legendHandles]
    fig.savefig(f"{out_prefix}_xz.png")

    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through y/x axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[1])
    ax.set_ylim(0,system.dimensions[0])
    ax.set_xlabel("y")
    ax.set_ylabel("x")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[1] for c in cluster_frames], [c.center_of_mass[0] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}", alpha = 0.1)
    leg = ax.legend()
    [lh.set_alpha(1) for lh in leg.legendHandles]
    fig.savefig(f"{out_prefix}_yx.png")

    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through y/z axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[1])
    ax.set_ylim(0,system.dimensions[2])
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[1] for c in cluster_frames], [c.center_of_mass[2] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}", alpha = 0.1)
    leg = ax.legend()
    [lh.set_alpha(1) for lh in leg.legendHandles]
    fig.savefig(f"{out_prefix}_yz.png")

def plot_z_membrane(out_prefix, membrane):
    fig, ax = plt.subplots()
    #ax.set_ylim(0, system.dimensions[2])
    ax.scatter([i for i in range(len(membrane.slices))],membrane.slices)

    plt.savefig(f"{out_prefix}_z_membrane.png")


if __name__ == "__main__":
    ARGS = args_gestion()

    start = time.time()

    logging.info("Load system...")

    system = mda.Universe(ARGS.topo, ARGS.traj)
    logging.info("Compute clusters through all frames...")

    #This lines are for profiling the time
    '''pr = cProfile.Profile()
    pr.enable()'''

    clusters = molecules_aggregate.load_clusters(system, "TO", ARGS.threshold, ARGS.to_keep, ARGS.frames, ARGS.method, ARGS.nb_corr)

    #This lines are for profiling the time
    '''pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())'''

    plot_size(ARGS.prefix, clusters)
    plot_com(ARGS.prefix, clusters)
    plot_2d(ARGS.prefix, clusters)
    logging.info(f"Analysis end in {time.time() - start} s")





