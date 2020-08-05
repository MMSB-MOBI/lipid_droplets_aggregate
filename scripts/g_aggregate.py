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
import os


def args_gestion():
    parser = argparse.ArgumentParser(description="Script to identify TO clusters size and position along the MD trajectory.")
    parser.add_argument("-f", "--traj", help = "Trajectory (all formats accepted by MDAnalysis)", required = True,  type=str, metavar = "FILE")
    parser.add_argument("-s", "--topo", help = "Topology (all formats accepted by MDAnalysis)", required = True, metavar ="FILE", type=str)
    parser.add_argument("-o", "--outdir", help = "Output directory (default : .)", default=".", metavar = "DIR", type=str)
    parser.add_argument("-p", "--prefix", help = "Output prefix (default : Trajectory file name)", metavar = "STR", type=str)
    parser.add_argument("-t", "--threshold", help = "Threshold for clustering (default : 13)", metavar = "NUMBER", type=float, default=13)
    parser.add_argument("-n", "--to-keep", help = "Number of largest clusters to keep for the results report (default : 2)", metavar = "NUMBER", type=int, default = 2)
    parser.add_argument("-m", "--method", help = "Method for clusters correspondance through frames (residue or position)", type=str, default = "residue", choices=['residue', 'position'])
    parser.add_argument("-c", "--nb-corr", help = "Number of largest clusters to take into account for clusters correspondance at time t (default : all)", type = int, metavar = "NUMBER", default = -1)
    parser.add_argument("--frames", help = "Number of frames to process (default : all)", metavar = "NUMBER", type = int, default=-1)

    args = parser.parse_args()

    if not args.prefix :
        args.prefix = ".".join(args.traj.split("/")[-1].split(".")[:-1]) # Get trajectory file name without all directories path and without file extension.
        
    if args.nb_corr != -1 and args.nb_corr < args.to_keep : 
        logging.error("Number of clusters for correspondance calculation (-c/--nb-corr) can't be lower than number of clusters to keep (-n/--to-keep).")
        exit()

    #Create output directory
    if not os.path.isdir(args.outdir):
        logging.info(f"Create {args.outdir} directory")
        os.mkdir(args.outdir)
    
    return args

def plot_size(out_prefix:str, system):
    """Plot size of clusters provided in ClusterIterator through time

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)
    ax.set_title(f"Size of the {system.nb_clusters} largest clusters through time ({ARGS.method} for cluster correspondance)")
    ax.set_ylim(0,300)
    ax.set_xlabel("Time")
    ax.set_ylabel("Cluster size")

    times = [frame.time for frame in system]
    clusters_sizes = []

    for i in range(system.nb_clusters):
        sizes = [frame.clusters[i].size for frame in system]
        clusters_sizes.append(sizes)

    for i, sizes in enumerate(clusters_sizes):
        ax.plot(times, sizes, label = f"Cluster {i}", alpha = 0.7)
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_size_clusters.svg", format = "svg")

def plot_com(out_prefix:str, clusters:molecules_aggregate.clusters.ClusterIterator):
    """Plot absolute center of masses of clusters for all frames along x,y and z axis on "3d" plot. 

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f"Absolute center of mass of clusters ({ARGS.method} for cluster correspondance)")
    #ax.set_xlim(0, system.box_dimension[0])
    #ax.set_ylim(0, system.box_dimension[1])
    #ax.set_zlim(0, system.box_dimension[2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    for i in range(system.nb_clusters):
        ax.scatter([frame.clusters[i].center_of_mass[0] for frame in system], [frame.clusters[i].center_of_mass[1] for frame in system], [frame.clusters[i].center_of_mass[2] for frame in system], label =  f"Cluster {i}", alpha = 0.5, s=10, depthshade = False)
        
    leg = ax.legend(loc="center left")
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([10])
    fig.savefig(f"{out_prefix}_absolute_com.svg", format="svg")  

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

def plot_relative_position(out_prefix, system):
    fig, ax = plt.subplots()
    ax.set_title(f"Relative z mean position")
    ax.set_ylim(0, system.box_dimension[2])
    ax.set_xlabel("Time")
    ax.set_ylabel("z mean position")

    times = [frame.time for frame in system]
    clusters_sizes = []

    #ax.plot(times, [frame.membrane.lowest_z_mean for frame in system], color="grey")


    ax.plot(times, [frame.membrane.highest.z_mean for frame in system], color = "grey", label = "Highest membrane point")
    #ax.plot(times, [frame.membrane.highest.z_mean + frame.membrane.highest.z_std for frame in system], color = "grey", linestyle="--", label = "Standard deviation")
    #ax.plot(times, [frame.membrane.highest.z_mean - frame.membrane.highest.z_std for frame in system], color = "grey", linestyle="--")

    ax.fill_between(times, [frame.membrane.highest.z_mean - frame.membrane.highest.z_std for frame in system], [frame.membrane.highest.z_mean + frame.membrane.highest.z_std for frame in system], color="grey", alpha = 0.5)

    ax.plot(times, [frame.membrane.lowest.z_mean for frame in system], color = "tan", label = "Lowest membrane point")
    #ax.plot(times, [frame.membrane.lowest.z_mean + frame.membrane.lowest.z_std for frame in system], color = "tan", linestyle="--", label = "Standard deviation")
    #ax.plot(times, [frame.membrane.lowest.z_mean - frame.membrane.lowest.z_std for frame in system], color = "tan", linestyle="--")

    ax.fill_between(times, [frame.membrane.lowest.z_mean - frame.membrane.lowest.z_std for frame in system], [frame.membrane.lowest.z_mean + frame.membrane.lowest.z_std for frame in system], color="tan", alpha = 0.5)

    for i in range(system.nb_clusters):
        ax.plot(times, [frame.clusters[i].center_of_mass[2] for frame in system], label = f"Cluster {i}")
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_relative_z_position.svg", format = "svg")

def plot_relative_com(out_prefix, system):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'"Relative" center of mass of clusters ({ARGS.method} for cluster correspondance)')

    #ax.set_title(f"Absolute center of mass of clusters ({ARGS.method} for cluster correspondance)")
    #ax.set_xlim(0, system.box_dimension[0])
    #ax.set_ylim(0, system.box_dimension[1])
    #ax.set_zlim(0, system.box_dimension[2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    for i in range(system.nb_clusters):
        ax.scatter([frame.clusters[i].x_relative_highest for frame in system], [frame.clusters[i].y_relative_highest for frame in system], [frame.clusters[i].z_relative_highest for frame in system], label =  f"Cluster {i}", depthshade = False, alpha = 0.5, s=10)
        
    leg = ax.legend(loc="center left")
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([10])
    fig.savefig(f"{out_prefix}_relative_com.svg", format = "svg")  


if __name__ == "__main__":
    ARGS = args_gestion()

    start = time.time()

    logging.info("Load system...")

    mda_system = mda.Universe(ARGS.topo, ARGS.traj)

    #This lines are for profiling the time
    '''pr = cProfile.Profile()
    pr.enable()'''

    #clusters = molecules_aggregate.load_clusters(system, "TO", ARGS.threshold, ARGS.to_keep, ARGS.frames, ARGS.method, ARGS.nb_corr)

    system = molecules_aggregate.load_system(mda_system, "TO", ARGS.threshold, ARGS.to_keep, "DOPC", 10, ARGS.frames, ARGS.method, ARGS.nb_corr)

    #This lines are for profiling the time
    '''pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())'''

    out_path = ARGS.outdir + "/" + ARGS.prefix

    plot_size(out_path, system)
    plot_com(out_path, system)
    plot_relative_position(out_path, system)
    plot_relative_com(out_path, system)
    #plot_2d(out_path, clusters)
    logging.info(f"Analysis end in {time.time() - start} s")





