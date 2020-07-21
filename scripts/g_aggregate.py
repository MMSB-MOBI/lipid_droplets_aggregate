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
    parser.add_argument("-f", "--traj", help = "Trajectory (all formats accepted by MDAnalysis)", required = True,  type=str)
    parser.add_argument("-s", "--topo", help = "Topology (all formats accepted by MDAnalysis)", required = True, metavar ="FILE", type=str)
    parser.add_argument("-o", "--outdir", help = "Output directory (default : .)", default=".", metavar = "DIR", type=str)
    parser.add_argument("-p", "--prefix", help = "Output prefix (default : Trajectory file name)", metavar = "STR", type=str)
    parser.add_argument("-t", "--threshold", help = "Threshold for clustering (default : 13)", metavar = "NUMBER", type=float, default=13)
    parser.add_argument("-n", "--to-keep", help = "Number of largest clusters to keep for the results report (default : 2)", metavar = "NUMBER", type=int, default = 2)
    parser.add_argument("--frames", help = "Number of frames to process (default : all)", metavar = "NUMBER", type = int, default=0)
    parser.add_argument("-m", "--method", help = "Method for clusters correspondance through frames (residue or position)", type=str, default = "residue", choices=['residue', 'position'])
    args = parser.parse_args()

    if not args.prefix :
        args.prefix = ".".join(args.traj.split("/")[-1].split(".")[:-1]) # Get trajectory file name without all directories path and without file extension.
    
    return args

def plot_size(out_prefix, clusters):
    """Plot clusters size results for the first X clusters (X is given by user), with one line per cluster. 

    Args:
        out_prefix ([type]): The prefix for output image name. 
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

def plot_com(out_prefix, clusters):
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

    ax.legend(loc="center left")
    fig.savefig(f"{out_prefix}_absolute_com.png")  

def plot_2d(out_prefix, clusters):
    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through x/z axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[0])
    ax.set_ylim(0,system.dimensions[2])
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[0] for c in cluster_frames], [c.center_of_mass[2] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}")
    
    ax.legend()
    fig.savefig(f"{out_prefix}_xz.png")

    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through y/x axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[1])
    ax.set_ylim(0,system.dimensions[0])
    ax.set_xlabel("y")
    ax.set_ylabel("x")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[1] for c in cluster_frames], [c.center_of_mass[0] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}")
    ax.legend()
    fig.savefig(f"{out_prefix}_yx.png")

    fig, ax = plt.subplots()
    ax.set_title(f"Absolute center of mass of clusters. Projection through y/z axis. ({ARGS.method} for cluster correspondance)")
    ax.set_xlim(0,system.dimensions[1])
    ax.set_ylim(0,system.dimensions[2])
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    for cluster_frames in clusters:
        ax.scatter([c.center_of_mass[1] for c in cluster_frames], [c.center_of_mass[2] for c in cluster_frames], s = 10, label =  f"Cluster {cluster_frames[0].idx}")
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

    '''pr = cProfile.Profile()
    pr.enable()'''

    clusters = molecules_aggregate.load_clusters(system, "TO", ARGS.threshold, ARGS.to_keep, ARGS.frames, ARGS.method)

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





