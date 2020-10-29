#!/usr/bin/env python3

import argparse
import MDAnalysis as mda
import molecules_aggregate
from molecules_aggregate import error, plot, serialize

import time
import logging
logging.basicConfig(level = logging.INFO, format='%(levelname)s\t%(message)s')
import cProfile, pstats, io
from pstats import SortKey
from mpl_toolkits.mplot3d import Axes3D
import os
import multiprocessing, psutil
from molecules_aggregate.utils import check_ram
import molecules_aggregate.plot as plot

class CheckFileExistence(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(CheckFileExistence, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string = None):
        if not os.path.exists(values):
            raise error.ArgumentError(f"{values} file doesn't exist")
        setattr(namespace, self.dest, values)

class CreateDirectory(argparse.Action):
    """Create directory given as argument if not exists"""
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(CreateDirectory, self).__init__(option_strings, dest, **kwargs)
    
    def __call__(self, parser, namespace, values, option_string = None):
        if not os.path.exists(values):
            os.makedirs(values)
        else:
            logging.warning(f"Output directory {values} already exists.")
        setattr(namespace, self.dest, values)


def args_gestion():
    parser = argparse.ArgumentParser(description="Script to identify TO clusters size and position along the MD trajectory.")
    group_io = parser.add_argument_group("General options")
    group_membrane = parser.add_argument_group("Membrane options")
    group_cluster = parser.add_argument_group("Clustering options")
    
    group_io.add_argument("-f", "--traj", help = "Trajectory (all formats accepted by MDAnalysis)", required = True,  type=str, metavar = "FILE", action = CheckFileExistence)
    group_io.add_argument("-s", "--topo", help = "Topology (all formats accepted by MDAnalysis)", required = True, metavar ="FILE", type=str, action = CheckFileExistence)
    group_io.add_argument("-o", "--outdir", help = "Output directory (default : .)", default=".", metavar = "DIR", type=str, action = CreateDirectory)
    group_io.add_argument("-p", "--prefix", help = "Output prefix (default : Trajectory file name)", metavar = "STR", type=str)
    group_cluster.add_argument("--c-mol", help = "Clusters molecule (default : TO)", metavar = "STR", type = str, default = "TO")
    group_cluster.add_argument("-t", "--threshold", help = "Threshold for clustering (default : 13)", metavar = "NUMBER", type=float, default=13)
    group_cluster.add_argument("-n", "--to-keep", help = "Number of largest clusters to keep (default : 2)", metavar = "NUMBER", type=int, default = 2)
    group_cluster.add_argument("-m", "--method", help = "Method for clusters correspondance through frames (residue or position) (default : residue)", type=str, default = "residue", choices=['residue', 'position'])
    group_cluster.add_argument("-c", "--nb-corr", help = "Number of largest clusters to take into account for clusters correspondance at time t (default : all)", type = int, metavar = "NUMBER", default = -1)
    group_io.add_argument("--frames", help = "Number of frames to process (default : all)", metavar = "NUMBER", type = int, default=-1)
    group_membrane.add_argument("--m-mol", help = "Membrane molecule (default : DOPC)", metavar = "STR", type = str, default = "DOPC")
    group_membrane.add_argument("--slice-size", help = "Membrane slice size (default : 10)", metavar = "FLOAT", type = float, default = 10)
    group_membrane.add_argument("--slice-axis", help = "Axis to slice membrane (default : x)", metavar = "STR", type = str, default = "x", choices = ["x", "y", "z"])
    group_io.add_argument("--process", help = "Number of process to run for membrane computation (default : available cpus)", metavar = "INT", type = int, default = multiprocessing.cpu_count())
    group_membrane.add_argument("--membrane-axis", help = "Axis to determine membrane highest and lowest points (default : z)", default = "z", type=str, choices=["x", "y", "z"])
    args = parser.parse_args()

    #Set axis
    axis_trans = {"x" : 0, "y" : 1, "z" : 2}
    args.slice_axis_idx = axis_trans[args.slice_axis]
    #args.membrane_axis_idx = axis_trans[args.membrane_axis]
    
    return args

def check_arguments():
    if ARGS.frames > len(UNIVERSE.trajectory):
        raise error.ArgumentError(f"Universe only contains {len(UNIVERSE.trajectory)} frames. You can't compute {ARGS.frames}")
    
    if ARGS.nb_corr != -1 and ARGS.nb_corr < ARGS.to_keep : 
        raise error.ArgumentError("Number of clusters for correspondance calculation (-c/--nb-corr) can't be lower than number of clusters to keep (-n/--to-keep).")
    
    if ARGS.slice_size > UNIVERSE.dimensions[ARGS.slice_axis_idx]:
        raise error.ArgumentError(f"Slice size can't be larger than box dimension ({ARGS.slice_size} > {UNIVERSE.dimensions[ARGS.slice_axis]})")

    ARGS.frames = len(UNIVERSE.trajectory) if ARGS.frames == -1 else ARGS.frames
    ARGS.prefix = ARGS.prefix if ARGS.prefix else ".".join(ARGS.traj.split("/")[-1].split(".")[:-1])

if __name__ == "__main__":
    ARGS = args_gestion()

    UNIVERSE = mda.Universe(ARGS.topo, ARGS.traj)

    check_arguments()

    logging.info(f'== Molecules aggregate configuration ==\n\
        trajectory : {ARGS.traj}\n\
        topology : {ARGS.topo}\n\
        outputs : {ARGS.outdir}/{ARGS.prefix}_*.svg\n\
        number of frames to compute : {ARGS.frames}\n\
        == Clustering\n\
        clusters molecule : {ARGS.c_mol}\n\
        clustering threshold : {ARGS.threshold}\n\
        number of clusters to keep : {ARGS.to_keep}\n\
        method for clusters correspondance : {ARGS.method}\n\
        number of cluster to look at for correspondance : {"all" if ARGS.nb_corr == -1 else ARGS.nb_corr}\n\
        == Membrane\n\
        membrane molecule : {ARGS.m_mol}\n\
        membrane slice size : {ARGS.slice_size}\n\
        membrane slice axis : {ARGS.slice_axis} ({ARGS.slice_axis_idx})\n\
        membrane highest and lowest point axis : {ARGS.membrane_axis}\n\
        process for computation : {ARGS.process}')


    start = time.time()
    #check_ram.start() #To check RAM usage and abort if it's too high 

    molecules_aggregate.config.init(UNIVERSE, ARGS) #Create global variable with user arguments

    computed_trajectory = molecules_aggregate.compute_membranes_and_clusters() #Compute membranes and clusters through frames

    plot_prefix = ARGS.outdir + "/" + ARGS.prefix

    serialize.write_results_size(plot_prefix, computed_trajectory)
    serialize.write_results_with_membrane(plot_prefix, computed_trajectory, ARGS.membrane_axis)
    plot.clusters_size(plot_prefix, computed_trajectory) #Plot clusters size through time
    plot.absolute_center_of_mass(plot_prefix, computed_trajectory) #Plot absolute center of mass position of each cluster, in 3d plot
    plot.absolute_coord_through_time(plot_prefix, computed_trajectory, "x") #Plot absolute x through time
    plot.absolute_coord_through_time(plot_prefix, computed_trajectory, "y") #Plot absolute y through time
    plot.absolute_coord_through_time(plot_prefix, computed_trajectory, "z") #Plot absolute z through time
    plot.relative_coord_through_time(plot_prefix, computed_trajectory, ARGS.membrane_axis) #Plot clusters z mean through time, and x membrane highest and lowest point mean 

    logging.info(f"Analysis end in {time.time() - start} s")
    #check_ram.stop()





