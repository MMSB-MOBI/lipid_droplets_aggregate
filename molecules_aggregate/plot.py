import matplotlib.pyplot as plt
from . import config
import logging

global COORD_IDX 
COORD_IDX = {"x": 0, "y": 1, "z": 2}

def clusters_size(out_prefix:str, system):
    """Plot size of clusters provided in ClusterIterator through time

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)
    ax.set_title(f"Size of the {config.NB_CLUSTERS} largest clusters through time ({config.CORRESPONDANCE} for cluster correspondance)")
    ax.set_ylim(0,300)
    ax.set_xlabel("Time")
    ax.set_ylabel("Cluster size")

    times = [frame.time for frame in system]
    clusters_sizes = []

    for i in range(config.NB_CLUSTERS):
        sizes = [frame.clusters[i].size for frame in system]
        clusters_sizes.append(sizes)

    for i, sizes in enumerate(clusters_sizes):
        ax.plot(times, sizes, label = f"Cluster {i}", alpha = 0.7)
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_size_clusters.svg", format = "svg")
    logging.info(f"Clusters size saved to {out_prefix}_size_clusters.svg")

def absolute_center_of_mass(out_prefix:str, system):
    """Plot absolute center of masses of clusters for all frames along x,y and z axis on "3d" plot. 

    Args:
        out_prefix (str): prefix for output png file
        clusters (molecules_aggregate.clusters.ClusterIterator): ClusterIterator object given by molecule_aggregate.load_clusters()
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f"Absolute center of mass of clusters ({config.CORRESPONDANCE} for cluster correspondance)")
    ax.set_xlim(0, config.UNIVERSE.dimensions[0])
    ax.set_ylim(0, config.UNIVERSE.dimensions[1])
    ax.set_zlim(0, config.UNIVERSE.dimensions[2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    for i in range(config.NB_CLUSTERS):
        ax.scatter([frame.clusters[i].center_of_mass[0] for frame in system], [frame.clusters[i].center_of_mass[1] for frame in system], [frame.clusters[i].center_of_mass[2] for frame in system], label =  f"Cluster {i}", alpha = 0.5, s=10, depthshade = False)
        
    leg = ax.legend(loc="center left")
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([10])
    fig.savefig(f"{out_prefix}_absolute_com.svg", format="svg")  
    logging.info(f"Clusters absolute center of mass saved to {out_prefix}_absolute_com.svg")

def absolute_coord_through_time(out_prefix, system, coord):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)
    ax.set_title(f"{coord} coordinate of the {config.NB_CLUSTERS} largest clusters through time ({config.CORRESPONDANCE} for cluster correspondance)")
    ax.set_ylim(0,config.UNIVERSE.dimensions[COORD_IDX[coord]])
    ax.set_xlabel("Time")
    ax.set_ylabel(f"Cluster {coord} coordinate")

    times = [frame.time for frame in system]

    for i in range(config.NB_CLUSTERS):
        ax.plot(times, [frame.clusters[i].center_of_mass[COORD_IDX[coord]] for frame in system], label = f"Cluster {i}", alpha = 0.7)
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_absolute_{coord}.svg", format = "svg")
    logging.info(f"Clusters absolute {coord} saved to {out_prefix}_absolute_{coord}.svg")

def relative_coord_through_time(out_prefix, system, coord):
    fig, ax = plt.subplots()
    ax.set_title(f"Relative {coord} mean position")
    ax.set_ylim(0, config.UNIVERSE.dimensions[COORD_IDX[coord]])
    ax.set_xlabel("Time")
    ax.set_ylabel(f"{coord} mean position")

    times = [frame.time for frame in system]

    ax.plot(times, [frame.membrane.highest_point[coord]["mean"] for frame in system], color = "grey", label = "Highest membrane point")
    ax.fill_between(times, [frame.membrane.highest_point[coord]["mean"] - frame.membrane.highest_point[coord]["std"] for frame in system], [frame.membrane.highest_point[coord]["mean"] + frame.membrane.highest_point[coord]["std"] for frame in system], color="grey", alpha = 0.5)

    ax.plot(times, [frame.membrane.lowest_point[coord]["mean"] for frame in system], color = "tan", label = "Lowest membrane point")
    ax.fill_between(times, [frame.membrane.lowest_point[coord]["mean"] - frame.membrane.lowest_point[coord]["std"] for frame in system], [frame.membrane.lowest_point[coord]["mean"] + frame.membrane.lowest_point[coord]["std"] for frame in system], color="tan", alpha = 0.5)

    for i in range(config.NB_CLUSTERS):
        ax.plot(times, [frame.clusters[i].center_of_mass[COORD_IDX[coord]] for frame in system], label = f"Cluster {i}")
    
    ax.legend()
    
    fig.savefig(f"{out_prefix}_relative_{coord}.svg", format = "svg")
    logging.info(f"Clusters relative {coord} saved to {out_prefix}_relative_{coord}.svg")

    