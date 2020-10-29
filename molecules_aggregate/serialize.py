import logging

def write_results_size(out_prefix:str, system):
    filename = out_prefix + "_clusters_size.tsv"
    headers = [f"Cluster {i} size" for i in range (system.nb_clusters)]
    with open(filename, "w") as o:
        o.write("#Time\t" + "\t".join(headers) + "\n")
        for frame in system:
            clusters_size_str = '\t'.join([str(cluster.size) for cluster in frame.clusters])
            o.write(f"{frame.time}\t{clusters_size_str}\n")

    logging.info(f"Size raw results saved to {filename}")
            

def write_results_position(out_prefix:str, system):
    filename = out_prefix + "_clusters_position.tsv"
    headers = [f"Cluster {i} {axis}" for i in range (system.nb_clusters) for axis in ['x', 'y', 'z']]
    with open(filename, "w") as o:
        o.write("#Time\t" + "\t".join(headers) + "\n")
        for frame in system: 
            position_str = "\t".join([str(coord) for cluster in frame.clusters for coord in cluster.center_of_mass])
            o.write(f"{frame.time}\t{position_str}\n")

    logging.info(f"Position raw results saved to {filename}")
        
def write_results_with_membrane(out_prefix:str, system, axis):
    filename = out_prefix + "_clusters_membrane_position.tsv"
    headers = [f"Cluster {i} {axis}" for i in range (system.nb_clusters) for axis in ['x', 'y', 'z']]
    with open(filename, "w") as o:
        o.write("#Time\t" + "\t".join(headers) + f"\tMembrane highest {axis} mean\tMembrane highest {axis} std\tMembrane lowest {axis} mean\tMembrane lowest {axis} std" + "\n")
        for frame in system:
            position_str = "\t".join([str(coord) for cluster in frame.clusters for coord in cluster.center_of_mass])
            o.write(f'{frame.time}\t{position_str}\t{frame.membrane.highest_point[axis]["mean"]}\t{frame.membrane.highest_point[axis]["std"]}\t{frame.membrane.lowest_point[axis]["mean"]}\t{frame.membrane.lowest_point[axis]["std"]}\n')
    
    logging.info(f"Position raw results with membrane saved to {filename}")