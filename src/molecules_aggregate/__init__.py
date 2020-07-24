from . import clusters
from . import membrane

def load_clusters(mda_system, molecules_to_cluster:str, clustering_threshold:float, clusters_to_keep:int, frames_to_process:int,cluster_correspondance:str, nb_cluster_for_corr:int):
    return clusters.ClusterIterator(mda_system, molecules_to_cluster, clustering_threshold, clusters_to_keep, frames_to_process, cluster_correspondance, nb_cluster_for_corr)

#def load_membrane(mda_system, membrane_molecule, slice_size):
#    return membrane.Membrane(mda_system, membrane_molecule, slice_size)