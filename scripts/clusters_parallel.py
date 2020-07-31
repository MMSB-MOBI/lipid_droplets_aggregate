import MDAnalysis as mda
from pmda.parallel import ParallelAnalysisBase 
import time
from MDAnalysis.lib import distances

def compute_clusters(box_dimension, mol_group, threshold):
    """Compute clusters.
    Center of mass is computed for each TO (or other molecule if other is given) residue. Then neighborhood is computed with MDAnalysis function and this neighborhood is reformated to have direct access to the list of neighbors for each residue, and clustering is done from these lists with clustering() function.  

    Args:
        system (MDAnalysis.universe.Universe): MDAnalysis universe
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


def format_neighborhood(neighbors):
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


def clustering(formated_neighbors, nb_residues):
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


def one_frame_clusters(atom_group, *args):
    threshold = 13
    clusters = compute_clusters(atom_group.universe.dimensions, atom_group, threshold)
    #print("B")
    return clusters

def com(atom_group, *args):
    return atom_group.center_of_mass()


class NewAnalysis(ParallelAnalysisBase):
    def __init__(self, atomgroup, parameter):
        self._ag = atomgroup
        super(NewAnalysis, self).__init__(atomgroup.universe,
                                          self._ag)

    def _single_frame(self, ts, agroups):
        # REQUIRED
        # called for every frame. ``ts`` contains the current time step
        # and ``agroups`` a tuple of atomgroups that are updated to the
        # current frame. Return result of `some_function` for a single
        # frame
        print("single frame")
        return com(agroups[0], self._parameter)

    def _conclude(self):
        # REQUIRED
        # Called once iteration on the trajectory is finished. Results
        # for each frame are stored in ``self._results`` in a per block
        # basis. Here those results should be moved and reshaped into a
        # sensible new variable.
        self.results = np.hstack(self._results)
        # Apply normalisation and averaging to results here if wanted.
        #self.results /= np.sum(self.results

if __name__ == "__main__":
    system = mda.Universe("../data/gro456.PO4.tpr", "../data/md.PO4.xtc")
    cluster_atoms = system.select_atoms(f"resname TO")
    a = []
    print("==serial")
    start = time.time()
    for trj in system.trajectory[:10]:
        if trj.frame % 1000 == 0:
            print(trj.frame)
        a.append(one_frame_clusters(cluster_atoms))
    print("END", time.time() - start)

    na = NewAnalysis(cluster_atoms, 35).run(start = 0, stop = 10)


