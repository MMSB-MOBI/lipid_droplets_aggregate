def init(universe, args):
    """Define global variable for library based on arguments provided by user. 

    Args:
        universe (MDAnalysis.Universe): MDAnalysis universe to work witb
        args (argparse.Namespace): argparse parsed arguments
    """
    global UNIVERSE
    global NB_FRAMES
    global NB_PROCESS
    global MEMBRANE_MOL
    global SLICE_AXIS
    global SLICE_AXIS_IDX
    global SLICE_SIZE
    global MEMBRANE_AXIS
    global CLUSTER_MOL
    global CLUSTER_THRESHOLD
    global NB_CORRESPONDANCE
    global CORRESPONDANCE
    global NB_CLUSTERS
    
    UNIVERSE = universe
    NB_FRAMES = args.frames
    NB_PROCESS = args.process
    MEMBRANE_MOL = args.m_mol
    SLICE_AXIS = args.slice_axis
    SLICE_AXIS_IDX = args.slice_axis_idx
    SLICE_SIZE = args.slice_size
    MEMBRANE_AXIS = args.membrane_axis
    CLUSTER_MOL = args.c_mol
    CLUSTER_THRESHOLD = args.threshold
    NB_CORRESPONDANCE = args.nb_corr
    CORRESPONDANCE = args.method
    NB_CLUSTERS = args.to_keep