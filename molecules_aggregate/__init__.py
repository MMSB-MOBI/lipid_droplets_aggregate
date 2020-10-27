from . import trajectory
import logging
from . import config


def compute_membranes_and_clusters():
    logging.info("== Compute clusters and membranes")
    frames = trajectory.TrajectoryIterator()
    frames.compute_membranes()
    frames.compute_clusters()
    return frames

