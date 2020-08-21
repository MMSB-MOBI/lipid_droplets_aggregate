from . import trajectory
import logging
from . import config

def compute_aggregation():
    logging.info("== Compute aggregation on universe ==")
    return trajectory.TrajectoryIterator()
