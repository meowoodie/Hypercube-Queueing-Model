import arrow
import numpy as np

class HypercubeQ(object):
    """
    Hypercube Queueing Model

    For the sake of simplication, 
    """


    def __init__(self, n_atoms):
        Lam = np.zeros((n_atoms))          # arrival rates vector: arrival rates for each atom
        Mu  = np.zeros((n_atoms))          # service rates vector: service rates for each response unit
        T   = np.zeros((n_atoms, n_atoms)) # traffic matrix:       average traffic time from atom i to atom j
        P   = np.zeros((n_atoms, n_atoms)) # preference matrix:    preference matrix indicates the priority of response units to each atom

        
        