import arrow
import numpy as np

class HypercubeQ(object):
    """
    Hypercube Queueing Model

    For the sake of simplication, 
    """


    def __init__(self, n_atoms, Lam, Mu, T, P):
        self.n_atoms = n_atoms                    # number of geographical atoms (response units)
        self.Lam     = np.array(Lam, dytpe=float) # arrival rates vector: arrival rates for each atom
        self.Mu      = np.array(Mu, dytpe=float)  # service rates vector: service rates for each response unit
        self.T       = np.array(T, dytpe=float)   # traffic matrix:       average traffic time from atom i to atom j
        self.P       = np.array(P, dtype=int)     # preference matrix:    preference matrix indicates the priority of response units to each atom

        assert self.n_atoms == self.Lam.shape[0] == self.Mu.shape[0] == \
               self.T.shape[0] == self.T.shape[1] == self.P.shape[0] == self.P.shape[1], \
               "Invalid shape of input parameters."
        
        self.S = np.zeros((2 ** n_atoms, n_atoms))

    def _tour(self):
        """
        In order to tour the hypercube in a unit-step manner.
        """
        # initialization
        m, i        = 2, 2 # m: number of states that needs step backwards; i: index of the state
        self.S[1,0] = 1    # S_0 = 0, S_1 = 1
        # add "one" at i-th position and step backwards
        for n in range(1, self.n_atoms):
            self.S[i:i+m,n]  = 1                            # add "one" at i-th position
            self.S[i:i+m,:n] = reverse(self.S[i-m/2,i, :n]) 
            




if __name__ == "__main__":
    n_atoms = 3
    Lam     = [1, 1, 2]
    Mu      = [1, 1, 1]
    T       = [[.5, .1, .1], 
               [.2, .2, .2],
               [.1, .5, .5]]
    P       = [[0, 1, 2],
               [1, 0, 2],
               [2, 0, 1]]

    HypercubeQ(n_atoms, Lam, Mu, T, P)