import arrow
import numpy as np
from collection import defaultdict

class HypercubeQ(object):
    """
    Hypercube Queueing Model

    * For the sake of simplication, 
    """


    def __init__(self, n_atoms, Lam=None, Mu=None, T=None, P=None):
        # Model configuration
        # - number of geographical atoms (response units)
        self.n_atoms = n_atoms                    
        # - arrival rates vector: arrival rates for each atom
        self.Lam     = np.array(Lam, dtype=float) if Lam is not None else np.random.rand(self.n_atoms)
        # - service rates vector: service rates for each response unit
        self.Mu      = np.array(Mu, dtype=float) if Mu is not None else np.random.rand(self.n_atoms)
        # - traffic matrix:       average traffic time from atom i to atom j
        self.T       = np.array(T, dtype=float) if T is not None else np.random.rand(self.n_atoms, self.n_atoms)
        # - preference matrix:    preference matrix indicates the priority of response units to each atom
        self.P       = np.array(P, dtype=int) if P is not None else np.random.rand(self.n_atoms, self.n_atoms).argsort()

        assert self.n_atoms == self.Lam.shape[0] == self.Mu.shape[0] == \
               self.T.shape[0] == self.T.shape[1] == self.P.shape[0] == self.P.shape[1], \
               "Invalid shape of input parameters."

        # Model status
        # - state space ordered as a sequence represents a complete unit-step tour of the hypercube
        self.S      = self._tour()
        # - upward transition rates matrix: a dictionary { (i,j) : lam_ij } (due to the sparsity of the matrix in nature)
        self.Lam_ij = self._upward_transition_rates()
        # - steady-state probability of states
        self.Pi     = np.zeros(2 ** self.n_atoms)

    def _tour(self):
        """
        Tour Algorithm

        In order to tour the hypercube in a unit-step manner, this function is able to generate a 
        complete sequence S_1, S_2, ... of N-digit binary numbers, with 2^N unique members in the 
        sequence and with adjacent members being exactly unit-Hamming distance apart. Such a 
        sequence represents a complete unit-step tour of the hypercube.
        """
        # initialization
        S      = np.zeros((2 ** self.n_atoms, self.n_atoms))
        S[1,0] = 1    # S_0 = 0, S_1 = 1
        m, i   = 2, 2 # m: number of states that needs step backwards; i: index of the state
        
        # add "one" at i-th position and step backwards
        for n in range(1, self.n_atoms):
            S[i:i+m,n]  = 1                                      # add "one" at i-th position
            S[i:i+m,:n] = np.flip(S[int(i-m):i,:n], axis=0) # step backwards
            i               += m                                      # update index of the state
            m               *= 2                                      # update the number of states that needs step backwards
        return S

    def _upward_transition_rates(self):
        """
        An efficient method for generating upward transition rates from state i to state j of the  
        hypercube queueing model.

        For each geographical atom j we shall tour the hypercube in a unit-step fashion. 
        """
        Upn    = self.__upward_neighbor()
        Lam_ij = defaultdict(lambda: 0)
        # iterative algorithm for generating upward transition rates
        for k in range(self.n_atoms):          # for each atom k
            for i in range(2 ** self.n_atoms): # for each state i
                for j in Upn[i]:               # for each adjacent state j that d_ij^+ = 1
                    Lam_ij[(i,j)] += self.Lam[k] 

        return Lam_ij

    def __upward_neighbor(self):
        """
        Helper function that collects the upward neighbors for each state of the hypercube, and 
        organizes them into a matrix where the key represents each state, and the value includes
        its upward neighbors.
        """
        Upn = defaultdict(lambda: [])
        for i in range(2 ** self.n_atoms):
            for j in range(2 ** self.n_atoms):
                if (self.S[j] - self.S[i]).sum() == 1:
                    Upn[i].append(j)
        return Upn

                    

    # def _steady_state_probabilities(self):
    #     """
    #     """
    #     pass

            


if __name__ == "__main__":
    n_atoms = 5
    Lam     = [1, 1, 2]
    Mu      = [1, 1, 1]
    T       = [[.5, .1, .1], 
               [.2, .2, .2],
               [.1, .5, .5]]
    P       = [[0, 1, 2],
               [1, 0, 2],
               [2, 0, 1]]

    hq = HypercubeQ(n_atoms)
    # hq._tour()
    print(hq.S)