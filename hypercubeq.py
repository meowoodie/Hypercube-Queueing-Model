import math
import arrow
import numpy as np
from collections import defaultdict

class HypercubeQ(object):
    """
    Hypercube Queueing Model with Infinite-line Capacity

    * For the sake of simplication, 
    """


    def __init__(self, n_atoms, Lam=None, Mu=None, T=None, P=None, max_iter=10):
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
        # - steady-state probability for unsaturate states
        self.Pi     = self._steady_state_probabilities(max_iter=max_iter)

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
            S[i:i+m,n]  = 1                                 # add "one" at i-th position
            S[i:i+m,:n] = np.flip(S[int(i-m):i,:n], axis=0) # step backwards
            i               += m                            # update index of the state
            m               *= 2                            # update the number of states that needs step backwards
        return S

    def _upward_transition_rates(self):
        """
        An efficient method for generating upward transition rates from state i to state j of the  
        hypercube queueing model.

        For each geographical atom j we shall tour the hypercube in a unit-step fashion. 
        """
        Upoptn = self.__upward_optimal_neighbor() # upward optimal neighbor matrix
        Lam_ij = defaultdict(lambda: 0)           # upward transition rates dictionary initialization

        # iterative algorithm for generating upward transition rates
        for k in range(self.n_atoms):          # for each atom k
            for i in range(2 ** self.n_atoms): # for each state i
                if self.S[i].sum() < self.n_atoms:
                    Lam_ij[(i,Upoptn[k,i])] += self.Lam[k]
        return Lam_ij

    def _steady_state_probabilities(self, cap="zero", max_iter=100):
        """
        A iterative procedure for obtaining the steady state probability on the hypercube. In a manner
        similar to point Jacobi iteration, we use the equation of detailed balance to determine the 
        values at successive iterations. 

        For all idle state, all busy state, and S_Q (more than N customers in the system), the steady   
        state probabilities can be calculated as a normal M/M/N queue with infinite-line capacity.
        """
        # steady state probabilities initialization
        def init_steady_state_prob(cap):
            Pi_0 = np.zeros(2 ** self.n_atoms)
            for n_busy in range(self.n_atoms + 1):
                denominator = sum([
                    self.Lam.sum() ** j / math.factorial(j)
                    for j in range(self.n_atoms + 1) ]) \
                    if cap == "zero" else \
                    sum([ (self.Lam.sum() ** j / math.factorial(j)) + 
                    (self.Lam.sum() ** self.n_atoms / math.factorial(self.n_atoms)) * 
                    (self.Lam.sum() / self.n_atoms / (1 - self.Lam.sum() / self.n_atoms)) 
                    for j in range(self.n_atoms + 1) ])
                n_states    = sum([ 1 for i in range(2 ** self.n_atoms) if self.S[i].sum() == n_busy ])
                init_P      = (self.Lam.sum() ** n_busy / math.factorial(n_busy)) / denominator / n_states
                for i in range(2 ** self.n_atoms):
                    if self.S[i].sum() == n_busy:
                        Pi_0[i] = init_P
            return Pi_0

        # point Jacobi iteration for all states except for all idle state, all busy state
        def iter_steady_state_prob(Pi_n):
            # initialize the steady state probabilities Pi_{n+1} for the next iteration
            Pi_n_1 = np.copy(Pi_n)
            # update the steady state probability of each state i
            for j in range(2 ** self.n_atoms):
                # except for all idle state, all busy state
                if self.S[j].sum() != 0 and self.S[j].sum() != self.n_atoms:
                    upward_sum   = 0
                    downward_sum = 0
                    for i in range(2 ** self.n_atoms):
                        # the k-th response unit that changed status (upward and downward in total)
                        changed_k = np.nonzero(self.S[i] - self.S[j])[0]
                        # only consider one-step changed states
                        if len(changed_k) == 1:
                            # upward transition state
                            if (self.S[j] - self.S[i]).sum() == 1:
                                upward_sum += Pi_n[i] * self.Lam_ij[(i,j)]
                            # downward transition state
                            elif (self.S[j] - self.S[i]).sum() == -1:
                                downward_sum += Pi_n[i] * self.Mu[changed_k[0]]
                    # update the j-th state in Pi_{n+1}
                    Pi_n_1[j] = (upward_sum + downward_sum) / (self.Lam.sum() + (self.Mu * self.S[j]).sum())
            return Pi_n_1
        
        # initialize all states
        Pi = init_steady_state_prob(cap)
        print(Pi)
        print(Pi.sum())
        # update all states except for all idle state, all busy state
        for n in range(max_iter):
            Pi = iter_steady_state_prob(Pi)
        return Pi

    def __upward_optimal_neighbor(self):
        """
        A function that collects the upward optimal neighbor given atom k at state i according 
        to the dispatch policy.
        """
        # helper function that returns the state index given the state
        def state_index(s):
            for i in range(2 ** self.n_atoms):
                if np.count_nonzero(s - self.S[i]) == 0:
                    return i

        # calculate upward optimal neighbors matrix
        Upopts = np.zeros((self.n_atoms, 2 ** self.n_atoms), dtype=int)
        for k in range(self.n_atoms):                # for each atom k
            for i in range(2 ** self.n_atoms):       # for each state i (last state is excluded)
                idle_s = np.where(self.S[i] == 0)[0] # indices of idle response units
                if len(idle_s) != 0:
                    disp_s = self.P[k]               # ordered indices of response units for atom k according to dispatch policy 
                    for s in disp_s:
                        if s in idle_s:
                            add         = np.zeros(self.n_atoms, dtype=int)
                            add[s]      = 1
                            upopts      = self.S[i] + add
                            Upopts[k,i] = state_index(upopts)
                            break
                else:
                    Upopts[k,i] = -1                 # for the all busy state, there is no upward neighbor
        return Upopts


            
if __name__ == "__main__":
    n_atoms = 3
    # Lam     = [1, 1, 1, 1, 1]
    # Mu      = [1, 1, 1, 1, 1]
    # P       = [[0, 1, 2, 3, 4],
    #            [1, 0, 2, 3, 4],
    #            [2, 0, 1, 3, 4],
    #            [3, 0, 1, 2, 4],
    #            [4, 0, 1, 2, 3]]
    Lam     = [1, 1, 1]
    Mu      = [1, 1, 1]
    P       = [[0, 1, 2],
               [1, 0, 2],
               [2, 0, 1]]


    # # * CHECK UPWARD TRANSITION RATES
    # hq = HypercubeQ(n_atoms, Lam=Lam, Mu=Mu, P=P)
    # print(hq.S)
    # state_order = [0, 7, 3, 1, 4, 6, 2, 5]
    # Lam_ij_dict = hq._upward_transition_rates()
    # Lam_ij_mat  = np.zeros((2 ** hq.n_atoms, 2 ** hq.n_atoms))
    # for key in Lam_ij_dict:
    #     i = state_order.index(key[0])
    #     j = state_order.index(key[1])
    #     Lam_ij_mat[i,j] = Lam_ij_dict[key]
    # print(Lam_ij_mat)

    # # * CHECK UPWARD OPTIMAL NEIGHBORS
    # hq  = HypercubeQ(n_atoms, Lam=Lam, Mu=Mu, P=P)
    # print(hq.S)
    # Upn = hq._upward_optimal_neighbor()
    # k   = 2
    # i   = 5
    # print(hq.S[i])
    # print(k)
    # print("neighbors")
    # print(hq.S[Upn[k, i]])

    # * CHECK STEADY STATE PROBABILITY
    hq = HypercubeQ(n_atoms=5, Mu=np.ones(5), max_iter=10)
    print(hq.Pi)
    print(hq.Pi.sum())