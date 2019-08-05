import math
import arrow
import numpy as np
from hypercubeq import HypercubeQ

np.random.seed(0)

# RANDOM INITIALIZED MODEL WITH ZERO-LINE CAPACITY
# - model initialization
start_t = arrow.now()
hq      = HypercubeQ(n_atoms=20, cap="zero", max_iter=10)
end_t   = arrow.now()
print("Calculation time: [%s]" % (end_t - start_t))
# - steady-state probability
print(hq.Pi)
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1)
# - average travel time per dispatch for each unit
print(hq.Tu)



# RANDOM INITIALIZED MODEL WITH ZERO-LINE CAPACITY
# - model initialization
hq = HypercubeQ(n_atoms=20, cap="inf", max_iter=10)
# - steady-state probability
print(hq.Pi)    # steady-state probability for unsaturate states
print(hq.Pi_Q)  # steady-state probability for saturate states (only for infinite-line capacity)
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1) # fraction of all dispatches that send unit n to atom j and incur no queue delay
print(hq.Rho_2) # fraction of all dispatches that send unit n to atom j and do incur a positive
# - average travel time per dispatch for each unit
print(hq.Tu)



# USER CUSTOMIZED MODEL WITH INFINITE-LINE CAPACITY
# - model configuration
n_atoms = 3
Lam     = [1, 1, 1]
P       = [[0, 1, 2],
           [1, 0, 2],
           [2, 0, 1]]
T       = np.random.rand(n_atoms, n_atoms)
# - model initialization
hq = HypercubeQ(n_atoms=3, Lam=Lam, P=P, T=T, cap="zero", max_iter=10)
# - steady-state probability
print(hq.Pi)
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1)
# - average travel time per dispatch for each unit
print(hq.Tu)