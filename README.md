Hypercube Queueing Model
===

This repository is an Python implementation of hypercube queueing model proposed by Richard C. Larson (1974). It was originally implemented for analyzing the operational dynamics of Atlanta patrol system in our research [Atlanta zone redesign](https://github.com/meowoodie/Zoning-Analysis). See more details in the implementation `hypercubeq.py`.

# Demo

There are three demos provided in `demo.py`:

1. Hypercube queueing model with random initialized parameters (including the arrival rates of each graphical atom, the traffic matrix, and the dispatch preference) and zero-line capacity:
```Python
# - model initialization
start_t = arrow.now()
hq      = HypercubeQ(n_atoms=5, cap="zero", max_iter=10)
end_t   = arrow.now()
print("Calculation time: [%s]" % (end_t - start_t))
# - steady-state probability
print(hq.Pi)
print(hq.Pi.sum())
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1)
print(hq.Rho_1.sum())
# - average travel time per dispatch for each unit
print(hq.Tu)
```

2. Hypercube queueing model with random initialized parameters (including the arrival rates of each graphical atom, the traffic matrix, and the dispatch preference) and infinite-line capacity:
```Python 
# - model initialization
hq = HypercubeQ(n_atoms=5, cap="inf", max_iter=10)
# - steady-state probability
print(hq.Pi)    # steady-state probability for unsaturate states
print(hq.Pi_Q)  # steady-state probability for saturate states (only for infinite-line capacity)
print(hq.Pi.sum() + hq.Pi_Q.sum())
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1) # fraction of all dispatches that send unit n to atom j and incur no queue delay
print(hq.Rho_2) # fraction of all dispatches that send unit n to atom j and do incur a positive
print(hq.Rho_1.sum() + hq.Rho_2.sum())
# - average travel time per dispatch for each unit
print(hq.Tu)
```

3. Hypercube queueing model given all parameters (including the arrival rates of each graphical atom, the traffic matrix, and the dispatch preference) and zero-line capacity:
```Python
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
print(hq.Pi.sum())
# - fraction of dispatches that send a unit n to a particular geographical atom j
print(hq.Rho_1)
print(hq.Rho_1.sum())
# - average travel time per dispatch for each unit
print(hq.Tu)
```

# References
- ["A hypercube queuing model for facility location and redistricting in urban emergency services". Richard C. Larson.](https://www.sciencedirect.com/science/article/pii/0305054874900768)