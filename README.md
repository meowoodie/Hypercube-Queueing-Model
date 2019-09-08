Hypercube Queueing Model
===

This repository is an Python implementation of hypercube queueing model proposed by Richard C. Larson (1974). It was originally implemented for analyzing the operational dynamics of Atlanta patrol system in our research [Atlanta zone redesign](https://github.com/meowoodie/Zoning-Analysis). See more details in the implementation `hypercubeq.py`.

![hypercube](https://github.com/meowoodie/Hypercube-Queueing-Model/blob/master/img/hypercube.png)

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
An example of output:
```bash
[2019-09-08T00:03:55.658167-04:00] calculating upward transition rates ...
[2019-09-08T00:03:55.665208-04:00] calculating steady-state probabilities ...
[2019-09-08T00:03:55.726440-04:00] calculating dispatch fraction ...
[2019-09-08T00:03:55.733735-04:00] calculating average travel time ...
Calculation time: [0:00:00.086639]
[0.36809816 0.01857178 0.00605335 0.08913722 0.03191012 0.00465053
 0.00506265 0.07306268 0.02369345 0.00385922 0.00326617 0.0121829
 0.02937934 0.00558559 0.01579855 0.1156378  0.02918155 0.00661827
 0.00327475 0.00895581 0.00422503 0.00306748 0.00225018 0.00638114
 0.01473322 0.00308174 0.00232129 0.00682475 0.0181038  0.00320974
 0.01013305 0.07168868]
1.0
[[0.03108382 0.01395226 0.04896924 0.00081446 0.00197556]
 [0.00043689 0.19428229 0.00608311 0.02546951 0.00665939]
 [0.00451146 0.04141935 0.00069833 0.15410511 0.00049082]
 [0.00104905 0.00058733 0.15495385 0.00881085 0.10890997]
 [0.15648301 0.00200309 0.0018876  0.00297811 0.03138554]]
0.9999999999999999
[0.73757693 0.5606971  0.60706605 0.33816496 0.59965686]
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