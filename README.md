# julamoto
Explicit Kuramoto model implemented in Julia, covering both the simple Kuramoto model and the Kuramoto model with noise (stochastic or finite-temperature) using different solvers.

Bindings to Python will be available, though they will be simple and we recommend writing Julia code directly. Functions for obtaining parameters such as instantaneous frequency and phase difference are available, in addition to functions for internal plotting and data output in a NumPy- and Matplotlib-friendly manner. Examples of setup usign Julia and plotting using Python (my preferred method) are provided.

## About the Kuramoto model

The Kuramoto model is a mathematical model used to describe coupled oscillators and the synchronization thereof. It is governed by the following system of differential equations on the following form:

$$\dot{\theta}_i = F_i + \frac{K}{N}\sum_j^N\text{sin}(\theta_j - \theta_i)$$

where $i$ and $j$ are oscillator indices, $K$ is a coupling term, $N$ is the number of oscillators and $F_i$ is $i$th oscillator natural frequency. In the absence of any interaction ($K=0$), the natural frequency is the actual frequency of the oscillators at all times. The oscillators will deviate from this natural frequency due to the coupling with all other oscillators, governed by a coupling term $K$. This may be a constant, resulting in uniform all-to-all coupling, but we allow it to be an NxN adjacency matrix such that the oscillators can be coupled in a variety of manners, e.g. uniform, non-uniform, local, unidirectional, self-coupled etc. In the typical all-to-all setup with no self-coupling is a symmetrix matrix with all diagonal elements set to 0. This coupling gives rise to strong synchronization for a certain range of natural frequencies, where the instantaneous frequencies become the same. However, an assumption of the Kuramoto model is weak coupling, and so any synchronization that happens is abrupt and only for a limited range of natural frequencies.
