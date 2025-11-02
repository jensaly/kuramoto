import numpy as np
import time
from pyjulamoto import Kuramoto

# Parameters
N = 10000  # number of oscillators
tstart = 0.0
tend = 10.0
dt = 0.01
abstol = 1e-6
reltol = 1e-3

# Initial conditions
u0 = np.random.rand(N) * 2 * np.pi  # random initial phases
omega = np.random.randn(N)          # natural frequencies

# Coupling matrix (all-to-all with zero diagonal)
K = np.ones((N, N)) / N
np.fill_diagonal(K, 0.0)

# Create Kuramoto model
model = Kuramoto(u0, omega, K, tstart, tend, dt)

# Run static solver
t0 = time.perf_counter()
model.run_static(abstol=abstol, reltol=reltol)
t1 = time.perf_counter()

# Access solution
sol = model.solution

print(f"Elapsed: {t1 - t0:.6f} seconds")