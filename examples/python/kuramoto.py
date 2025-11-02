import numpy as np
from pyjulamoto import Kuramoto

# Parameters
N = 10  # number of oscillators
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
model.run_static(abstol=abstol, reltol=reltol)

# Access solution
sol = model.solution

# Print final phases
print("Final oscillator phases:", sol.u[end])
