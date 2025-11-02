import networkx as nx
import time

from kuramoto import Kuramoto, plot_activity

# Interactions are represented as an adjacency matrix _A_, a 2D numpy ndarray.
# Instantiate a random graph and transform into an adjacency matrix
graph_nx = nx.erdos_renyi_graph(n=10000, p=1)  # p=1 -> all-to-all connectivity
adj_mat = nx.to_numpy_array(graph_nx)

# Instantiate model with parameters
model = Kuramoto(coupling=3, dt=0.01, T=10, n_nodes=len(adj_mat))

# Run simulation - output is time series for all nodes (node vs time)
t0 = time.perf_counter()
activity = model.run(adj_mat=adj_mat)
t1 = time.perf_counter()

print(f"Elapsed: {t1 - t0:.6f} seconds")
