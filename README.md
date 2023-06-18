# julamoto
Explicit Kuramoto model implemented in Julia, covering both the simple Kuramoto model and the Kuramoto model with noise (stochastic or finite-temperature) using different solvers.

Bindings to Python are available, though they are simple and we recommend writing Julia code directly. Functions for obtaining parameters such as instantaneous frequency and phase difference are available, in addition to functions for internal plotting and data output in a NumPy- and Matplotlib-friendly manner. Examples of setup usign Julia and plotting using Python (my preferred method) are provided.

## About the Kuramoto model

The Kuramoto model is a mathematical model used to describe coupled oscillators and the synchronization thereof. It is governed by the following system of differential equations on the following form:
