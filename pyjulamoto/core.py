from julia.api import Julia
jl = Julia(compiled_modules=False)

# Unnecessary is requiring global install
#from julia import Pkg
#Pkg.activate("Kuramoto")  # relative to project root
#Pkg.instantiate()

from julia import Kuramoto as jK

class Kuramoto:
    def __init__(self, u0, omega, K, tstart, tend, dt):
        self._model = jK.KuramotoModel(u0, omega, K, tstart, tend, dt)
    def run_static(self, abstol=1e-6, reltol=1e-3):
        return jK.run_static(self._model, abstol=abstol, reltol=reltol)
    @property
    def solution(self):
        return self._model.sol
