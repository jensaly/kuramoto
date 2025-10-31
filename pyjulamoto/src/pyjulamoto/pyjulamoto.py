from julia.api import Julia
jl = Julia(compiled_modules=False)

from julia import Main

Main.include("src/Kuramoto.jl")
jlmod = Main.Kuramoto

class Kuramoto:
    def __init__(self, u0, omega, K, tstart, tend, dt, D=0.0):
        self._model = jlmod.KuramotoModel(u0, omega, K, tstart, tend, dt, D=D)

    def run_static(self, abstol=1e-6, reltol=1e-3):
        return jlmod.run_static!(self._model, abstol, reltol)

    def run_dynamic(self, abstol=1e-6, reltol=1e-3):
        return jlmod.run_dynamic!(self._model, abstol, reltol)

    def run_stochastic(self, abstol=1e-6, reltol=1e-3):
        return jlmod.run_static_stochastic!(self._model, abstol, reltol)

    @property
    def solution(self):
        return self._model.sol
