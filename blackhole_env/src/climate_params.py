from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Callable
from constants import SIGMA_SB



@dataclass
class ClimateParams:
    albedo: float                  # A
    emissivity: float              # epsilon
    column_density_mol: float      # chi [mol m^-2]
    molar_heat_capacity: float     # cp [J mol^-1 K^-1]

    @property
    def heat_capacity_area(self) -> float:
        """Capacidade térmica efetiva por área [J m^-2 K^-1]."""
        return self.column_density_mol * self.molar_heat_capacity


def rk4_step(func: Callable[[float, float], float], t: float, y: float, dt: float) -> float:
    k1 = func(t, y)
    k2 = func(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = func(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = func(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

# função para integrar a temperatura média
def integrate_temperature(
    t_end: float,
    dt: float,
    T0: float,
    flux_in: Callable[[float], float],
    params: ClimateParams,
):
    """
    Modelo 0D de temperatura média (Pinotti; TCC RGGFarias, eq. 3.7).

    dT/dt = (1-A) F(t) / (4 χ c_p) − ε σ T^4 / (χ c_p),

    com χ densidade colunar molar [mol m^-2], c_p capacidade térmica molar
    [J mol^-1 K^-1], e F(t) o fluxo radiativo total na órbita (W m^-2), como
    na eq. 3.2 do mesmo trabalho (análogo à “constante solar” instantânea).
    """
    C = params.heat_capacity_area
    A = params.albedo
    eps = params.emissivity

    def rhs(t: float, T: float) -> float:
        incoming = (1.0 - A) * flux_in(t) / 4.0
        outgoing = eps * SIGMA_SB * T**4
        return (incoming - outgoing) / C

    ts = []
    temps = []

    t = 0.0
    T = T0

    while t <= t_end:
        ts.append(t)
        temps.append(T)
        T = rk4_step(rhs, t, T, dt)
        t += dt

    return ts, temps