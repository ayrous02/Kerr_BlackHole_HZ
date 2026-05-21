from __future__ import annotations
from dataclasses import dataclass

from bakala_cmb_flux import HZ_CALIBRATION, hz_calibrated_flux_W_m2
from constants import SIGMA_SB


@dataclass
class BlackHoleFluxProxy:
    """
    Compatibilidade com código antigo: delega para a calibração HZ de
    Bakala et al. (2020) em bakala_cmb_flux.hz_calibrated_flux_W_m2.
    """

    r_mars: float = HZ_CALIBRATION[0][0]
    r_earth: float = HZ_CALIBRATION[1][0]
    r_venus: float = HZ_CALIBRATION[2][0]

    def flux_from_radius(self, r_orb: float) -> float:
        """
        Φ bolométrico interpolado (W m^-2), válido perto da HZ quase-extrema.
        Para cenário completo (marés, Γ, spin vs r), use BlackSunScenario.
        """
        return hz_calibrated_flux_W_m2(r_orb)

    @staticmethod
    def equilibrium_temperature_from_flux(flux: float) -> float:
        """T_eq com A = 0 e ε = 1: (1-A)F/4 = εσT^4 ⇒ T = (F/(4σ))^1/4."""
        return (flux / (4.0 * SIGMA_SB)) ** 0.25