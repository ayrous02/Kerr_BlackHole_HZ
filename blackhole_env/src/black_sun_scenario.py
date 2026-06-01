from __future__ import annotations

"""
Cenário “Black Sun” (Bakala et al. 2020): planeta em órbita circular corotante
quase na ISCO de um Kerr quase extremo, aquecido pelo CMB blueshitzado/amplificado.

Fornece Φ em W m⁻² (calibração HZ do paper), verificação grosso modo de marés
(Rees 1988), fator Γ de dilatação temporal e um callable flux_in(t) para o
modelo atmosférico 0D (planet_climate).
"""

from dataclasses import dataclass
from typing import Callable, Optional

from bakala_cmb_flux import (
    equilibrium_temperature_black_body,
    hz_calibrated_flux_W_m2,
    make_constant_flux_callable,
    suggested_spin_for_hz_orbit,
    time_dilation_gamma_kepler,
    integrate_cmb_flux_mollweide_naive,
    integrate_cmb_flux_mollweide_with_shadow,
)
from constants import C, G, M_SUN
from kerr_orbit import KerrOrbit


# Bakala et al. (2020), Sec. 5 — ordem de grandeza para Terra na HZ.
REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL = 1.63e8


@dataclass
class BlackSunScenario:
    """
    mass_solar: massa do BH em massas solares.
    spin_a: momento angular adimensional a ∈ [0, 1).
    r_orb_geom: raio orbital em GM/c² (órbita circular corotante).
    """

    mass_solar: float
    spin_a: float
    r_orb_geom: float

    def __post_init__(self) -> None:
        if self.mass_solar <= 0.0:
            raise ValueError("mass_solar deve ser positiva")
        if not (0.0 <= self.spin_a < 1.0):
            raise ValueError("spin_a deve estar em [0, 1)")
        KerrOrbit(spin=self.spin_a, r_orb=self.r_orb_geom).check_stable()

    @property
    def mass_kg(self) -> float:
        return self.mass_solar * M_SUN

    def orbit_radius_si_m(self) -> float:
        """r_orb em metros: r_geom GM/c²."""
        return self.r_orb_geom * G * self.mass_kg / (C * C)

    def tidal_disruption_radius_m(
        self,
        m_planet_kg: float = 5.97e24,
        r_planet_m: float = 6.37e6,
    ) -> float:
        """Eq. (26) Rees (1988): R_t ≈ (M/m)^(1/3) r_*."""
        return (self.mass_kg / m_planet_kg) ** (1.0 / 3.0) * r_planet_m

    def roche_margin_ratio(
        self,
        m_planet_kg: float = 5.97e24,
        r_planet_m: float = 6.37e6,
    ) -> float:
        """R_orbit / R_t (Rees 1988); valores ~1 são marginalmente estáveis."""
        ro = self.orbit_radius_si_m()
        rt = self.tidal_disruption_radius_m(m_planet_kg, r_planet_m)
        return ro / rt

    def survives_tidal_approximation(
        self,
        m_planet_kg: float = 5.97e24,
        r_planet_m: float = 6.37e6,
        *,
        safety_factor: float = 1.0,
    ) -> bool:
        """True se R_orbit > safety_factor × R_tidal (Newtoniano; perto da HZ ≈ limite)."""
        ro = self.orbit_radius_si_m()
        rt = self.tidal_disruption_radius_m(m_planet_kg, r_planet_m)
        return ro > safety_factor * rt

    def meets_paper_tidal_scale(self) -> bool:
        """Ordem de grandeza do paper: M ≳ 1.63×10⁸ M☉ para Terra-like na HZ."""
        return self.mass_solar >= REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL * 0.99

    def suggested_spin_hz(self) -> float:
        return suggested_spin_for_hz_orbit(self.r_orb_geom)

    def spin_hz_delta(self) -> float:
        """|a − a_sugerido(r)| — deve ser ~0 se seguires a curva HZ do artigo."""
        return abs(self.spin_a - self.suggested_spin_hz())

    def flux_W_m2(
        self,
        *,
        mode: str = "hz_calibrated",
        numeric_grid_n: int = 56,
        shadow_max_workers: Optional[int] = None,
        shadow_backend: str = "einsteinpy",
        orbit_mode: str = "free",
    ) -> float:
        """
        mode='hz_calibrated': Φ(r) dos três pontos Bakala et al. (Fig. 3).

        mode='numeric_naive': integral Mollweide g⁴ sem sombra (teste).

        mode='numeric_shadow': integral com sombra (EinsteinPy + fallback Carter por omissão).

        mode='paper_zoom': eq. (25) + sombra + zoom Brent (lsdplus_flux); reprodução do artigo.

        shadow_backend: 'einsteinpy' | 'carter' (numeric_shadow e paper_zoom).
        orbit_mode: 'free' | 'isco' (só paper_zoom; isco força r = r_ISCO(a)).
        """
        if mode == "hz_calibrated":
            return hz_calibrated_flux_W_m2(self.r_orb_geom)
        if mode == "numeric_naive":
            return integrate_cmb_flux_mollweide_naive(self.r_orb_geom, self.spin_a, n=numeric_grid_n)
        if mode == "numeric_shadow":
            return integrate_cmb_flux_mollweide_with_shadow(
                self.r_orb_geom,
                self.spin_a,
                n=numeric_grid_n,
                max_workers=shadow_max_workers,
                shadow_backend=shadow_backend,
            )
        if mode == "paper_zoom":
            from lsdplus_flux import FluxIntegrationConfig, OrbitMode, integrate_cmb_flux

            om = OrbitMode.ISCO if orbit_mode == "isco" else OrbitMode.FREE
            cfg = FluxIntegrationConfig(
                r_orb_geom=self.r_orb_geom,
                spin_a=self.spin_a,
                orbit_mode=om,
                n_grid=numeric_grid_n,
                shadow_backend=shadow_backend,
                max_workers=shadow_max_workers,
            )
            return integrate_cmb_flux(cfg).phi_W_m2
        raise ValueError(f"mode desconhecido: {mode}")

    def equilibrium_temperature_K(self, phi_W_m2: Optional[float] = None) -> float:
        """Temperatura de equilíbrio corpo negro (eq. 3 Bakala); sem albedo."""
        phi = phi_W_m2 if phi_W_m2 is not None else self.flux_W_m2()
        return equilibrium_temperature_black_body(phi)

    def time_dilation_gamma(self) -> float:
        return time_dilation_gamma_kepler(self.r_orb_geom, self.spin_a)

    def flux_in(
        self,
        *,
        mode: str = "hz_calibrated",
        numeric_grid_n: int = 56,
        shadow_max_workers: Optional[int] = None,
        shadow_backend: str = "carter",
        orbit_mode: str = "free",
    ) -> Callable[[float], float]:
        """
        flux_in(t) constante (órbita circular) em W m⁻², pronto para
        planet_climate.integrate_temperature (lembrar o fator 1/4 na absorção média).
        """
        return make_constant_flux_callable(
            self.flux_W_m2(
                mode=mode,
                numeric_grid_n=numeric_grid_n,
                shadow_max_workers=shadow_max_workers,
                shadow_backend=shadow_backend,
                orbit_mode=orbit_mode,
            )
        )

    def flux_integration_paper(
        self,
        *,
        numeric_grid_n: int = 500,
        shadow_max_workers: Optional[int] = None,
        shadow_backend: str = "carter",
        orbit_mode: str = "isco",
    ):
        """
        Resultado completo da integração (Φ, T, g_max, zoom k, offsets ISCO).
        """
        from lsdplus_flux import FluxIntegrationConfig, OrbitMode, integrate_cmb_flux

        om = OrbitMode.ISCO if orbit_mode == "isco" else OrbitMode.FREE
        return integrate_cmb_flux(
            FluxIntegrationConfig(
                r_orb_geom=self.r_orb_geom,
                spin_a=self.spin_a,
                orbit_mode=om,
                n_grid=numeric_grid_n,
                shadow_backend=shadow_backend,
                max_workers=shadow_max_workers,
            )
        )

    def risco_corotating_geom(self) -> float:
        return KerrOrbit(spin=self.spin_a, r_orb=self.r_orb_geom).risco_corotating()
