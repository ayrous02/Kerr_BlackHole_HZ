from black_sun_scenario import BlackSunScenario, REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL, equilibrium_temperature_K
from lsdplus_flux import PAPER_HZ_ANCHORS
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Callable
from constants import SIGMA_SB
from climate_params import ClimateParams


# pegar o output devolvido no Teste 2
print("Dados Terra no paper: \n\n", PAPER_HZ_ANCHORS[1])
print(" ----------------------------- ")
print("ALVO: ", PAPER_HZ_ANCHORS[1][0])
print("MASSA SOLAR BURACO NEGRO: ", REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL)
print("SPIN: ", PAPER_HZ_ANCHORS[1][2])
print("R_ORB: ", PAPER_HZ_ANCHORS[1][1])
print(" ----------------------------- ")

# cria o objeto do cenario do artigo
scenario_blackSun = BlackSunScenario(
    mass_solar = REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL,
    spin_a = PAPER_HZ_ANCHORS[1][2],
    r_orb_geom = PAPER_HZ_ANCHORS[1][1],
)

# vamos usar o fluxo bolometrico calibrado. interpolações, dados mockados
phi_flux = scenario_blackSun.flux_W_m2(mode="hz_calibrated")
print("\n\nFluxo bolometrico calibrado: ", phi_flux)

# função F(t) para o integrador de clima chamar a cada passo de tempo:
flux_in = scenario_blackSun.flux_in(shadow_backend="carter")

#temperatura de equilibrio
T_eq = equilibrium_temperature_K(phi_flux)

# parametros iniciais do clima
