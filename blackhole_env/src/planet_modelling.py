from __future__ import annotations
from black_sun_scenario import BlackSunScenario, REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL
from lsdplus_flux import PAPER_HZ_ANCHORS

from dataclasses import dataclass
import math

from typing import Callable
from constants import SIGMA_SB, G
from climate_params import ClimateParams


@dataclass(frozen=True)
class PlanetParams:
    name: str
    mass_kg: float
    radius_m: float
    

@dataclass(frozen=True)
class PlanetRadiativeCase:
    name: str
    planet_name: str
    albedo: float=0.0 #blackbody
    emissivity: float = 1.0 #blackbody
    redistribution_factor: float = 4.0


def density_kg_m3(planet: PlanetParams) -> float:
    volume = (4.0 / 3.0) * math.pi * planet.radius_m**3
    return planet.mass_kg / volume

def surface_gravity_m_s2(planet: PlanetParams) -> float:
    return G * planet.mass_kg / planet.radius_m**2

def espace_velocity_m_s(planet: PlanetParams) -> float:
    return math.sqrt(2.0 * G * planet.mass_kg / planet.radius_m)

def equilibrium_temperature_general(
    phi_W_m2: float,
    albedo: float = 0.0,
    emissivity: float = 1.0,
    redistribution_factor: float = 4.0
) -> float:
    """
    Temperatura de equilíbrio generalizada.

    redistribution_factor = 4 -> média global, como no Bakala.
    redistribution_factor = 2 -> média no hemisfério iluminado.
    redistribution_factor = 1 -> ponto/região diretamente aquecida.
    """
    absorvido = (1.0 - albedo) * phi_W_m2
    emitido = redistribution_factor * emissivity * SIGMA_SB
    return (absorvido / emitido)**0.25

def print_planet_all(scenario: BlackSunScenario, planet: PlanetParams, phi_W_m2: float):
    print("-----------------------------")
    print("PLANET: ", planet.name)
    print("Massa: [kg]", planet.mass_kg)
    print("RAIO [m]", planet.radius_m)
    print("DENSIDADE [kg/m3]", density_kg_m3(planet))
    print("GRAVIDADE [m/s2]", surface_gravity_m_s2(planet))
#    print("VELOCIDADE DE ESCAPE [m/s]", espace_velocity_m_s(planet))
    print("VELOCIDADE DE ESCAPE [km/s]", espace_velocity_m_s(planet) / 1000)

    print("\n ============= CENÁRIO DE CORPO NEGRO IDEAL =================")
    print("MASSA BURACO NEGRO [Msol]: ", scenario.mass_solar)
    print("SPIN BURACO NEGRO [adimensional]: ", scenario.spin_a)
    print("R_ORB [GM/c2]: ", scenario.r_orb_geom)
    print("RAIO ORBITAL [m]: ", scenario.orbit_radius_si_m())
    print("FLUXO BOLOMÉTRICO CMB CALIBRADO: ", phi_W_m2)
    print("FATOR DE DILATAÇÃO TEMPORAL: ", scenario.time_dilation_gamma())

    print("\n ============= MARÉS =================")
    print("RAIO DE DISRUPÇÃO TIDAL APPROX. [m]: ", scenario.tidal_disruption_radius_m(planet.mass_kg, planet.radius_m))
    print("MARGEM DE ROCHE [R_orb/R_t]: ", scenario.roche_margin_ratio(planet.mass_kg, planet.radius_m))
    print("Sobrevive a aproximação tidal? ", scenario.survives_tidal_approximation(planet.mass_kg, planet.radius_m))

    print("\n ============= TEMPERATURAS PRÉ ATMOSFERA =================")
    T_bakala = scenario.equilibrium_temperature_K(phi_W_m2)
    print("Albedo : 0.0")
    print("Emissividade 1.0")
    print("Redistribuição global: ", T_bakala)
    print("T_Bakala [ºC]: ", T_bakala - 273.15)

    #explorar variação de albedo (não mais corpo negro)
    for albedo in [0.0, 0.1, 0.2, 0.3]:
        T_global = equilibrium_temperature_general(
            phi_W_m2,
            albedo = albedo,
            emissivity = 1.0,
            redistribution_factor = 4.0

        )
        T_hotsky = equilibrium_temperature_general(
            phi_W_m2,
            albedo=albedo,
            emissivity=1.0,
            redistribution_factor=2.0
        )
        T_substellar = equilibrium_temperature_general(
            phi_W_m2,
            albedo=albedo,
            emissivity=1.0,
            redistribution_factor=1.0
        )

        print("Albedo variado: ", albedo)
        print("T_media global: ", T_global)
        print("T hemisfério iluminado: ", T_hotsky)
        print("T substelar local: ", T_substellar)
        print("\n")




# Valores tirados do artigo de Bakala e da vida real
BLACKBODY_CASE = PlanetRadiativeCase(
    name="Blackbody Ideal",
    planet_name= "Earth-like",
    albedo=0.0,
    emissivity=1.0,
    redistribution_factor=4.0,
)

EARTH_LIKE = PlanetParams(
    name="Earth-like",
    mass_kg=5.97e24, #calculado pelo paper
    radius_m=6.37e6, #calculado pelo paper
    #albedo=0.306, #albedo da Terra
    #emissivity=0.61, #emissividade da Terra
)

EARTH_EFFECTIVE_TEMPERATURE_CASE = PlanetRadiativeCase(
    name="earth_effective_temperature",
    planet_name="Earth-like",
    albedo=0.306,
    emissivity=1.0,
    redistribution_factor=4.0,
)

EARTH_SURFACE_GREENHOUSE_CASE = PlanetRadiativeCase(
    name="earth_surface_greenhouse_effective",
    planet_name="Earth-like",
    albedo=0.306,
    emissivity=0.61,
    redistribution_factor=4.0,
)

MARS_LIKE = PlanetParams(
    name="Mars-like",
    mass_kg=0.642e24,
    radius_m=3.390e6,
)
MARS_LIKE_CASE = PlanetRadiativeCase(
    name="mars_like_radiative_case",
    planet_name="Mars-like",
    albedo= 0.25,
    emissivity=0.95,
    redistribution_factor=4.0,
)

VENUS_LIKE = PlanetParams(
    name="Venus-like",
    mass_kg=4.867e24,
    radius_m=6.052e6,
)
VENUS_GREENHOUSE_EFFECTIVE_CASE = PlanetRadiativeCase(
    name="venus_greenhouse_effective",
    planet_name="Venus-like",
    albedo= 0.756,
    emissivity=0.013,
    redistribution_factor=4.0,
)

PLANETS = [
    EARTH_LIKE,
    MARS_LIKE,
    VENUS_LIKE,
]

RADIATIVE_CASES = [
    BLACKBODY_CASE,
    EARTH_EFFECTIVE_TEMPERATURE_CASE,
    EARTH_SURFACE_GREENHOUSE_CASE,
    MARS_LIKE_CASE,
    VENUS_GREENHOUSE_EFFECTIVE_CASE,
]
PLANETS_BY_NAME = {
    planet.name: planet
    for planet in PLANETS
}


def main() -> None:

    #começando pela Terra
    target_T = PAPER_HZ_ANCHORS[1]
    name = target_T[0]
    r_orbit = target_T[1]
    spin = target_T[2]
    phi_ref = target_T[3]
    g_ref = target_T[4]


    print("================= ANCORA PRINCIPAL DO PAPER =========================== ")
    print("ALVO 01: ", name)
    print("R_ORB: ", r_orbit)
    print("SPIN: ", spin)
    print("Fluxo de referência: ", phi_ref)
    print("g_max de referencia: ", g_ref)

    scenario_BlackSun = BlackSunScenario(
        mass_solar=REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL,
        spin_a=spin,
        r_orb_geom=r_orbit,
    )

    phi_flux = scenario_BlackSun.flux_W_m2(mode="hz_calibrated")

    print_planet_all(
        scenario_BlackSun,
        EARTH_LIKE,
        phi_flux
    )

    T_bakala = scenario_BlackSun.equilibrium_temperature_K(phi_flux)

    print("\n")
    print("================= CASOS RADIATIVOS =========================== ")

    for case in RADIATIVE_CASES:
        planet = PLANETS_BY_NAME[case.planet_name]
        T = equilibrium_temperature_general(
            phi_flux,
            albedo = case.albedo,
            emissivity = case.emissivity,
            redistribution_factor = case.redistribution_factor
        )

        print("CASO: ", case.name)
        print("Albedo: ", case.albedo)
        print("Emissividade: ", case.emissivity)
        print("T [K]: ", T)
    

        densidade = density_kg_m3(planet)
        gravidade = surface_gravity_m_s2(planet)
        velocidade_escape = espace_velocity_m_s(planet)
        print("Densidade: ", densidade)
        print("Gravidade: ", gravidade)
        print("Velocidade de espace: ", velocidade_escape)

        if case.name == "blackbody_ideal":
            print("Diferença para T_Bakala [K]:", T - T_bakala)
    
        print("\n\n")


        

if __name__ == "__main__":
    main()
