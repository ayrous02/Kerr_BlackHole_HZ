from __future__ import annotations
from black_sun_scenario import BlackSunScenario, REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL
from lsdplus_flux import PAPER_HZ_ANCHORS

from dataclasses import dataclass
import math

from typing import Callable
from constants import SIGMA_SB, G
from climate_params import ClimateParams

import csv

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


# interpretações
def classify_temperature_K(T: float) -> str:
    """
    Classificação simples da temperatura média global.
    Não é uma avaliação completa de habitabilidade, apenas um diagnóstico térmico.
    """
    if T < 250.0:
        return "muito frio"
    if T < 273.15:
        return "frio, abaixo do ponto de congelamento da água"
    if T <= 310.0:
        return "temperado, compatível com água líquida em princípio"
    if T <= 373.15:
        return "quente, mas ainda abaixo da ebulição da água a 1 atm"
    return "muito quente, acima da faixa clássica de água líquida"


def compute_radiative_temperatures(phi_flux: float) -> dict[str, float]:
    """
    Calcula a temperatura de cada caso radiativo para o fluxo da âncora atual.
    Retorna um dicionário:
    {
        "Blackbody Ideal": T,
        "earth_effective_temperature": T,
        ...
    }
    """
    temperatures = {}

    for case in RADIATIVE_CASES:
        T = equilibrium_temperature_general(
            phi_flux,
            albedo=case.albedo,
            emissivity=case.emissivity,
            redistribution_factor=case.redistribution_factor
        )

        temperatures[case.name] = T

    return temperatures

def interpret_anchor(
    anchor_name: str,
    phi_flux: float,
    scenario: BlackSunScenario,
    T_bakala: float,
    radiative_temperatures: dict[str, float],
    reference_planet: PlanetParams = EARTH_LIKE,
) -> None:
    roche_margin = scenario.roche_margin_ratio(
        reference_planet.mass_kg,
        reference_planet.radius_m
    )

    survives_tidal = scenario.survives_tidal_approximation(
        reference_planet.mass_kg,
        reference_planet.radius_m
    )

    print("\n================= INTERPRETAÇÃO DA ÂNCORA =================")

    print(f"Âncora analisada: {anchor_name}")
    print(f"Fluxo recebido do CMB blueshiftado: {phi_flux:.2f} W/m²")
    print(f"Temperatura de corpo negro ideal: {T_bakala:.2f} K ({T_bakala - 273.15:.2f} °C)")
    print(f"Classificação térmica do corpo negro: {classify_temperature_K(T_bakala)}")

    if survives_tidal:
        print(
            f"Marés: o planeta sobrevive na aproximação usada "
            f"(margem de Roche = {roche_margin:.4f})."
        )
    else:
        print(
            f"Marés: o planeta NÃO sobrevive na aproximação usada "
            f"(margem de Roche = {roche_margin:.4f})."
        )

    earth_eff = radiative_temperatures.get("earth_effective_temperature")
    earth_greenhouse = radiative_temperatures.get("earth_surface_greenhouse_effective")
    mars_like = radiative_temperatures.get("mars_like_radiative_case")
    venus_like = radiative_temperatures.get("venus_greenhouse_effective")

    if earth_eff is not None:
        print(
            f"Com albedo terrestre e sem efeito estufa efetivo, a temperatura cai para "
            f"{earth_eff:.2f} K ({earth_eff - 273.15:.2f} °C), indicando o efeito resfriador "
            f"da reflexão planetária."
        )

    if earth_greenhouse is not None:
        print(
            f"Com uma emissividade efetiva tipo Terra, a temperatura sobe para "
            f"{earth_greenhouse:.2f} K ({earth_greenhouse - 273.15:.2f} °C), mostrando que "
            f"a retenção radiativa pode deslocar o planeta para uma faixa mais quente."
        )

    if mars_like is not None:
        print(
            f"O caso radiativo tipo Marte resulta em {mars_like:.2f} K "
            f"({mars_like - 273.15:.2f} °C), sugerindo uma atmosfera pouco eficiente em reter calor."
        )

    if venus_like is not None:
        print(
            f"O caso tipo Vênus resulta em {venus_like:.2f} K "
            f"({venus_like - 273.15:.2f} °C), mostrando que baixa emissividade/forte efeito estufa "
            f"pode superaquecer o planeta mesmo quando o fluxo orbital não é extremo."
        )

    
    print("\nSíntese específica:")
    if anchor_name == "mars_cold":
        print(
            "Esta âncora representa a borda fria da zona habitável do modelo. "
            "Mesmo com o CMB blueshiftado fornecendo energia ao planeta, a temperatura de corpo negro "
            "fica muito abaixo de 273 K. Isso sugere que, sem um efeito estufa forte, o planeta tenderia "
            "a permanecer congelado. O caso tipo Vênus mostra que uma atmosfera extremamente opaca ao "
            "infravermelho poderia aquecer muito o planeta, mas esse é um cenário limite e ainda não "
            "representa uma atmosfera fisicamente modelada."
        )

    elif anchor_name == "earth":
        print(
            "Esta âncora reproduz o ponto intermediário do paper, com fluxo equivalente ao recebido pela Terra. "
            "O corpo negro ideal fica em torno de 278 K, próximo ao valor de referência de Bakala. "
            "Quando se adiciona albedo terrestre, a temperatura cai abaixo do congelamento, mas com uma "
            "emissividade efetiva tipo Terra ela sobe para aproximadamente 288 K. Isso mostra que esta "
            "âncora é a mais parecida com uma situação terrestre simplificada."
        )

    elif anchor_name == "venus_hot":
        print(
            "Esta âncora representa a borda quente da zona habitável do modelo. "
            "O corpo negro ideal já fica acima de 327 K, indicando um ambiente quente. "
            "Com albedo terrestre e sem efeito estufa efetivo, a temperatura ainda permanece moderada, "
            "mas com retenção radiativa tipo Terra ela sobe bastante. O caso tipo Vênus entra em regime "
            "extremamente quente, sugerindo risco de superaquecimento caso a atmosfera seja muito eficiente "
            "em reter radiação infravermelha."
        )

    else:
        print(
            "Nesta etapa, a âncora define quanta energia chega ao planeta; os casos radiativos mostram "
            "como albedo e emissividade modificam a temperatura média. Ainda não é um modelo climático "
            "completo, pois não inclui dinâmica atmosférica, composição química detalhada, pressão, "
            "circulação ou variação temporal."
        )



def print_planet_catalog() -> None:
    print("================= PLANETAS USADOS ===========================")
    print(f"{'Planeta':<15} {'Massa [kg]':>14} {'Raio [m]':>12} {'rho [kg/m3]':>14} {'g [m/s2]':>12} {'v_escape [km/s]':>16}")

    for planet in PLANETS:
        print(
            f"{planet.name:<15} "
            f"{planet.mass_kg:>14.3e} "
            f"{planet.radius_m:>12.3e} "
            f"{density_kg_m3(planet):>14.2f} "
            f"{surface_gravity_m_s2(planet):>12.2f} "
            f"{espace_velocity_m_s(planet) / 1000:>16.2f}"
        )


def print_radiative_cases_table(
    phi_flux: float,
    T_bakala: float,
    radiative_temperatures: dict[str, float],
) -> None:
    print(f"\nCasos radiativos recebendo {phi_flux:.2f} W/m²:")
    print(
        f"{'Caso':<38} "
        f"{'Planeta':<12} "
        f"{'A':>6} "
        f"{'eps':>6} "
        f"{'T [K]':>10} "
        f"{'T [°C]':>10} "
        f"{'Delta Bakala':>14}"
    )

    for case in RADIATIVE_CASES:
        planet = PLANETS_BY_NAME[case.planet_name]
        T = radiative_temperatures[case.name]

        delta_txt = "-"
        if case.name == "Blackbody Ideal":
            delta_txt = f"{T - T_bakala:.3f}"

        print(
            f"{case.name:<38} "
            f"{planet.name:<12} "
            f"{case.albedo:>6.3f} "
            f"{case.emissivity:>6.3f} "
            f"{T:>10.2f} "
            f"{T - 273.15:>10.2f} "
            f"{delta_txt:>14}"
        )


def print_albedo_sweep(phi_flux: float) -> None:
    print("\nVariação simples de albedo, sem atmosfera:")
    print(f"{'Albedo':>8} {'T global [K]':>14} {'T hemisf. [K]':>14} {'T substelar [K]':>16}")

    for albedo in [0.0, 0.1, 0.2, 0.3]:
        T_global = equilibrium_temperature_general(
            phi_flux,
            albedo=albedo,
            emissivity=1.0,
            redistribution_factor=4.0
        )

        T_hotsky = equilibrium_temperature_general(
            phi_flux,
            albedo=albedo,
            emissivity=1.0,
            redistribution_factor=2.0
        )

        T_substellar = equilibrium_temperature_general(
            phi_flux,
            albedo=albedo,
            emissivity=1.0,
            redistribution_factor=1.0
        )

        print(
            f"{albedo:>8.2f} "
            f"{T_global:>14.2f} "
            f"{T_hotsky:>14.2f} "
            f"{T_substellar:>16.2f}"
        )


def main() -> None:
    print_planet_catalog()

    results = []  # fica FORA do loop

    print("\n================= ÂNCORAS DO PAPER ===========================")

    for name, r_orbit, spin, phi_ref, g_ref in PAPER_HZ_ANCHORS:

        scenario_black_sun = BlackSunScenario(
            mass_solar=REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL,
            spin_a=spin,
            r_orb_geom=r_orbit,
        )

        phi_flux = scenario_black_sun.flux_W_m2(mode="hz_calibrated")
        T_bakala = scenario_black_sun.equilibrium_temperature_K(phi_flux)

        roche_margin = scenario_black_sun.roche_margin_ratio(
            EARTH_LIKE.mass_kg,
            EARTH_LIKE.radius_m
        )

        survives_tidal = scenario_black_sun.survives_tidal_approximation(
            EARTH_LIKE.mass_kg,
            EARTH_LIKE.radius_m
        )

        print("\n--------------------------------------------------------------")
        print(f"Âncora: {name}")
        print(f"R_orb [GM/c²]: {r_orbit:.6f}")
        print(f"Spin a*: {spin:.12f}")
        print(f"Fluxo CMB calibrado [W/m²]: {phi_flux:.2f} | referência paper: {phi_ref:.2f}")
        print(f"g_max referência paper: {g_ref:.2f}")
        print(f"Raio orbital [m]: {scenario_black_sun.orbit_radius_si_m():.3e}")
        print(f"Dilatação temporal gamma: {scenario_black_sun.time_dilation_gamma():.2f}")
        print(f"Margem de Roche [R_orb/R_t]: {roche_margin:.4f}")
        print(f"Sobrevive à aproximação tidal? {survives_tidal}")
        print(f"T_Bakala corpo negro [K]: {T_bakala:.2f}")
        print(f"T_Bakala corpo negro [°C]: {T_bakala - 273.15:.2f}")

        print_albedo_sweep(phi_flux)

        radiative_temperatures = compute_radiative_temperatures(phi_flux)

        print_radiative_cases_table(
            phi_flux=phi_flux,
            T_bakala=T_bakala,
            radiative_temperatures=radiative_temperatures,
        )

        interpret_anchor(
            anchor_name=name,
            phi_flux=phi_flux,
            scenario=scenario_black_sun,
            T_bakala=T_bakala,
            radiative_temperatures=radiative_temperatures,
        )

        for case in RADIATIVE_CASES:
            planet = PLANETS_BY_NAME[case.planet_name]
            T = radiative_temperatures[case.name]

            results.append({
                "anchor": name,
                "r_orb_GM_c2": r_orbit,
                "spin": spin,
                "phi_flux_W_m2": phi_flux,
                "phi_ref_W_m2": phi_ref,
                "g_ref": g_ref,
                "orbit_radius_m": scenario_black_sun.orbit_radius_si_m(),
                "time_dilation_gamma": scenario_black_sun.time_dilation_gamma(),
                "roche_margin": roche_margin,
                "survives_tidal": survives_tidal,
                "T_bakala_K": T_bakala,
                "T_bakala_C": T_bakala - 273.15,
                "radiative_case": case.name,
                "planet": planet.name,
                "planet_mass_kg": planet.mass_kg,
                "planet_radius_m": planet.radius_m,
                "planet_density_kg_m3": density_kg_m3(planet),
                "planet_gravity_m_s2": surface_gravity_m_s2(planet),
                "planet_escape_velocity_km_s": espace_velocity_m_s(planet) / 1000,
                "albedo": case.albedo,
                "emissivity": case.emissivity,
                "redistribution_factor": case.redistribution_factor,
                "T_K": T,
                "T_C": T - 273.15,
                "thermal_classification": classify_temperature_K(T),
            })

    # fica FORA do loop das âncoras
    with open("planet_modelling_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print("\nArquivo salvo: planet_modelling_results.csv")
    print(f"Total de cenários salvos: {len(results)}")

if __name__ == "__main__":
    main()
