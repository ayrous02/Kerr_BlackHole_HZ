from __future__ import annotations

import csv
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from constants import SIGMA_SB
from black_sun_scenario import BlackSunScenario, REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL
from lsdplus_flux import PAPER_HZ_ANCHORS
from atmosphere_cases import ATMOSPHERE_CASES


def equilibrium_temperature(phi, albedo, emissivity, redistribution_factor):
    return ((1 - albedo) * phi / (redistribution_factor * emissivity * SIGMA_SB)) ** 0.25


def classify_temperature(T):
    if T < 240:
        return "muito frio"
    if 240 <= T < 273.15:
        return "frio/marginal"
    if 273.15 <= T <= 373.15:
        return "faixa água líquida"
    if 373.15 < T <= 450:
        return "quente/marginal"
    return "muito quente"


def main():
    anchor = PAPER_HZ_ANCHORS[1]  # earth
    name, r_orb, spin, phi_ref, g_ref = anchor

    scenario = BlackSunScenario(
        mass_solar=REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL,
        spin_a=spin,
        r_orb_geom=r_orb,
    )

    phi = scenario.flux_W_m2(mode="hz_calibrated")

    rows = []

    for atm in ATMOSPHERE_CASES:
        for A in atm["albedos"]:
            for eps in atm["emissivities"]:
                for f in atm["redistribution"]:
                    T = equilibrium_temperature(phi, A, eps, f)
                    rows.append({
                        "anchor": name,
                        "atmosphere_case": atm["name"],
                        "composition": atm["composition"],
                        "phi_W_m2": phi,
                        "albedo": A,
                        "emissivity": eps,
                        "redistribution_factor": f,
                        "T_K": T,
                        "T_C": T - 273.15,
                        "classification": classify_temperature(T),
                    })

    out = Path("climate_grid_results.csv")
    with out.open("w", newline="", encoding="utf-8") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Arquivo salvo: {out.resolve()}")
    print(f"Total de cenários: {len(rows)}")


if __name__ == "__main__":
    main()