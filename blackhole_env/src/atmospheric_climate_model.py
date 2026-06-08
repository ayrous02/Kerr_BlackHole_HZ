from __future__ import annotations

import pandas as pd
import numpy as np

from constants import SIGMA_SB


def equilibrium_temperature_atmosphere(
    phi_flux_W_m2: float,
    albedo: float,
    emissivity: float,
    redistribution_factor: float = 4.0,
) -> float:
    """
    Temperatura de equilíbrio com parâmetros atmosféricos efetivos.

    redistribution_factor = 4 -> média global
    redistribution_factor = 2 -> hemisfério aquecido
    redistribution_factor = 1 -> região/subponto aquecido
    """
    absorbed_flux = (1.0 - albedo) * phi_flux_W_m2 / redistribution_factor
    emitted_factor = emissivity * SIGMA_SB
    return (absorbed_flux / emitted_factor) ** 0.25


def column_heat_capacity_J_m2_K(
    column_density_mol_m2: float,
    cp_mix_J_mol_K: float,
) -> float:
    """
    Capacidade térmica da coluna atmosférica por unidade de área.
    """
    return column_density_mol_m2 * cp_mix_J_mol_K


def radiative_timescale_seconds(
    column_heat_capacity: float,
    equilibrium_temperature_K: float,
    emissivity: float,
) -> float:
    """
    Escala de tempo radiativa aproximada.

    tau ~ C / (4 epsilon sigma T^3)
    """
    return column_heat_capacity / (
        4.0 * emissivity * SIGMA_SB * equilibrium_temperature_K**3
    )


def classify_climate(T_K: float) -> str:
    if T_K < 240.0:
        return "muito frio"
    if T_K < 273.15:
        return "frio/marginal"
    if T_K <= 310.0:
        return "temperado"
    if T_K <= 373.15:
        return "quente/marginal"
    return "muito quente"


def main() -> None:
    input_csv = "atmosphere_physical_parameters.csv"
    output_csv = "atmospheric_climate_results.csv"

    df = pd.read_csv(input_csv)

    rows = []

    for _, row in df.iterrows():
        redistribution_factor = 4.0

        T_eq_atm = equilibrium_temperature_atmosphere(
            phi_flux_W_m2=row["phi_flux_W_m2"],
            albedo=row["albedo_assumed"],
            emissivity=row["emissivity_assumed"],
            redistribution_factor=redistribution_factor,
        )

        C_col = column_heat_capacity_J_m2_K(
            column_density_mol_m2=row["column_density_mol_m2"],
            cp_mix_J_mol_K=row["cp_mix_J_mol_K"],
        )

        tau_s = radiative_timescale_seconds(
            column_heat_capacity=C_col,
            equilibrium_temperature_K=T_eq_atm,
            emissivity=row["emissivity_assumed"],
        )

        rows.append({
            **row.to_dict(),
            "redistribution_factor_atm": redistribution_factor,
            "T_eq_atm_K": T_eq_atm,
            "T_eq_atm_C": T_eq_atm - 273.15,
            "climate_classification": classify_climate(T_eq_atm),
            "column_heat_capacity_J_m2_K": C_col,
            "radiative_timescale_seconds": tau_s,
            "radiative_timescale_days": tau_s / 86400.0,
            "radiative_timescale_years": tau_s / (86400.0 * 365.25),
        })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de linhas: {len(out)}")

    print(
        out[
            [
                "anchor",
                "planet",
                "candidate_atmosphere",
                "mixture_name",
                "T_eq_atm_K",
                "T_eq_atm_C",
                "climate_classification",
                "column_heat_capacity_J_m2_K",
                "radiative_timescale_days",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()