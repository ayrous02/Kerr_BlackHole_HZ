from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np


AU_M = 1.495978707e11


def first_non_null(series: pd.Series):
    values = series.dropna()
    if values.empty:
        return np.nan
    return values.iloc[0]


def deduplicate_by_keys(df: pd.DataFrame, keys: list[str], cols: list[str]) -> pd.DataFrame:
    """
    Agrupa por chaves e pega o primeiro valor não-nulo de cada coluna.
    Útil porque alguns CSVs têm duplicatas herdadas de radiative_case.
    """
    return (
        df[keys + cols]
        .groupby(keys, as_index=False)
        .agg(first_non_null)
    )


def main() -> None:
    selected_csv = Path("psg_selected_cases.csv")
    planet_csv = Path("planet_modelling_results.csv")
    candidate_csv = Path("atmosphere_candidate_screening.csv")
    physical_csv = Path("atmosphere_physical_parameters.csv")
    climate_csv = Path("atmospheric_climate_results.csv")
    output_csv = Path("psg_ready_cases.csv")

    selected = pd.read_csv(selected_csv)
    planet = pd.read_csv(planet_csv)
    candidate = pd.read_csv(candidate_csv)
    physical = pd.read_csv(physical_csv)
    climate = pd.read_csv(climate_csv)

    # 1) Dados físicos do planeta e do cenário Bakala
    planet_keys = ["anchor", "planet"]
    planet_cols = [
        "r_orb_GM_c2",
        "spin",
        "phi_flux_W_m2",
        "phi_ref_W_m2",
        "g_ref",
        "orbit_radius_m",
        "time_dilation_gamma",
        "roche_margin",
        "survives_tidal",
        "T_bakala_K",
        "T_bakala_C",
        "planet_mass_kg",
        "planet_radius_m",
        "planet_density_kg_m3",
        "planet_gravity_m_s2",
        "planet_escape_velocity_km_s",
    ]

    planet_unique = deduplicate_by_keys(
        planet,
        keys=planet_keys,
        cols=planet_cols,
    )

    # 2) Dados da família atmosférica candidata
    candidate_keys = ["anchor", "planet", "candidate_atmosphere"]
    candidate_cols = [
        "category",
        "priority",
        "required_gases",
        "optional_gases",
        "retained_gases",
        "scientific_reason",
        "comment",
    ]

    candidate_unique = deduplicate_by_keys(
        candidate,
        keys=candidate_keys,
        cols=candidate_cols,
    )

    # 3) Dados físicos da mistura atmosférica
    physical_keys = ["anchor", "planet", "candidate_atmosphere", "mixture_name"]
    physical_cols = [
        "mole_fractions",
        "surface_pressure_bar",
        "mean_molar_mass_g_mol",
        "cp_mix_J_mol_K",
        "column_density_mol_m2",
        "albedo_assumed",
        "emissivity_assumed",
        "mixture_notes",
    ]

    physical_unique = deduplicate_by_keys(
        physical,
        keys=physical_keys,
        cols=physical_cols,
    )

    # 4) Dados climáticos finais
    climate_keys = ["anchor", "planet", "candidate_atmosphere", "mixture_name"]
    climate_cols = [
        "redistribution_factor_atm",
        "T_eq_atm_K",
        "T_eq_atm_C",
        "climate_classification",
        "column_heat_capacity_J_m2_K",
        "radiative_timescale_seconds",
        "radiative_timescale_days",
        "radiative_timescale_years",
    ]

    climate_unique = deduplicate_by_keys(
        climate,
        keys=climate_keys,
        cols=climate_cols,
    )

    # 5) Merge geral
    df = selected.copy()

    df = df.merge(
        planet_unique,
        on=["anchor", "planet"],
        how="left",
        suffixes=("", "_planet_model"),
    )

    df = df.merge(
        candidate_unique,
        on=["anchor", "planet", "candidate_atmosphere"],
        how="left",
    )

    df = df.merge(
        physical_unique,
        on=["anchor", "planet", "candidate_atmosphere", "mixture_name"],
        how="left",
        suffixes=("", "_physical"),
    )

    df = df.merge(
        climate_unique,
        on=["anchor", "planet", "candidate_atmosphere", "mixture_name"],
        how="left",
        suffixes=("", "_climate"),
    )

    # 6) Corrigir/derivar campos úteis para PSG
    # Se psg_selected_cases tinha planet_radius_m vazio, usa o do planet_modelling.
    if "planet_radius_m_planet_model" in df.columns:
        df["planet_radius_m"] = df["planet_radius_m"].combine_first(
            df["planet_radius_m_planet_model"]
        )

    if "planet_gravity_m_s2_planet_model" in df.columns:
        df["planet_gravity_m_s2"] = df["planet_gravity_m_s2"].combine_first(
            df["planet_gravity_m_s2_planet_model"]
        )

    df["planet_diameter_km"] = 2.0 * df["planet_radius_m"] / 1000.0
    df["orbit_radius_AU"] = df["orbit_radius_m"] / AU_M

    # Temperatura para PSG: usar a temperatura climática final se existir.
    df["T_surface_K_for_PSG"] = df["T_surface_K_for_PSG"].combine_first(
        df["T_eq_atm_K"]
    )

    df["T_surface_C_for_PSG"] = df["T_surface_C_for_PSG"].combine_first(
        df["T_eq_atm_C"]
    )

    # Pressão em Pa, caso ainda não esteja preenchida
    if "surface_pressure_Pa" not in df.columns:
        df["surface_pressure_Pa"] = df["surface_pressure_bar"] * 1e5
    else:
        df["surface_pressure_Pa"] = df["surface_pressure_Pa"].combine_first(
            df["surface_pressure_bar"] * 1e5
        )

    # 7) Selecionar colunas finais mais úteis
    final_cols = [
        "case_id",
        "anchor",
        "planet",
        "candidate_atmosphere",
        "mixture_name",
        "climate_priority",

        # Bakala / Black Sun
        "r_orb_GM_c2",
        "spin",
        "phi_flux_W_m2",
        "g_ref",
        "orbit_radius_m",
        "orbit_radius_AU",
        "time_dilation_gamma",
        "roche_margin",
        "survives_tidal",
        "T_bakala_K",

        # Planeta
        "planet_mass_kg",
        "planet_radius_m",
        "planet_diameter_km",
        "planet_density_kg_m3",
        "planet_gravity_m_s2",
        "planet_escape_velocity_km_s",

        # Atmosfera
        "surface_pressure_bar",
        "surface_pressure_Pa",
        "mean_molar_mass_g_mol",
        "cp_mix_J_mol_K",
        "column_density_mol_m2",
        "column_heat_capacity_J_m2_K",
        "albedo_assumed",
        "emissivity_assumed",
        "mole_fractions",
        "molecules_for_psg",

        # Clima
        "T_surface_K_for_PSG",
        "T_surface_C_for_PSG",
        "T_eq_atm_K",
        "T_eq_atm_C",
        "climate_classification",
        "radiative_timescale_days",

        # Justificativas
        "category",
        "required_gases",
        "optional_gases",
        "retained_gases",
        "scientific_reason",
        "comment",
        "notes",
    ]

    existing_final_cols = [col for col in final_cols if col in df.columns]
    df = df[existing_final_cols].copy()

    df.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de casos: {len(df)}")

    print(
        df[
            [
                "case_id",
                "planet",
                "planet_diameter_km",
                "planet_gravity_m_s2",
                "T_surface_K_for_PSG",
                "surface_pressure_bar",
                "molecules_for_psg",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()