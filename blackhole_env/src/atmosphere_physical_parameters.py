from __future__ import annotations

from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class Gas:
    name: str
    molar_mass_g_mol: float
    cp_J_mol_K: float


GAS_DATABASE = {
    "H2": Gas("H2", 2.016, 28.84),
    "He": Gas("He", 4.003, 20.79),
    "CH4": Gas("CH4", 16.04, 35.69),
    "NH3": Gas("NH3", 17.03, 35.06),
    "H2O": Gas("H2O", 18.015, 33.58),
    "N2": Gas("N2", 28.014, 29.12),
    "CO": Gas("CO", 28.01, 29.14),
    "O2": Gas("O2", 31.998, 29.38),
    "H2S": Gas("H2S", 34.08, 34.60),
    "CO2": Gas("CO2", 44.01, 37.14),
    "SO2": Gas("SO2", 64.07, 40.45),
}


@dataclass(frozen=True)
class AtmosphereMixture:
    name: str
    candidate_atmosphere: str
    mole_fractions: dict[str, float]
    surface_pressure_bar: float
    albedo: float
    emissivity: float
    notes: str


ATMOSPHERE_MIXTURES = [
    AtmosphereMixture(
        name="N2_CO2_H2O_baseline",
        candidate_atmosphere="N2_CO2_H2O",
        mole_fractions={
            "N2": 0.90,
            "CO2": 0.09,
            "H2O": 0.01,
        },
        surface_pressure_bar=1.0,
        albedo=0.30,
        emissivity=0.85,
        notes="Atmosfera secundária de alto peso molecular, baseline rochoso/water world.",
    ),

    AtmosphereMixture(
        name="CO2_dominated",
        candidate_atmosphere="CO2_N2_H2O",
        mole_fractions={
            "CO2": 0.90,
            "N2": 0.09,
            "H2O": 0.01,
        },
        surface_pressure_bar=2.0,
        albedo=0.35,
        emissivity=0.75,
        notes="Atmosfera CO2-dominante por outgassing, efeito estufa mais forte.",
    ),

    AtmosphereMixture(
        name="reducing_CH4",
        candidate_atmosphere="N2_CO2_CH4_H2O",
        mole_fractions={
            "N2": 0.84,
            "CO2": 0.10,
            "CH4": 0.05,
            "H2O": 0.01,
        },
        surface_pressure_bar=1.0,
        albedo=0.25,
        emissivity=0.80,
        notes="Atmosfera redutora tipo Terra arqueana simplificada.",
    ),

    AtmosphereMixture(
        name="sulfur_volcanic",
        candidate_atmosphere="sulfur_volcanic",
        mole_fractions={
            "CO2": 0.80,
            "N2": 0.15,
            "H2O": 0.03,
            "SO2": 0.01,
            "H2S": 0.01,
        },
        surface_pressure_bar=1.0,
        albedo=0.40,
        emissivity=0.80,
        notes="Atmosfera secundária vulcânica/sulfurada.",
    ),
]

def mean_molar_mass_g_mol(mole_fractions: dict[str, float]) -> float:
    return sum(
        x * GAS_DATABASE[gas].molar_mass_g_mol
        for gas, x in mole_fractions.items()
    )


def cp_mixture_J_mol_K(mole_fractions: dict[str, float]) -> float:
    return sum(
        x * GAS_DATABASE[gas].cp_J_mol_K
        for gas, x in mole_fractions.items()
    )


def column_density_mol_m2(
    surface_pressure_bar: float,
    gravity_m_s2: float,
    mean_molar_mass_g_mol_value: float,
) -> float:
    """
    Coluna molar aproximada:
    massa_coluna = P/g
    mol_coluna = massa_coluna / M
    """
    pressure_Pa = surface_pressure_bar * 1e5
    mean_molar_mass_kg_mol = mean_molar_mass_g_mol_value / 1000.0

    mass_column_kg_m2 = pressure_Pa / gravity_m_s2
    return mass_column_kg_m2 / mean_molar_mass_kg_mol


MIXTURES_BY_CANDIDATE = {}

for mixture in ATMOSPHERE_MIXTURES:
    MIXTURES_BY_CANDIDATE.setdefault(mixture.candidate_atmosphere, []).append(mixture)


def main() -> None:
    input_csv = "atmosphere_candidate_screening.csv"
    planet_csv = "planet_modelling_results.csv"
    output_csv = "atmosphere_physical_parameters.csv"

    candidates_df = pd.read_csv(input_csv)
    planet_df = pd.read_csv(planet_csv)

    planet_cols = [
        "anchor",
        "radiative_case",
        "planet",
        "phi_flux_W_m2",
        "planet_gravity_m_s2",
        "planet_mass_kg",
        "planet_radius_m",
        "planet_density_kg_m3",
    ]

    planet_unique = planet_df[planet_cols].drop_duplicates()

    df = candidates_df.merge(
        planet_unique,
        on=["anchor", "radiative_case", "planet"],
        how="left",
    )

    useful_priorities = {"alta", "média", "média-baixa"}

    df = df[df["priority"].isin(useful_priorities)].copy()

    rows = []

    for _, row in df.iterrows():
        candidate_name = row["candidate_atmosphere"]

        mixtures = MIXTURES_BY_CANDIDATE.get(candidate_name, [])

        if not mixtures:
            continue

        for mixture in mixtures:
            mu = mean_molar_mass_g_mol(mixture.mole_fractions)
            cp_mix = cp_mixture_J_mol_K(mixture.mole_fractions)

            column_mol_m2 = column_density_mol_m2(
                surface_pressure_bar=mixture.surface_pressure_bar,
                gravity_m_s2=row["planet_gravity_m_s2"],
                mean_molar_mass_g_mol_value=mu,
            )

            rows.append({
                "anchor": row["anchor"],
                "radiative_case": row["radiative_case"],
                "planet": row["planet"],
                "T_K": row["T_K"],
                "T_C": row["T_C"],
                "escape_velocity_km_s": row["escape_velocity_km_s"],
                "phi_flux_W_m2": row["phi_flux_W_m2"],
                "planet_gravity_m_s2": row["planet_gravity_m_s2"],

                "candidate_atmosphere": candidate_name,
                "candidate_priority": row["priority"],
                "mixture_name": mixture.name,
                "mole_fractions": mixture.mole_fractions,

                "surface_pressure_bar": mixture.surface_pressure_bar,
                "mean_molar_mass_g_mol": mu,
                "cp_mix_J_mol_K": cp_mix,
                "column_density_mol_m2": column_mol_m2,

                "albedo_assumed": mixture.albedo,
                "emissivity_assumed": mixture.emissivity,
                "mixture_notes": mixture.notes,
            })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de linhas: {len(out)}")

    if not out.empty:
        print(
            out[
                [
                    "anchor",
                    "radiative_case",
                    "planet",
                    "candidate_atmosphere",
                    "mixture_name",
                    "mean_molar_mass_g_mol",
                    "cp_mix_J_mol_K",
                    "column_density_mol_m2",
                    "albedo_assumed",
                    "emissivity_assumed",
                ]
            ].to_string(index=False)
        )



if __name__ == "__main__":
    main()