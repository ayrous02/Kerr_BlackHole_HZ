# leitura dos arquivos gerados anteriormente para agrupar em um arquivo de fongração para dar upload no PSG

from __future__ import annotations
import ast
from pathlib import Path
import pandas as pd


def parse_mole_fractions(value: str) -> dict[str, float]:
    """
    Converte string de dicionário vinda do CSV em dict Python.

    Exemplo:
    "{'N2': 0.84, 'CO2': 0.10, 'CH4': 0.05, 'H2O': 0.01}"
    """
    if isinstance(value, dict):
        return value

    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return {}

    if not isinstance(parsed, dict):
        return {}

    return parsed


def format_molecules_for_psg_table(mole_fractions: dict[str, float]) -> str:
    """
    Formata as abundâncias para leitura humana.
    Por enquanto, isso NÃO é ainda o formato final do config PSG.
    """
    return "; ".join(
        f"{gas}={fraction:.6g}"
        for gas, fraction in mole_fractions.items()
    )


def select_cases(df: pd.DataFrame) -> pd.DataFrame:
    """
    Seleciona cenários prioritários para teste piloto no PSG.

    Estratégia:
    - prioridade climática alta ou média
    - preferir cenários Earth-like
    - manter alguns controles científicos
    """

    selected = df[
        df["climate_priority"].isin(["alta", "média"])
    ].copy()

    # Mantém só cenários únicos: nesta etapa radiative_case antigo pode gerar duplicatas.
    selected = selected.drop_duplicates(
        subset=[
            "anchor",
            "planet",
            "candidate_atmosphere",
            "mixture_name",
            "T_eq_atm_K",
        ]
    )

    # Ordena os melhores primeiro
    priority_order = {
        "alta": 0,
        "média": 1,
        "média-baixa": 2,
        "baixa": 3,
    }

    selected["priority_order"] = selected["climate_priority"].map(priority_order)

    selected = selected.sort_values(
        by=[
            "priority_order",
            "anchor",
            "planet",
            "T_eq_atm_K",
        ]
    )

    return selected


def make_psg_notes(row: pd.Series) -> str:
    """
    Cria uma nota curta explicando por que esse cenário foi selecionado.
    """
    atm = row["candidate_atmosphere"]
    anchor = row["anchor"]
    T = row["T_eq_atm_K"]

    if atm == "N2_CO2_CH4_H2O":
        return (
            f"Cenário redutor temperado na âncora {anchor}; "
            f"bom para testar CH4, CO2 e H2O. "
            f"T_eq_atm={T:.2f} K. Atenção: CH4 pode ser vulnerável à fotólise UV."
        )

    if atm == "sulfur_volcanic":
        return (
            f"Cenário vulcânico/sulfurado na âncora {anchor}; "
            f"bom para testar SO2, H2S, CO2 e aerossóis. "
            f"T_eq_atm={T:.2f} K."
        )

    if atm == "N2_CO2_H2O":
        return (
            f"Baseline de alto peso molecular na âncora {anchor}; "
            f"bom para testar CO2 e H2O. "
            f"T_eq_atm={T:.2f} K."
        )

    if atm == "CO2_N2_H2O":
        return (
            f"Cenário CO2-dominante na âncora {anchor}; "
            f"bom para testar efeito estufa e bandas fortes de CO2. "
            f"T_eq_atm={T:.2f} K."
        )

    return f"Cenário selecionado para teste espectral. T_eq_atm={T:.2f} K."


def main() -> None:
    input_csv = "ranked_climate_scenarios.csv"
    output_csv = "psg_selected_cases.csv"
    output_dir = Path("psg_case_cards")
    output_dir.mkdir(exist_ok=True)

    df = pd.read_csv(input_csv)

    selected = select_cases(df)

    rows = []

    for idx, row in selected.iterrows():
        mole_fractions = parse_mole_fractions(row["mole_fractions"])

        case_id = (
            f"{row['anchor']}_"
            f"{row['planet']}_"
            f"{row['candidate_atmosphere']}_"
            f"{row['mixture_name']}"
        )
        case_id = (
            case_id.replace(" ", "_")
            .replace("/", "_")
            .replace("–", "_")
            .replace("-", "_")
        )

        psg_notes = make_psg_notes(row)

        rows.append({
            "case_id": case_id,
            "anchor": row["anchor"],
            "planet": row["planet"],
            "candidate_atmosphere": row["candidate_atmosphere"],
            "mixture_name": row["mixture_name"],
            "climate_priority": row["climate_priority"],
            "T_surface_K_for_PSG": row["T_eq_atm_K"],
            "T_surface_C_for_PSG": row["T_eq_atm_C"],
            "surface_pressure_bar": row["surface_pressure_bar"],
            "surface_pressure_Pa": row["surface_pressure_bar"] * 1e5,
            "planet_gravity_m_s2": row["planet_gravity_m_s2"],
            "planet_radius_m": row.get("planet_radius_m", None),
            "mean_molar_mass_g_mol": row["mean_molar_mass_g_mol"],
            "cp_mix_J_mol_K": row["cp_mix_J_mol_K"],
            "albedo_assumed": row["albedo_assumed"],
            "emissivity_assumed": row["emissivity_assumed"],
            "molecules_for_psg": format_molecules_for_psg_table(mole_fractions),
            "notes": psg_notes,
        })

        # Cria um cartão de caso em .txt para você preencher manualmente no PSG
        case_card = output_dir / f"{case_id}.txt"

        with case_card.open("w", encoding="utf-8") as f:
            f.write(f"CASE ID: {case_id}\n")
            f.write(f"Anchor Bakala: {row['anchor']}\n")
            f.write(f"Planet: {row['planet']}\n")
            f.write(f"Atmosphere candidate: {row['candidate_atmosphere']}\n")
            f.write(f"Mixture: {row['mixture_name']}\n")
            f.write(f"Climate priority: {row['climate_priority']}\n\n")

            f.write("=== Values to enter in PSG ===\n")
            f.write(f"Surface / atmospheric temperature [K]: {row['T_eq_atm_K']:.3f}\n")
            f.write(f"Surface pressure [bar]: {row['surface_pressure_bar']:.3f}\n")
            f.write(f"Gravity [m/s2]: {row['planet_gravity_m_s2']:.5f}\n")
            f.write(f"Mean molar mass [g/mol]: {row['mean_molar_mass_g_mol']:.5f}\n")
            f.write(f"Albedo assumed: {row['albedo_assumed']:.5f}\n")
            f.write(f"Emissivity assumed: {row['emissivity_assumed']:.5f}\n")
            f.write("\nMolecular abundances:\n")

            for gas, fraction in mole_fractions.items():
                f.write(f"  {gas}: {fraction:.8f}\n")

            f.write("\nScientific note:\n")
            f.write(psg_notes + "\n")

            f.write("\nImportant limitation:\n")
            f.write(
                "This is a simplified atmosphere derived from the Black Sun pipeline. "
                "The incident source in Bakala is blueshifted CMB, not a standard stellar spectrum. "
                "Use PSG first to inspect spectral features, not as final proof of atmospheric stability.\n"
            )

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Cartões salvos na pasta: {output_dir.resolve()}")
    print(f"Total de casos selecionados: {len(out)}")

    print(
        out[
            [
                "case_id",
                "climate_priority",
                "T_surface_K_for_PSG",
                "surface_pressure_bar",
                "molecules_for_psg",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()