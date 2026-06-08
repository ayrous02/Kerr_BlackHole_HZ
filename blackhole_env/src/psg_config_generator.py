from __future__ import annotations

from pathlib import Path
import re
import math
import pandas as pd


TEMPLATE_PATH = Path("psg_cfg.txt")
READY_CASES_CSV = Path("psg_ready_cases.csv")
OUTPUT_DIR = Path("psg_configs")

SOLAR_CONSTANT_W_M2 = 1366.0


PSG_GAS_TYPE = {
    # Mapeamento inicial baseado em HITRAN/ExoMol usados pelo PSG.
    # Se o PSG reclamar de algum gás, ajuste pela interface e baixe um novo template.
    "H2": "HIT[45]",
    "He": "HIT[0]",
    "H2O": "EXO[1]",
    "CH4": "EXO[6]",
    "NH3": "HIT[11]",
    "N2": "HIT[22]",
    "CO": "HIT[5]",
    "O2": "HIT[7]",
    "CO2": "HIT[2]",
    "SO2": "HIT[9]",
    "H2S": "HIT[31]",
}


def parse_molecules(molecules_string: str) -> dict[str, float]:
    """
    Converte:
    'N2=0.84; CO2=0.1; CH4=0.05; H2O=0.01'
    em:
    {'N2': 0.84, 'CO2': 0.1, 'CH4': 0.05, 'H2O': 0.01}
    """
    molecules = {}

    if not isinstance(molecules_string, str):
        return molecules

    for part in molecules_string.split(";"):
        part = part.strip()
        if not part:
            continue

        gas, value = part.split("=")
        molecules[gas.strip()] = float(value.strip())

    return molecules


def replace_tag(config: str, tag: str, value: str) -> str:
    """
    Substitui uma tag do tipo:
    <TAG>valor
    por:
    <TAG>novo_valor

    Se a tag não existir, adiciona no final.
    """
    pattern = rf"(<{re.escape(tag)}>)[^\n\r]*"

    if re.search(pattern, config):
        return re.sub(
            pattern,
            lambda match: f"{match.group(1)}{value}",
            config,
            count=1,
        )

    return config + f"\n<{tag}>{value}"


def remove_atmosphere_layers(config: str) -> str:
    """
    Remove o perfil vertical detalhado do template GJ 1214b.

    Isso evita que o PSG continue usando as camadas antigas do template
    em vez da mistura simplificada gerada pelo nosso pipeline.
    """
    lines = config.splitlines()
    cleaned_lines = []

    for line in lines:
        if line.startswith("<ATMOSPHERE-LAYER-"):
            continue
        if line.startswith("<ATMOSPHERE-LAYERS>"):
            continue
        if line.startswith("<ATMOSPHERE-LAYERS-MOLECULES>"):
            continue

        cleaned_lines.append(line)

    return "\n".join(cleaned_lines) + "\n"


def make_safe_case_id(case_id: str) -> str:
    return (
        str(case_id)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("–", "_")
        .replace("-", "_")
        .replace("__", "_")
    )


def equivalent_solar_distance_AU(phi_flux_W_m2: float) -> float:
    """
    Distância equivalente a uma estrela solar para produzir o mesmo fluxo bolométrico.

    Isso NÃO é a distância física do planeta ao buraco negro.
    É apenas um placeholder útil para o PSG, que espera uma fonte estelar convencional.
    """
    return math.sqrt(SOLAR_CONSTANT_W_M2 / phi_flux_W_m2)


def build_config_for_case(template: str, row: pd.Series) -> str:
    config = template

    case_id = make_safe_case_id(row["case_id"])

    molecules = parse_molecules(row["molecules_for_psg"])
    gases = list(molecules.keys())
    abundances = [molecules[gas] for gas in gases]

    gas_types = [
        PSG_GAS_TYPE.get(gas, "HIT[0]")
        for gas in gases
    ]

    diameter_km = float(row["planet_diameter_km"])
    gravity_m_s2 = float(row["planet_gravity_m_s2"])
    temperature_K = float(row["T_surface_K_for_PSG"])
    pressure_bar = float(row["surface_pressure_bar"])
    mean_molar_mass = float(row["mean_molar_mass_g_mol"])

    phi_flux = float(row["phi_flux_W_m2"])
    pseudo_star_distance_AU = equivalent_solar_distance_AU(phi_flux)

    config = remove_atmosphere_layers(config)

    # =========================
    # Objeto
    # =========================
    config = replace_tag(config, "OBJECT", "Exoplanet")
    config = replace_tag(config, "OBJECT-NAME", case_id)
    config = replace_tag(config, "OBJECT-DIAMETER", f"{diameter_km:.6f}")
    config = replace_tag(config, "OBJECT-GRAVITY", f"{gravity_m_s2:.6f}")
    config = replace_tag(config, "OBJECT-GRAVITY-UNIT", "g")

    # Placeholder estelar para o PSG.
    # O cenário real é CMB blueshiftada, então isso serve só para piloto espectral.
    config = replace_tag(config, "OBJECT-STAR-TYPE", "G")
    config = replace_tag(config, "OBJECT-STAR-TEMPERATURE", "5777")
    config = replace_tag(config, "OBJECT-STAR-RADIUS", "1.0")
    config = replace_tag(config, "OBJECT-STAR-DISTANCE", f"{pseudo_star_distance_AU:.6f}")

    # =========================
    # Atmosfera
    # =========================
    config = replace_tag(config, "ATMOSPHERE-DESCRIPTION", f"BlackSun candidate: {case_id}")
    config = replace_tag(config, "ATMOSPHERE-STRUCTURE", "Equilibrium")
    config = replace_tag(config, "ATMOSPHERE-WEIGHT", f"{mean_molar_mass:.6f}")
    config = replace_tag(config, "ATMOSPHERE-PRESSURE", f"{pressure_bar:.6f}")
    config = replace_tag(config, "ATMOSPHERE-PUNIT", "bar")
    config = replace_tag(config, "ATMOSPHERE-TEMPERATURE", f"{temperature_K:.6f}")
    config = replace_tag(config, "ATMOSPHERE-NGAS", str(len(gases)))
    config = replace_tag(config, "ATMOSPHERE-GAS", ",".join(gases))
    config = replace_tag(config, "ATMOSPHERE-TYPE", ",".join(gas_types))
    config = replace_tag(config, "ATMOSPHERE-ABUN", ",".join(f"{x:.8g}" for x in abundances))
    config = replace_tag(config, "ATMOSPHERE-UNIT", ",".join(["ratio"] * len(gases)))
    config = replace_tag(config, "ATMOSPHERE-TAU", ",".join(["1"] * len(gases)))

    # =========================
    # Superfície
    # =========================
    config = replace_tag(config, "SURFACE-TEMPERATURE", f"{temperature_K:.6f}")
    config = replace_tag(config, "SURFACE-ALBEDO", f"{float(row['albedo_assumed']):.6f}")
    config = replace_tag(config, "SURFACE-EMISSIVITY", f"{float(row['emissivity_assumed']):.6f}")

    # =========================
    # Gerador espectral
    # =========================
    config = replace_tag(config, "GENERATOR-RANGE1", "0.5")
    config = replace_tag(config, "GENERATOR-RANGE2", "17")
    config = replace_tag(config, "GENERATOR-RANGEUNIT", "um")
    config = replace_tag(config, "GENERATOR-RESOLUTION", "200")
    config = replace_tag(config, "GENERATOR-RESOLUTIONUNIT", "RP")

    # =========================
    # Comentários/metadados
    # =========================
    config += "\n# BLACK SUN PIPELINE NOTES\n"
    config += f"# Anchor: {row['anchor']}\n"
    config += f"# Planet: {row['planet']}\n"
    config += f"# Candidate atmosphere: {row['candidate_atmosphere']}\n"
    config += f"# Mixture: {row['mixture_name']}\n"
    config += f"# Climate priority: {row['climate_priority']}\n"
    config += f"# Phi flux [W/m2]: {row['phi_flux_W_m2']}\n"
    config += f"# Spin: {row['spin']}\n"
    config += f"# r_orb [GM/c2]: {row['r_orb_GM_c2']}\n"
    config += f"# orbit_radius_m: {row['orbit_radius_m']}\n"
    config += f"# orbit_radius_AU: {row['orbit_radius_AU']}\n"
    config += f"# pseudo_star_distance_AU_for_same_flux: {pseudo_star_distance_AU:.6f}\n"
    config += f"# time_dilation_gamma: {row['time_dilation_gamma']}\n"
    config += f"# roche_margin: {row['roche_margin']}\n"
    config += f"# survives_tidal: {row['survives_tidal']}\n"
    config += f"# Molecules: {row['molecules_for_psg']}\n"
    config += f"# Retained gases: {row.get('retained_gases', '')}\n"
    config += f"# Scientific reason: {row.get('scientific_reason', '')}\n"
    config += f"# Notes: {row.get('notes', '')}\n"
    config += "# Limitation: PSG template assumes conventional illumination geometry.\n"
    config += "# The Bakala source is blueshifted CMB, so this config is a spectral-feature pilot.\n"

    return config


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    template = TEMPLATE_PATH.read_text(encoding="utf-8")
    df = pd.read_csv(READY_CASES_CSV)

    # Para o primeiro teste, mantenha só os casos Earth-like.
    # Depois remova esta linha para gerar todos.
    df = df[df["planet"] == "Earth-like"].copy()

    generated_files = []

    for _, row in df.iterrows():
        case_id = make_safe_case_id(row["case_id"])
        config = build_config_for_case(template, row)

        out_path = OUTPUT_DIR / f"{case_id}.txt"
        out_path.write_text(config, encoding="utf-8")

        generated_files.append(out_path)

    print(f"Configs gerados: {len(generated_files)}")
    print(f"Pasta: {OUTPUT_DIR.resolve()}")

    for path in generated_files:
        print(path)


if __name__ == "__main__":
    main()