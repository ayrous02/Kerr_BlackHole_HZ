from __future__ import annotations
from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class AtmosphereCandidate:
    name: str
    required_gases: tuple[str, ...]
    optional_gases: tuple[str, ...]
    category: str
    scientific_reason: str
'''
    composition: str
    category: str
    priority: str
    min_reasonable_T_K: float
    max_reasonable_T_K: float
    notes: str
    possible_biosignatures: tuple[str, ...]
    main_false_positives: tuple[str, ...]
'''



ATMOSPHERE_CANDIDATES = [
    AtmosphereCandidate(
        name="N2_CO2_H2O",
        required_gases=("N2", "CO2", "H2O"),
        optional_gases=("CH4", "O2", "O3"),
        category="atmosfera secundária rochosa / alto peso molecular",
        scientific_reason=(
            "Cenário baseline para planeta rochoso ou water world; inspirado em atmosferas "
            "N2/CO2 de alto peso molecular discutidas por Damiano et al."
        ),
    ),
    AtmosphereCandidate(
        name="CO2_N2_H2O",
        required_gases=("CO2", "N2", "H2O"),
        optional_gases=("SO2", "H2S", "CO"),
        category="CO2-dominante / outgassing",
        scientific_reason=(
            "Cenário plausível para atmosfera secundária por outgassing, com efeito estufa forte."
        ),
    ),
    AtmosphereCandidate(
        name="N2_CO2_CH4_H2O",
        required_gases=("N2", "CO2", "CH4", "H2O"),
        optional_gases=("H2", "CO"),
        category="redutora / tipo Terra arqueana",
        scientific_reason=(
            "Interessante para metanogênese, desequilíbrio CH4+CO2 e formação de haze; "
            "precisará passar depois por filtro UV/fotólise."
        ),
    ),
    AtmosphereCandidate(
        name="H2_rich",
        required_gases=("H2",),
        optional_gases=("He", "CH4", "NH3", "H2O", "CO2"),
        category="primordial / Hycean / super-Terra",
        scientific_reason=(
            "Só deve ser prioridade se H2/He forem bem retidos e o planeta for massivo o suficiente."
        ),
    ),
    AtmosphereCandidate(
        name="sulfur_volcanic",
        required_gases=("CO2", "N2"),
        optional_gases=("SO2", "H2S", "H2O"),
        category="secundária vulcânica / sulfurada",
        scientific_reason=(
            "Cenário inspirado em atmosferas secundárias com SO2/H2S; útil para testar vulcanismo "
            "e falsos positivos."
        ),
    ),
    AtmosphereCandidate(
        name="airless_or_thin",
        required_gases=(),
        optional_gases=(),
        category="controle negativo",
        scientific_reason=(
            "Controle para cenários em que a atmosfera é fina ou ausente."
        ),
    ),
]

'''
ATMOSPHERE_CANDIDATES = [
    AtmosphereCandidate(
        name="N2_CO2_H2O",
        composition="N2-dominante com CO2 e H2O",
        category="atmosfera secundária rochosa / alto peso molecular",
        priority="alta",
        min_reasonable_T_K=250.0,
        max_reasonable_T_K=330.0,
        notes=(
            "Boa atmosfera baseline para planeta rochoso ou water world. "
            "Serve para testar efeito estufa moderado, água líquida e proteção parcial contra UV."
        ),
        possible_biosignatures=("O2", "O3", "CH4", "N2O"),
        main_false_positives=("O2 abiótico por fotólise de H2O/CO2", "O3 abiótico", "CO acumulado"),
    ),

    AtmosphereCandidate(
        name="CO2_N2_H2O",
        composition="CO2-dominante com N2 e H2O",
        category="outgassing / atmosfera secundária CO2-rich",
        priority="alta",
        min_reasonable_T_K=220.0,
        max_reasonable_T_K=350.0,
        notes=(
            "Importante para borda fria e para mundos com vulcanismo/outgassing. "
            "Pode aquecer planetas frios, mas sob UV forte pode gerar CO, O e O2/O3 abióticos."
        ),
        possible_biosignatures=("CH4 + CO2", "O2/O3 com cautela"),
        main_false_positives=("CO2 photolysis", "O2/O3 abiótico", "CO alto"),
    ),

    AtmosphereCandidate(
        name="N2_CO2_CH4_H2O",
        composition="N2, CO2, CH4 e H2O",
        category="redutora / análoga à Terra arqueana",
        priority="média-alta",
        min_reasonable_T_K=250.0,
        max_reasonable_T_K=330.0,
        notes=(
            "Boa para testar metanogênese, desequilíbrio CH4+CO2 e formação de haze orgânico. "
            "No cenário Bakala, o UV pode destruir CH4 rapidamente."
        ),
        possible_biosignatures=("CH4 + CO2", "haze orgânico", "variação sazonal de CH4"),
        main_false_positives=("CH4 abiótico por serpentinização", "haze abiótico"),
    ),

    AtmosphereCandidate(
        name="CO2_N2_SO2_H2S_H2O",
        composition="CO2/N2 com SO2, H2S e H2O",
        category="vulcânica / sulfurada",
        priority="média",
        min_reasonable_T_K=220.0,
        max_reasonable_T_K=360.0,
        notes=(
            "Útil para cenários de vulcanismo ativo e química de enxofre. "
            "Pode produzir aerossóis e alterar fortemente o albedo."
        ),
        possible_biosignatures=("H2S", "DMS", "DMDS", "desequilíbrio sulfurado"),
        main_false_positives=("SO2 vulcânico", "H2S vulcânico", "aerossóis abióticos"),
    ),

    AtmosphereCandidate(
        name="H2_CO2_CH4_H2O",
        composition="H2-dominante com CO2, CH4 e H2O",
        category="H2-rich / Hycean / super-Terra",
        priority="média-baixa",
        min_reasonable_T_K=230.0,
        max_reasonable_T_K=320.0,
        notes=(
            "Interessante para super-Terras ou mini-Netunos, mas menos defensável para planeta tipo Terra. "
            "Precisa passar por filtro de retenção atmosférica e risco de runaway greenhouse."
        ),
        possible_biosignatures=("CH4", "DMS", "DMDS", "NH3 com cautela"),
        main_false_positives=("CH4 abiótico", "química atmosférica em H2", "degenerescência espectral"),
    ),

    AtmosphereCandidate(
        name="airless_or_thin",
        composition="sem atmosfera ou atmosfera muito fina",
        category="controle negativo",
        priority="controle",
        min_reasonable_T_K=0.0,
        max_reasonable_T_K=1000.0,
        notes=(
            "Serve como comparação. Não sustenta clima estável nem proteção UV significativa, "
            "mas ajuda a mostrar o papel da atmosfera."
        ),
        possible_biosignatures=(),
        main_false_positives=("assinaturas minerais/superfície confundidas com atmosfera",),
    ),
]
'''


GOOD_RETENTION = {
    "retenção plausível",
    "retenção moderada",
}


def get_retained_gases(group: pd.DataFrame) -> set[str]:
    retained = group[group["retention_class"].isin(GOOD_RETENTION)]
    return set(retained["gas"].tolist())


def classify_atmosphere_candidate(
    candidate: AtmosphereCandidate,
    retained_gases: set[str],
    T_K: float,
    planet: str,
) -> tuple[str, str]:
    """
    Retorna prioridade e comentário.
    """

    if candidate.name == "airless_or_thin":
        return "controle", "Cenário controle, não depende de retenção de gases específicos."

    missing_required = [
        gas for gas in candidate.required_gases
        if gas not in retained_gases
    ]

    if missing_required:
        return (
            "baixa",
            f"Faltam gases essenciais no filtro de retenção: {', '.join(missing_required)}."
        )

    # Regras térmicas simples
    if T_K < 240:
        thermal_note = "cenário frio; exigiria efeito estufa forte."
    elif T_K < 273.15:
        thermal_note = "cenário frio/marginal; efeito estufa pode compensar."
    elif T_K <= 310:
        thermal_note = "faixa temperada."
    elif T_K <= 373.15:
        thermal_note = "cenário quente; avaliar perda de água e efeito estufa."
    else:
        thermal_note = "cenário muito quente; risco de superaquecimento."

    # Regra especial para H2-rich
    if candidate.name == "H2_rich":
        if planet in {"Earth-like", "Mars-like"}:
            return (
                "média-baixa",
                "H2 aparece retível no filtro simples, mas este filtro usa T_eq e pode superestimar retenção de H2; tratar como cenário exploratório. "
                + thermal_note
            )
        return (
            "média",
            "H2-rich pode ser testado, mas precisa de etapa posterior com escape hidrodinâmico/fotólise. "
            + thermal_note
        )

    # Regra especial para sulfurada
    if candidate.name == "sulfur_volcanic":
        sulfur_gases = {"SO2", "H2S"} & retained_gases
        if sulfur_gases:
            return (
                "média",
                f"Gases sulfurados retidos preliminarmente ({', '.join(sorted(sulfur_gases))}); depois testar fotólise e vulcanismo. "
                + thermal_note
            )
        return (
            "média-baixa",
            "Base CO2/N2 retida, mas SO2/H2S não aparecem como robustos; manter como cenário secundário. "
            + thermal_note
        )

    return (
        "alta",
        "Gases essenciais passam no filtro de retenção. " + thermal_note
    )


def main() -> None:
    input_csv = "gas_retention_screening.csv"
    output_csv = "atmosphere_candidate_screening.csv"

    df = pd.read_csv(input_csv)

    scenario_cols = [
        "anchor",
        "radiative_case",
        "planet",
        "T_K",
        "T_C",
        "escape_velocity_km_s",
    ]

    rows = []

    for scenario_values, group in df.groupby(scenario_cols):
        scenario = dict(zip(scenario_cols, scenario_values))
        retained_gases = get_retained_gases(group)

        for candidate in ATMOSPHERE_CANDIDATES:
            priority, comment = classify_atmosphere_candidate(
                candidate=candidate,
                retained_gases=retained_gases,
                T_K=scenario["T_K"],
                planet=scenario["planet"],
            )

            rows.append({
                **scenario,
                "candidate_atmosphere": candidate.name,
                "category": candidate.category,
                "priority": priority,
                "required_gases": ", ".join(candidate.required_gases) or "-",
                "optional_gases": ", ".join(candidate.optional_gases) or "-",
                "retained_gases": ", ".join(sorted(retained_gases)),
                "scientific_reason": candidate.scientific_reason,
                "comment": comment,
            })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de linhas: {len(out)}")

    # Tabela-resumo rápida no terminal
    summary = out[out["priority"].isin(["alta", "média", "média-baixa"])]
    print("\nAtmosferas candidatas priorizadas:")
    print(
        summary[
            ["anchor", "radiative_case", "planet", "candidate_atmosphere", "priority", "comment"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()