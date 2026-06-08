from __future__ import annotations
from atmospheres_candidates import ATMOSPHERE_CANDIDATES


def classify_thermal_plausibility(T_K: float) -> str:
    if T_K < 240.0:
        return "baixa: frio demais sem efeito estufa forte"
    if T_K < 273.15:
        return "média: frio, mas efeito estufa pode compensar"
    if T_K <= 310.0:
        return "alta: faixa temperada"
    if T_K <= 373.15:
        return "média: quente, avaliar perda de água"
    return "baixa: risco de superaquecimento"


def suggest_atmospheres(anchor_name: str, radiative_case: str, T_K: float) -> list[dict]:
    suggestions = []

    for atm in ATMOSPHERE_CANDIDATES:
        thermal_match = atm.min_T_K <= T_K <= atm.max_T_K

        suggestions.append({
            "anchor": anchor_name,
            "radiative_case": radiative_case,
            "T_K": round(T_K, 2),
            "T_C": round(T_K - 273.15, 2),
            "atmosphere": atm.name,
            "composition": atm.composition,
            "category": atm.category,
            "base_priority": atm.priority,
            "thermal_match": thermal_match,
            "thermal_classification": classify_thermal_plausibility(T_K),
            "notes": atm.notes,
            "possible_biosignatures": ", ".join(atm.possible_biosignatures),
            "main_false_positives": ", ".join(atm.main_false_positives),
        })

    return suggestions