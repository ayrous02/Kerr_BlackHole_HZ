from __future__ import annotations
import pandas as pd




'''
lista dos cenários mais interessantes para seguir para a etapa seguinte
'''



def rank_climate_scenario(T_K: float) -> tuple[str, str]:
    if T_K < 240:
        return "baixa", "muito frio; exigiria aquecimento/efeito estufa muito forte"
    if T_K < 273.15:
        return "média-baixa", "frio/marginal; pode ser salvo por maior pressão ou efeito estufa"
    if T_K <= 310:
        return "alta", "faixa temperada, compatível com água líquida em princípio"
    if T_K <= 340:
        return "média", "quente, mas ainda interessante; avaliar perda de água e UV"
    if T_K <= 373.15:
        return "média-baixa", "muito quente; risco de estresse térmico e perda de água"
    return "baixa", "quente demais no modelo 0D"


def main() -> None:
    input_csv = "atmospheric_climate_results.csv"
    output_csv = "ranked_climate_scenarios.csv"

    df = pd.read_csv(input_csv)

    df = df.drop_duplicates(
        subset=[
            "anchor",
            "planet",
            "candidate_atmosphere",
            "mixture_name",
            "T_eq_atm_K",
        ]
    ).copy()

    priorities = []
    comments = []

    for _, row in df.iterrows():
        priority, comment = rank_climate_scenario(row["T_eq_atm_K"])
        priorities.append(priority)
        comments.append(comment)

    df["climate_priority"] = priorities
    df["ranking_comment"] = comments

    priority_order = {
        "alta": 0,
        "média": 1,
        "média-baixa": 2,
        "baixa": 3,
    }

    df["priority_order"] = df["climate_priority"].map(priority_order)

    df = df.sort_values(
        by=["priority_order", "anchor", "planet", "T_eq_atm_K"]
    )

    df.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de cenários únicos: {len(df)}")

    print(
        df[
            [
                "anchor",
                "planet",
                "candidate_atmosphere",
                "mixture_name",
                "T_eq_atm_K",
                "T_eq_atm_C",
                "climate_priority",
                "ranking_comment",
                "radiative_timescale_days",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()