#1) Tabela. com a velocidade de escape do planeta e a temperatura de equilíbrio, quais gases ele poderia reter?
from __future__ import annotations
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from constants import AMU_TO_KG, K_B


GASES = [
    {"name": "H2", "molecular_mass_amu": 2.016},
    {"name": "He", "molecular_mass_amu": 4.003},
    {"name": "CH4", "molecular_mass_amu": 16.04},
    {"name": "NH3", "molecular_mass_amu": 17.03},
    {"name": "H2O", "molecular_mass_amu": 18.015},
    {"name": "N2", "molecular_mass_amu": 28.014},
    {"name": "CO", "molecular_mass_amu": 28.01},
    {"name": "O2", "molecular_mass_amu": 31.998},
    {"name": "H2S", "molecular_mass_amu": 34.08},
    {"name": "CO2", "molecular_mass_amu": 44.01},
    {"name": "SO2", "molecular_mass_amu": 64.07},
#    {"name": "C", "molecular_mass_amu": 12.011},
 #   {"name": "O", "molecular_mass_amu": 15.99},
]

#calcular o escape térmico de cada um dos GASES. 
''' parâmetro de Jeans para cada gás
λ=( (m*v_esc**2) / (2*k_b*T) )
'''
def jeans_parameter(
    molecular_mass_amu: float,
    escape_velocity_km_s: float,
    temperature_K: float,
) -> float:
    """
    - v_escape entra em m/s
    - massa molecular entra em kg
    - T entra em K
    """
    molecular_mass_kg = molecular_mass_amu * AMU_TO_KG
    escape_velocity_m_s = escape_velocity_km_s * 1000.0

    return molecular_mass_kg * escape_velocity_m_s**2 / (2.0 * K_B * temperature_K)


'''
Calcula a velocidade de escape de cada gás
Isola da função acima
'''
def escape_velocity_for_lambda(
    molecular_mass_amu: float,
    temperature_K: float,
    lambda_threshold: float = 30.0,
) -> float:

    molecular_mass_kg = molecular_mass_amu * AMU_TO_KG
    escape_velocity_m_s = math.sqrt(
        (2.0 * lambda_threshold * K_B * temperature_K) / molecular_mass_kg
    )
    return escape_velocity_m_s / 1000.0

''' classificar retenção
λ < 10      → gás dificilmente retido
10-20       → retenção frágil / depende de reposição
20-30       → retenção moderada
λ > 30      → retenção mais plausível
'''
def classify_retention(lambda_value: float) -> str:
    """
    Classificação preliminar baseada no parâmetro de Jeans.

    Este filtro avalia apenas a tendência de retenção térmica
    usando velocidade de escape e temperatura média/equilíbrio.

    Ele NÃO substitui um modelo completo de escape atmosférico,
    pois não inclui:
    - temperatura da exosfera;
    - fluxo UV/XUV;
    - fotólise;
    - escape hidrodinâmico;
    - escape não térmico;
    - vento/partículas;
    - reposição por vulcanismo ou outgassing.
    """
    if lambda_value < 10:
        return "improvável reter"
    if lambda_value < 20:
        return "retenção frágil"
    if lambda_value < 30:
        return "retenção moderada"
    return "retenção plausível"


''' 
Ler planet_modelling_results.csv
'''
def main() -> None:
    input_csv = "planet_modelling_results.csv"
    output_csv = "gas_retention_screening.csv"

    df = pd.read_csv(input_csv)

    rows = []

    for _, row in df.iterrows():
        T_K = row["T_K"]
        escape_velocity_km_s = row["planet_escape_velocity_km_s"]

        for gas in GASES:
            lambda_value = jeans_parameter(
                molecular_mass_amu=gas["molecular_mass_amu"],
                escape_velocity_km_s=escape_velocity_km_s,
                temperature_K=T_K,
            )

            rows.append({
                "anchor": row["anchor"],
                "radiative_case": row["radiative_case"],
                "planet": row["planet"],
                "T_K": T_K,
                "T_C": row["T_C"],
                "escape_velocity_km_s": escape_velocity_km_s,
                "gas": gas["name"],
                "molecular_mass_amu": gas["molecular_mass_amu"],
                "jeans_lambda": lambda_value,
                "retention_class": classify_retention(lambda_value),
            })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)

    print(f"Arquivo salvo: {output_csv}")
    print(f"Total de linhas: {len(out)}")

    # Gráfico de posição dos cenários planetários + linhas de retenção dos gases
    plt.figure(figsize=(10, 7))

    from adjustText import adjust_text
    texts = []
#    gases_ordenados = sorted(GASES, key=lambda x: x["molecular_mass_amu"])

    # eixo de temperatura em escala log
    T_min = max(50.0, df["T_K"].min() * 0.5)
    T_max = df["T_K"].max() * 2.0
    temperature_grid = np.logspace(
        math.log10(T_min),
        math.log10(T_max),
        300,
    )

    lambda_threshold = 30.0

    gases_ordenados = sorted(GASES, key=lambda x: x["molecular_mass_amu"])
    fontes_posicoes_x = np.linspace(0.4, 0.85, len(gases_ordenados))
    i = 0

    for gas in GASES:

        print("gas: ", gas)
        v_line = [
            escape_velocity_for_lambda(
                molecular_mass_amu=gas["molecular_mass_amu"],
                temperature_K=T,
                lambda_threshold=lambda_threshold,
            )
            for T in temperature_grid
        ]

        line, = plt.plot(
            temperature_grid,
            v_line,
            linewidth=1.0,
            alpha=0.6,          
            label=f"{gas['name']} λ={lambda_threshold:.0f}",
        )

        idx = int(len(temperature_grid) * fontes_posicoes_x[i])
        x_pos = temperature_grid[idx]
        y_pos = v_line[idx]

        print("pos x: ", x_pos)
        print("pos y: ", y_pos)

        angulo_inclinacao = 19.5 

        plt.text(
            x_pos,
            y_pos * 1,                      
            gas['name'],
            color=line.get_color(),            
            fontsize=9,
            fontweight='bold',
            verticalalignment='bottom',       
            horizontalalignment='center',     
            rotation=angulo_inclinacao,        
            rotation_mode='anchor'
        )
        i+=1


    for planet_name, group in df.groupby("planet"):
        plt.scatter(
            group["T_K"],
            group["planet_escape_velocity_km_s"],
            s=60,
            edgecolor="black",
            label=planet_name,
            zorder=5
        )

    plt.xscale("log")   
    plt.yscale("log")

    plt.xlim(T_min, T_max * 1.8)

    plt.xlabel("Temperatura de equilíbrio / efetiva [K]")
    plt.ylabel("Velocidade de escape [km/s]")
#quais gases são mais fáceis ou difíceis de reter?”
    plt.title("Retenção atmosférica preliminar: v_escape vs temperatura")

    adjust_text(texts, only_move={'text': 'y', 'objects': 'y'}, force_points=0.0)

    plt.legend(fontsize=9, loc="upper left")
    plt.tight_layout()
    plt.savefig("escape_velocity_vs_temperature_with_gases.png", dpi=300)
    plt.show()

    print("Gráfico salvo: escape_velocity_vs_temperature_with_gases.png")


if __name__ == "__main__":
    main()


#2) Artigo de ref: Damiano et al. 2024

#dado um planeta com massa, raio, densidade e temperatura de equilíbrio, quais atmosferas são fisicamente plausíveis?
#massa/raio/fluxo/T_eq
#→ planeta consegue reter atmosfera?
#→ quais famílias atmosféricas fazem sentido?
#→ N2/CO2/H2O? H2-rich? CO2-dominante? sem atmosfera?


#3) Artigo de Banerjee: possibilidade de atmosfera secundária. ele mostra que a temperatura atmosférica e o peso molecular médio são degenerados no espectro, porque a altura de escala depende de T/μ. Por isso eles impõem limites físicos na temperatura usando temperatura diurna teórica, albedo e eficiência de circulação. Isso pode ajudar depois, quando for definir intervalos de temperatura e pressão para atmosferas candidatas.

#4) qual cp, albedo e emissividade essa atmosfera tem? Eles dependem de superfície, nuvens, aerossóis, gelo, oceanos, composição atmosférica e absorção no infravermelho.
