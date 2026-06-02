from __future__ import annotations

ATMOSPHERE_CASES = [
    {
        "name": "bare_rock",
        "composition": "sem atmosfera ou atmosfera muito fina",
        "albedos": [0.05, 0.1, 0.2],
        "emissivities": [0.9, 1.0],
        "redistribution": [1, 2],
    },
    {
        "name": "N2_CO2_H2O",
        "composition": "N2 dominante, CO2 e H2O",
        "albedos": [0.2, 0.3, 0.4],
        "emissivities": [0.7, 0.85, 1.0],
        "redistribution": [2, 4],
    },
    {
        "name": "CO2_dominated",
        "composition": "CO2 dominante, N2 e H2O",
        "albedos": [0.3, 0.5, 0.7],
        "emissivities": [0.6, 0.8, 0.95],
        "redistribution": [2, 4],
    },
    {
        "name": "H2_rich",
        "composition": "H2 dominante, CO2, CH4 e H2O",
        "albedos": [0.3, 0.5, 0.7, 0.8],
        "emissivities": [0.5, 0.7, 0.9],
        "redistribution": [4],
    },
    {
        "name": "sulfur_volcanic",
        "composition": "CO2/N2 com SO2 e H2S",
        "albedos": [0.3, 0.5, 0.7],
        "emissivities": [0.6, 0.8, 0.95],
        "redistribution": [2, 4],
    },
]