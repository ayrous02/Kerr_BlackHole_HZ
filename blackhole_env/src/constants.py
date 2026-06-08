from __future__ import annotations


# constants found at Bakala et al., 2020
SIGMA_SB = 5.670374419e-8       # Stefan–Boltzmann [W m^-2 K^-4]
T_CMB = 2.726                   # K
AU = 1.495978707e11             # m
G = 6.67430e-11                 # SI
C = 299_792_458.0               # m/s
M_SUN = 1.98847e30              # kg

# Fluxos de referência usados pra os limites, frio, medio, quente noi artigo de Bakala
FLUX_MARS = 589.0               # W m^-2
FLUX_EARTH = 1366.0             # W m^-2
FLUX_VENUS = 2611.0             # W m^-2


# Constantes para serem utilizadas na modelagem atmosférica
K_B = 1.380649e-23 
AMU_TO_KG = 1.66053906660e-27