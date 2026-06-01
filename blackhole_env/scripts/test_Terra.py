import sys
sys.path.insert(0, "blackhole_env/src")

from lsdplus_flux import config_for_paper_anchor, integrate_cmb_flux

cfg = config_for_paper_anchor("earth", n_grid=500)  # ISCO, spin do paper
res = integrate_cmb_flux(cfg)
# Esperado: Phi ~ 1366 W/m², T ~ 5 °C, g_max ~ 18260