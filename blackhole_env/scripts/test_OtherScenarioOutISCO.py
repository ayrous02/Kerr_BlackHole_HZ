from lsdplus_flux import OrbitMode, FluxIntegrationConfig, integrate_cmb_flux

# Ex.: 0,05% acima da ISCO, mesmo spin “Terra”
cfg = FluxIntegrationConfig(
    r_orb_geom=1.00110,
    spin_a=1.0 - 4.5e-11,
    orbit_mode=OrbitMode.FREE,
    n_grid=500,
    shadow_backend="carter",
)
res = integrate_cmb_flux(cfg)
print(res.phi_W_m2, res.isco_offset, res.on_isco)