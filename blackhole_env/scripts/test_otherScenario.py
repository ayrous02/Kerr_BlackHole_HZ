from black_sun_scenario import BlackSunScenario

s = BlackSunScenario(
    mass_solar=1.63e8,
    spin_a=1.0 - 4.5e-11,
    r_orb_geom=1.00057,  # ignorado se orbit_mode="isco"
)
detail = s.flux_integration_paper(numeric_grid_n=500, orbit_mode="isco")