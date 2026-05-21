# blackhole_env — CMB em órbitas Kerr (Bakala et al. 2020)

Reimplementação **aberta** do pipeline numérico de *Black Sun Revisited* (arXiv:2001.10991), sem o pacote proprietário LSDPlus (Bakala et al. 2015).

## O que está implementado

| Componente | Módulo |
|------------|--------|
| g(χ,ψ), tetrad ZAMO + Kepler (eqs. 14–21) | `bakala_cmb_flux.py` |
| Sombra Kerr (Carter R(r), Apêndice A) | `kerr_shadow_ray.py` |
| Ray tracing opcional | `einsteinpy_visibility.py` |
| **Integral eq. (25) + zoom Brent** | `lsdplus_flux.py` |
| Cenário + marés + Γ | `black_sun_scenario.py` |
| Tabela HZ interpolada (Fig. 3) | `HZ_CALIBRATION` em `bakala_cmb_flux.py` |

## Modos de órbita

- **`orbit_mode='isco'`** — reproduz o artigo: usa `r = r_ISCO(a)` para o spin dado.
- **`orbit_mode='free'`** — testa `r` mais perto ou mais longe da ISCO com o mesmo pipeline.

## Uso rápido

```python
import sys
sys.path.insert(0, "src")

from black_sun_scenario import BlackSunScenario
from lsdplus_flux import FluxIntegrationConfig, OrbitMode, integrate_cmb_flux

# Ponto “Terra” do paper (ISCO, spin quase 1)
cfg = FluxIntegrationConfig(
    r_orb_geom=1.00057,
    spin_a=1.0 - 4.5e-11,
    orbit_mode=OrbitMode.ISCO,
    n_grid=500,
    shadow_backend="carter",
)
res = integrate_cmb_flux(cfg)
print(res.phi_W_m2, res.temperature_K - 273.15)

# Órbita livre: 0,1% acima da ISCO
cfg_free = FluxIntegrationConfig(
    r_orb_geom=1.0010,
    spin_a=1.0 - 4.5e-11,
    orbit_mode=OrbitMode.FREE,
    n_grid=500,
)
```

Via `BlackSunScenario`:

```python
s = BlackSunScenario(mass_solar=1.63e8, spin_a=1 - 4.5e-11, r_orb_geom=1.00057)
phi = s.flux_W_m2(mode="paper_zoom", numeric_grid_n=500, orbit_mode="isco")
detail = s.flux_integration_paper(numeric_grid_n=500)
```

## Validação

```bash
cd blackhole_env
pip install -r requirements.txt
python scripts/validate_paper_hz.py --anchor earth --n 128
python scripts/validate_paper_hz.py --anchor earth --n 500   # reprodução completa (lento)
```

Teste fora da ISCO:

```bash
python scripts/validate_paper_hz.py --free --r 1.0012 --a 0.99999999995 --n 64
```

## Modos de fluxo em `flux_W_m2`

| `mode` | Descrição |
|--------|-----------|
| `hz_calibrated` | Interpolação nos 3 pontos da Fig. 3 (rápido, fiel aos números publicados) |
| `paper_zoom` | Integração numérica + zoom Brent (`lsdplus_flux`) |
| `numeric_shadow` | Integral sem zoom (exploratório) |
| `numeric_naive` | Sem sombra (sanity check) |

## Limitações honestas

- LSDPlus original não é público; esta é uma reimplementação com máscara Carter + Brent.
- `n=500` e zoom podem levar **minutos** por ponto; use `max_workers` ou reduza `n` para desenvolvimento.
- Espectro multibody (Figs. 5–6) ainda não está no repositório — próximo passo natural.

## Referência

Bakala, Dočekal & Turoňová 2020, ApJL 889, L37 — [arXiv:2001.10991](https://arxiv.org/abs/2001.10991)
