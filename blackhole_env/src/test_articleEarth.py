"""
Valida o ponto Terra (e opcionalmente os três âncoras HZ) vs Bakala et al. (2020).

No Windows é obrigatório rodar como script principal (spawn do multiprocessing):
  python test_articleEarth.py

Para teste rápido:
  python test_articleEarth.py --n 64 --no-zoom
"""

from __future__ import annotations

import argparse

from lsdplus_flux import (
    compare_to_paper_anchor,
    config_for_paper_anchor,
    integrate_cmb_flux,
)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=500, help="grade Mollweide (paper: 500)")
    p.add_argument("--no-zoom", action="store_true", help="sem otimização Brent em k")
    p.add_argument("--workers", type=int, default=None, help="processos paralelos (omitir = auto)")
    p.add_argument("--all-anchors", action="store_true", help="mars_cold, earth, venus_hot")
    args = p.parse_args()

    names = ["mars_cold", "earth", "venus_hot"] if args.all_anchors else ["earth"]

    for name in names:
        cfg = config_for_paper_anchor(name, n_grid=args.n)
        cfg.max_workers = args.workers
        cfg.skip_zoom_optimization = args.no_zoom
        print(f"\n=== {name} (n={cfg.n_grid}, zoom={not cfg.skip_zoom_optimization}) ===")
        result = integrate_cmb_flux(cfg)
        print(f"  Phi={result.phi_W_m2:.1f} W/m²  T={result.temperature_K - 273.15:.1f} °C")
        print(f"  g_max={result.g_max:.0f}  k_zoom={result.zoom_k_opt:.2f}")
        cmp_ = compare_to_paper_anchor(result, name)
        print(
            f"  vs paper: Phi {cmp_['phi_paper']:.0f} "
            f"(err {100 * cmp_['phi_rel_err']:+.1f}%)  "
            f"g {cmp_['g_paper']:.0f} (err {100 * cmp_['g_rel_err']:+.1f}%)"
        )


if __name__ == "__main__":
    main()
