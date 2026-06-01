#!/usr/bin/env python3
"""
Valida a integração Φ contra os três pontos da HZ (Bakala et al. 2020, Fig. 3).

Uso (a partir de blackhole_env/):
  python scripts/validate_paper_hz.py
  python scripts/validate_paper_hz.py --anchor earth --n 500
  python scripts/validate_paper_hz.py --free --r 1.0006 --a 0.99999999995

Com n=500 e zoom Brent, cada ponto pode levar vários minutos (paralelização ajuda).
Para teste rápido: --n 64 --no-zoom-opt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from lsdplus_flux import (  # noqa: E402
    FluxIntegrationConfig,
    OrbitMode,
    PAPER_HZ_ANCHORS,
    compare_to_paper_anchor,
    config_for_paper_anchor,
    integrate_cmb_flux,
)


def main() -> int:
    p = argparse.ArgumentParser(description="Validar Φ CMB vs Bakala et al. (2020)")
    p.add_argument("--anchor", choices=["mars_cold", "earth", "venus_hot", "all"], default="earth")
    p.add_argument("--n", type=int, default=128, help="resolução Mollweide (paper: 500)")
    p.add_argument("--backend", choices=["carter", "einsteinpy"], default="carter")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--no-zoom-opt", action="store_true", help="céu inteiro k=1, sem Brent")
    p.add_argument("--free", action="store_true", help="órbita livre em vez de forçar ISCO")
    p.add_argument("--r", type=float, default=None, help="r_orb em GM/c² (modo --free)")
    p.add_argument("--a", type=float, default=None, help="spin a")
    args = p.parse_args()

    anchors = [a[0] for a in PAPER_HZ_ANCHORS] if args.anchor == "all" else [args.anchor]

    for name in anchors:
        if args.free and args.r is not None and args.a is not None:
            cfg = FluxIntegrationConfig(
                r_orb_geom=args.r,
                spin_a=args.a,
                orbit_mode=OrbitMode.FREE,
                n_grid=args.n,
                shadow_backend=args.backend,
                max_workers=args.workers,
                skip_zoom_optimization=args.no_zoom_opt,
            )
            label = f"free r={args.r} a={args.a}"
        else:
            cfg = config_for_paper_anchor(name, n_grid=args.n)
            cfg.shadow_backend = args.backend
            cfg.max_workers = args.workers
            cfg.skip_zoom_optimization = args.no_zoom_opt
            if args.free:
                cfg.orbit_mode = OrbitMode.FREE
            label = name

        print(f"\n=== {label} (n={cfg.n_grid}, backend={cfg.shadow_backend}) ===")
        res = integrate_cmb_flux(cfg)
        print(f"  r_used={res.r_orb_used:.8f}  r_ISCO={res.r_isco:.8f}  offset={res.isco_offset:+.2e}")
        print(f"  on_isco={res.on_isco}  k_zoom={res.zoom_k_opt:.2f}")
        print(f"  Phi={res.phi_W_m2:.1f} W/m^2  T={res.temperature_K - 273.15:.1f} C")
        print(f"  g_max={res.g_max:.0f}")

        if not args.free or (args.r is None):
            cmp_ = compare_to_paper_anchor(res, name)
            print(
                f"  vs paper: Phi {cmp_['phi_paper']:.0f} "
                f"(err {100 * cmp_['phi_rel_err']:+.1f}%)  "
                f"g {cmp_['g_paper']:.0f} (err {100 * cmp_['g_rel_err']:+.1f}%)"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
