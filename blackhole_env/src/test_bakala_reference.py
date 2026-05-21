from __future__ import annotations

"""
test_bakala_reference.py

Arquivo de teste para comparar os módulos do cenário Black Sun com os valores
de referência do artigo de Bakala et al. (2020).

Para comparar integrais numéricas (sem zoom LSDPlus):
    python test_bakala_reference.py --run-numeric --numeric-n-grid 24

Para rodar também a integração mais fiel/experimental via LSDPlusFlux:
    python test_bakala_reference.py --run-paper-zoom --n-grid 80 --g-search-n 64

Observação:
- O modo "calibrated" usa diretamente os três pontos do artigo e deve bater quase
  exatamente por construção.
- numeric_naive / numeric_shadow são exploratórios (grade fixa k=1, sem zoom Brent).
- O modo "paper_zoom" tenta reproduzir a integração numérica do artigo; pode ser
  lento e pode não bater enquanto a máscara da sombra/zoom não estiver validada.
"""

import argparse
import math
import warnings
from dataclasses import dataclass
from typing import Iterable

from bakala_cmb_flux import (
    equilibrium_temperature_black_body,
    hz_calibrated_flux_W_m2,
    suggested_spin_for_hz_orbit,
)
from black_sun_scenario import BlackSunScenario, REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL
from kerr_orbit import KerrOrbit
from lsdplus_flux import (
    PAPER_HZ_ANCHORS,
    FluxIntegrationConfig,
    OrbitMode,
    compare_to_paper_anchor,
    diagnose_paper_anchor,
    integrate_cmb_flux,
)


# por enquanto apenas o Teste1 usa a classe Row
@dataclass(frozen=True)
class Row:
    label: str
    r_orb: float
    spin_a: float
    r_plus: float
    r_isco: float
    isco_offset: float
    phi_ref: float
    phi_calc: float
    phi_rel_err: float
    t_eq_K: float
    t_eq_C: float
    gamma: float
    roche_ratio: float
    stable: bool


def rel_err(value: float, ref: float) -> float:
    if ref == 0:
        return float("nan")
    return (value - ref) / ref


def pct(value: float) -> float:
    return 100.0 * value


def fmt_float(x: float, digits: int = 6) -> str:
    if not math.isfinite(x):
        return str(x)
    if abs(x) >= 1e5 or (0 < abs(x) < 1e-4):
        return f"{x:.{digits}e}"
    return f"{x:.{digits}f}"


def print_table(headers: list[str], rows: Iterable[list[str]]) -> None:
    rows = list(rows)
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def line(parts: list[str]) -> str:
        return " | ".join(part.ljust(widths[i]) for i, part in enumerate(parts))

    print(line(headers))
    print("-+-".join("-" * w for w in widths))
    for row in rows:
        print(line(row))


def run_calibrated_tests(mass_solar: float) -> list[Row]:
    """
    Teste rápido: usa os três pontos do artigo em hz_calibrated_flux_W_m2.
    Este teste deve reproduzir Φ de Marte/Terra/Vênus por construção (os fluxos calibrados)
    """
    rows: list[Row] = []

    for label, r_orb, spin_a, phi_ref, _g_ref in PAPER_HZ_ANCHORS:
        orbit = KerrOrbit(spin=spin_a, r_orb=r_orb)
        r_plus = orbit.r_plus
        r_isco = orbit.risco_corotating()
        isco_offset = r_orb - r_isco

        stable = True
        try:
            orbit.check_stable()
        except ValueError:
            stable = False

        phi_calc = hz_calibrated_flux_W_m2(r_orb)
        t_eq_K = equilibrium_temperature_black_body(phi_calc)

        scenario = BlackSunScenario(
            mass_solar=mass_solar,
            spin_a=spin_a,
            r_orb_geom=r_orb,
        )

        rows.append(
            Row(
                label=label,
                r_orb=r_orb,
                spin_a=spin_a,
                r_plus=r_plus,
                r_isco=r_isco,
                isco_offset=isco_offset,
                phi_ref=phi_ref,
                phi_calc=phi_calc,
                phi_rel_err=rel_err(phi_calc, phi_ref),
                t_eq_K=t_eq_K,
                t_eq_C=t_eq_K - 273.15,
                gamma=scenario.time_dilation_gamma(),
                roche_ratio=scenario.roche_margin_ratio(),
                stable=stable,
            )
        )

    return rows


def print_calibrated_report(rows: list[Row]) -> None:
    print("\n=== TESTE 1 — valores calibrados da Fig. 3 / pontos A, B, C ===")
    print("Este teste usa hz_calibrated_flux_W_m2(r). Deve bater com Φ por construção.\n")

    print_table(
        [
            "âncora",
            "r_orb",
            "a",
            "r+",
            "r_ISCO",
            "r-r_ISCO",
            "Φ ref",
            "Φ calc",
            "erro Φ",
            "T_eq (K)",
            "T_eq (°C)",
            "Γ",
            "R/Rt",
            "estável?",
        ],
        [
            [
                row.label,
                fmt_float(row.r_orb, 8),
                f"{row.spin_a:.12f}",
                fmt_float(row.r_plus, 8),
                fmt_float(row.r_isco, 8),
                fmt_float(row.isco_offset, 3),
                fmt_float(row.phi_ref, 3),
                fmt_float(row.phi_calc, 3),
                f"{pct(row.phi_rel_err):+.3e}%",
                fmt_float(row.t_eq_K, 3),
                fmt_float(row.t_eq_C, 3),
                fmt_float(row.gamma, 3),
                fmt_float(row.roche_ratio, 6),
                "sim" if row.stable else "não",
            ]
            for row in rows
        ],
    )


def print_scenario_consistency(mass_solar: float) -> None:
    print("\n=== TESTE 2 — checagem do BlackSunScenario ===")
    print(
        "Confere se BlackSunScenario reproduz o fluxo calibrado, calcula marés "
        "e compara o spin sugerido pela interpolação.\n"
    )

    table_rows = []
    for label, r_orb, spin_a, phi_ref, _g_ref in PAPER_HZ_ANCHORS:
        scenario = BlackSunScenario(
            mass_solar=mass_solar,
            spin_a=spin_a,
            r_orb_geom=r_orb,
        )

        phi = scenario.flux_W_m2(mode="hz_calibrated")
        t_eq = scenario.equilibrium_temperature_K(phi)
        spin_sugg = scenario.suggested_spin_hz()

        table_rows.append(
            [
                label,
                fmt_float(phi_ref, 3),
                fmt_float(phi, 3),
                f"{pct(rel_err(phi, phi_ref)):+.3e}%",
                fmt_float(t_eq, 3),
                fmt_float(scenario.orbit_radius_si_m(), 3),
                fmt_float(scenario.tidal_disruption_radius_m(), 3),
                fmt_float(scenario.roche_margin_ratio(), 6),
                "sim" if scenario.survives_tidal_approximation() else "limite/não",
                f"{spin_sugg:.12f}",
                fmt_float(abs(spin_a - spin_sugg), 3),
            ]
        )

    print_table(
        [
            "âncora",
            "Φ ref",
            "Φ cenário",
            "erro",
            "T_eq K",
            "R_orb m",
            "R_tidal m",
            "R/Rt",
            "marés?",
            "a sugerido",
            "|a-a_sug|",
        ],
        table_rows,
    )


def run_numeric_integral_tests(
    mass_solar: float,
    n_grid: int,
    shadow_backend: str,
    max_workers: int,
) -> None:
    """
    Compara Φ da Fig. 3 com integrais Mollweide:
    numeric_naive (sem sombra) e numeric_shadow (com máscara).
    """
    print("\n=== TESTE 3 — integrais numéricas: naive vs shadow (sem zoom Brent) ===")
    print(
        "Compara os três âncoras HZ com Φ do paper. Espera-se divergência forte "
        "em relação ao calibrado; serve para diagnóstico do pipeline numérico.\n"
    )
    print(
        f"n_grid={n_grid}, shadow_backend={shadow_backend}, max_workers={max_workers}\n"
    )

    if shadow_backend not in ("carter", "einsteinpy"):
        raise ValueError("shadow_backend deve ser 'carter' ou 'einsteinpy'")

    table_rows = []
    for label, r_orb, spin_a, phi_ref, _g_ref in PAPER_HZ_ANCHORS:
        scenario = BlackSunScenario(
            mass_solar=mass_solar,
            spin_a=spin_a,
            r_orb_geom=r_orb,
        )
        phi_cal = scenario.flux_W_m2(mode="hz_calibrated")
        phi_naive = scenario.flux_W_m2(mode="numeric_naive", numeric_grid_n=n_grid)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            phi_shadow = scenario.flux_W_m2(
                mode="numeric_shadow",
                numeric_grid_n=n_grid,
                shadow_backend=shadow_backend,
                shadow_max_workers=max_workers,
            )

        table_rows.append(
            [
                label,
                fmt_float(phi_ref, 3),
                fmt_float(phi_cal, 3),
                f"{pct(rel_err(phi_cal, phi_ref)):+.3e}%",
                fmt_float(phi_naive, 3),
                f"{pct(rel_err(phi_naive, phi_ref)):+.3f}%",
                fmt_float(phi_shadow, 3),
                f"{pct(rel_err(phi_shadow, phi_ref)):+.3f}%",
                fmt_float(equilibrium_temperature_black_body(phi_naive), 3),
                fmt_float(equilibrium_temperature_black_body(phi_shadow), 3),
            ]
        )

    print_table(
        [
            "âncora",
            "Φ paper",
            "Φ calibrado",
            "erro cal",
            "Φ naive",
            "erro naive",
            "Φ shadow",
            "erro shadow",
            "T naive K",
            "T shadow K",
        ],
        table_rows,
    )


def run_lsdplus_diagnosis(n_grid: int, g_search_n: int, max_workers: int) -> None:
    """
    Engenharia reversa do lsdplus_flux: estágios naive → sombra k=1 → pipeline completo.
    """
    print("\n=== TESTE 4a — diagnóstico LSDPlus (engenharia reversa) ===")
    print(
        "Compara Φ do paper com integral sem sombra, sombra k=1 e pipeline "
        "(zoom por convergência). Ajuda a ver onde o erro entra.\n"
    )
    print(f"n_grid={n_grid}, g_search_n={g_search_n}, max_workers={max_workers}\n")

    rows = []
    for label, _r, _a, phi_ref, g_ref in PAPER_HZ_ANCHORS:
        d = diagnose_paper_anchor(
            label,
            n_grid=n_grid,
            g_search_n=g_search_n,
            max_workers=max_workers,
        )
        rows.append(
            [
                label,
                fmt_float(phi_ref, 3),
                fmt_float(d["phi_naive"], 3),
                f"{d['err_naive_pct']:+.1f}%",
                fmt_float(d["phi_shadow_k1"], 3),
                f"{d['err_shadow_k1_pct']:+.1f}%",
                fmt_float(d["phi_pipeline_zoom"], 3),
                f"{d['err_pipeline_zoom_pct']:+.1f}%",
                fmt_float(g_ref, 0),
                fmt_float(d["g_max_found"], 0),
                f"{d['err_g_pct']:+.1f}%",
                fmt_float(d["zoom_k"], 2),
            ]
        )

    print_table(
        [
            "âncora",
            "Φ paper",
            "Φ naive",
            "erro naive",
            "Φ shadow k=1",
            "erro sh.k=1",
            "Φ pipeline",
            "erro pipe",
            "g paper",
            "g found",
            "erro g",
            "k zoom",
        ],
        rows,
    )
    print(
        "\nLeitura rápida: se Φ naive ≈ paper mas Φ shadow << paper, o gargalo é a máscara Carter. "
        "Se g_found >> g_paper, refine exclude_denominator ou g_search_n."
    )


def run_paper_zoom_tests(n_grid: int, g_search_n: int, max_workers: int, use_free_r: bool) -> None:
    """
    Teste experimental: tenta recalcular Φ com integração de céu local + sombra + zoom.
    Pode ser lento. Use n_grid pequeno para teste rápido; n_grid=500 se quiser tentar
    aproximar o protocolo do artigo.
    """
    print("\n=== TESTE 4 — integração LSDPlusFlux / paper_zoom experimental ===")
    print(
        "Este teste NÃO é a interpolação calibrada. Ele tenta recomputar Φ integrando "
        "o céu local com sombra e zoom. Etapa em progresso....\n"
    )
    print(f"n_grid={n_grid}, g_search_n={g_search_n}, max_workers={max_workers}\n")

    table_rows = []

    for label, r_orb, spin_a, _phi_ref, _g_ref in PAPER_HZ_ANCHORS:
        cfg = FluxIntegrationConfig(
            r_orb_geom=r_orb,
            spin_a=spin_a,
            orbit_mode=OrbitMode.FREE if use_free_r else OrbitMode.ISCO,
            n_grid=n_grid,
            shadow_backend="carter",
            max_workers=max_workers,
            g_search_n=g_search_n,
            optimize_zoom=True,
            skip_zoom_optimization=False,
        )

        result = integrate_cmb_flux(cfg)
        comp = compare_to_paper_anchor(result, label)
        temperature_C = result.temperature_K - 273.15

        table_rows.append(
            [
                label,
                "FREE r_paper" if use_free_r else "ISCO(a)",
                fmt_float(result.r_orb_used, 8),
                fmt_float(result.r_isco, 8),
                fmt_float(result.isco_offset, 3),
                fmt_float(comp["phi_paper"], 3),
                fmt_float(comp["phi_computed"], 3),
                f"{pct(comp['phi_rel_err']):+.3f}%",
                fmt_float(comp["g_paper"], 3),
                fmt_float(comp["g_computed"], 3),
                f"{pct(comp['g_rel_err']):+.3f}%",
                fmt_float(result.temperature_K, 3),
                fmt_float(temperature_C, 3),
                fmt_float(result.zoom_k_opt, 3),
            ]
        )

    print_table(
        [
            "âncora",
            "modo órbita",
            "r usado",
            "r_ISCO",
            "offset",
            "Φ paper",
            "Φ calc",
            "erro Φ",
            "g paper",
            "g calc",
            "erro g",
            "T K",
            "T ºC"
            "zoom k",
        ],
        table_rows,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mass-solar",
        type=float,
        default=REFERENCE_MIN_MASS_SOLAR_FOR_TIDAL,
        help="massa do buraco negro em massas solares",
    )
    parser.add_argument(
        "--run-numeric",
        action="store_true",
        help="compara Φ paper vs numeric_naive e numeric_shadow (sem zoom Brent)",
    )
    parser.add_argument(
        "--numeric-n-grid",
        type=int,
        default=24,
        help="resolução Mollweide para numeric_naive/shadow (mín. 12 com sombra)",
    )
    parser.add_argument(
        "--shadow-backend",
        choices=["carter", "einsteinpy"],
        default="carter",
        help="máscara de sombra no numeric_shadow",
    )
    parser.add_argument(
        "--run-lsdplus-diagnose",
        action="store_true",
        help="tabela de engenharia reversa (naive / sombra / zoom vs paper)",
    )
    parser.add_argument(
        "--run-paper-zoom",
        action="store_true",
        help="também roda a integração mais pesada/experimental LSDPlusFlux",
    )
    parser.add_argument(
        "--n-grid",
        type=int,
        default=80,
        help="resolução da grade Mollweide para paper_zoom. Paper usa ~500.",
    )
    parser.add_argument(
        "--g-search-n",
        type=int,
        default=64,
        help="grade para procurar g_max.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="número de processos. No Windows, comece com 1.",
    )
    parser.add_argument(
        "--use-free-r",
        action="store_true",
        help="no paper_zoom, usa o r_orb impresso no artigo em vez de forçar r=r_ISCO(a).",
    )

    args = parser.parse_args()

    print("\nTeste de referência — cenário Black Sun / Bakala et al. (2020)")
    print(f"Massa usada para teste de marés: {args.mass_solar:.3e} M_sun")

    rows = run_calibrated_tests(args.mass_solar)
    print_calibrated_report(rows)
    print_scenario_consistency(args.mass_solar)

    print("\nInterpretação rápida:")
    print("- Se TESTE 1 e TESTE 2 batem em Φ, o caminho calibrado está funcionando.")
    print("- TESTE 3 (numeric) costuma divergir do paper: naive superestima; shadow depende da grade.")
    print("- R/Rt perto de 1 significa que a massa mínima do paper está no limite de marés.")
    print("- O paper_zoom (TESTE 4) é o mais ambicioso, mas também o mais sensível.")

    if args.run_numeric:
        run_numeric_integral_tests(
            mass_solar=args.mass_solar,
            n_grid=args.numeric_n_grid,
            shadow_backend=args.shadow_backend,
            max_workers=args.max_workers,
        )
    else:
        print("\nTeste 3 (numeric) não foi executado. Para rodar:")
        print("python test_bakala_reference.py --run-numeric --numeric-n-grid 24")

    if args.run_lsdplus_diagnose:
        run_lsdplus_diagnosis(
            n_grid=args.n_grid,
            g_search_n=args.g_search_n,
            max_workers=args.max_workers,
        )

    if args.run_paper_zoom:
        run_paper_zoom_tests(
            n_grid=args.n_grid,
            g_search_n=args.g_search_n,
            max_workers=args.max_workers,
            use_free_r=args.use_free_r,
        )
    elif not args.run_lsdplus_diagnose:
        print("\nTestes LSDPlus não executados. Sugestão:")
        print("python test_bakala_reference.py --run-lsdplus-diagnose --n-grid 128 --g-search-n 512")
        print("python test_bakala_reference.py --run-paper-zoom --n-grid 500 --g-search-n 512")


if __name__ == "__main__":
    main()
