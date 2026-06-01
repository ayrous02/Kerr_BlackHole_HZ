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
import csv
import math
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, List, Optional

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


def default_test_runs_root() -> Path:
    """Pasta base: blackhole_env/test_runs/."""
    return Path(__file__).resolve().parent.parent / "test_runs"


def create_run_output_dir(base_dir: Optional[Path] = None) -> Path:
    """Cria test_runs/AAAA-MM-DD_HH-MM-SS/ e devolve o caminho."""
    root = base_dir if base_dir is not None else default_test_runs_root()
    stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out = root / stamp
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_csv(path: Path, headers: List[str], rows: Iterable[Iterable[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        for row in rows:
            writer.writerow(row)


def save_run_metadata(out_dir: Path, args: argparse.Namespace) -> None:
    save_csv(
        out_dir / "run_config.csv",
        ["param", "value"],
        [[k, str(v)] for k, v in sorted(vars(args).items())],
    )
    save_csv(
        out_dir / "run_info.csv",
        ["key", "value"],
        [
            ["timestamp", datetime.now().isoformat(timespec="seconds")],
            ["python", sys.version.replace("\n", " ")],
            ["cwd", str(Path.cwd())],
        ],
    )


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


def calibrated_rows_to_csv(rows: list[Row]) -> tuple[List[str], List[List[Any]]]:
    headers = [
        "anchor",
        "r_orb",
        "spin_a",
        "r_plus",
        "r_isco",
        "isco_offset",
        "phi_ref_W_m2",
        "phi_calc_W_m2",
        "phi_rel_err",
        "phi_err_pct",
        "t_eq_K",
        "t_eq_C",
        "gamma",
        "roche_ratio",
        "stable",
    ]
    data = []
    for row in rows:
        data.append(
            [
                row.label,
                row.r_orb,
                row.spin_a,
                row.r_plus,
                row.r_isco,
                row.isco_offset,
                row.phi_ref,
                row.phi_calc,
                row.phi_rel_err,
                pct(row.phi_rel_err),
                row.t_eq_K,
                row.t_eq_C,
                row.gamma,
                row.roche_ratio,
                int(row.stable),
            ]
        )
    return headers, data


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


def run_scenario_consistency(mass_solar: float) -> tuple[List[str], List[List[Any]]]:
    headers = [
        "anchor",
        "phi_ref_W_m2",
        "phi_scenario_W_m2",
        "phi_rel_err",
        "phi_err_pct",
        "t_eq_K",
        "r_orb_m",
        "r_tidal_m",
        "roche_ratio",
        "survives_tidal",
        "spin_suggested",
        "spin_abs_delta",
    ]
    data: List[List[Any]] = []
    for label, r_orb, spin_a, phi_ref, _g_ref in PAPER_HZ_ANCHORS:
        scenario = BlackSunScenario(
            mass_solar=mass_solar,
            spin_a=spin_a,
            r_orb_geom=r_orb,
        )

        phi = scenario.flux_W_m2(mode="hz_calibrated")
        t_eq = scenario.equilibrium_temperature_K(phi)
        spin_sugg = scenario.suggested_spin_hz()
        err = rel_err(phi, phi_ref)

        data.append(
            [
                label,
                phi_ref,
                phi,
                err,
                pct(err),
                t_eq,
                scenario.orbit_radius_si_m(),
                scenario.tidal_disruption_radius_m(),
                scenario.roche_margin_ratio(),
                int(scenario.survives_tidal_approximation()),
                spin_sugg,
                abs(spin_a - spin_sugg),
            ]
        )
    return headers, data


def print_scenario_consistency(
    mass_solar: float,
    data: Optional[List[List[Any]]] = None,
) -> None:
    print("\n=== TESTE 2 — checagem do BlackSunScenario ===")
    print(
        "Confere se BlackSunScenario reproduz o fluxo calibrado, calcula marés "
        "e compara o spin sugerido pela interpolação.\n"
    )

    if data is None:
        _, data = run_scenario_consistency(mass_solar)
    table_rows = [
        [
            row[0],
            fmt_float(row[1], 3),
            fmt_float(row[2], 3),
            f"{row[4]:+.3e}%",
            fmt_float(row[5], 3),
            fmt_float(row[6], 3),
            fmt_float(row[7], 3),
            fmt_float(row[8], 6),
            "sim" if row[9] else "limite/não",
            f"{row[10]:.12f}",
            fmt_float(row[11], 3),
        ]
        for row in data
    ]

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
) -> tuple[List[str], List[List[Any]]]:
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

    headers = [
        "anchor",
        "phi_paper_W_m2",
        "phi_calibrated_W_m2",
        "phi_cal_err_pct",
        "phi_naive_W_m2",
        "phi_naive_err_pct",
        "phi_shadow_W_m2",
        "phi_shadow_err_pct",
        "t_naive_K",
        "t_shadow_K",
        "n_grid",
        "shadow_backend",
        "max_workers",
    ]
    data: List[List[Any]] = []
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

        data.append(
            [
                label,
                phi_ref,
                phi_cal,
                pct(rel_err(phi_cal, phi_ref)),
                phi_naive,
                pct(rel_err(phi_naive, phi_ref)),
                phi_shadow,
                pct(rel_err(phi_shadow, phi_ref)),
                equilibrium_temperature_black_body(phi_naive),
                equilibrium_temperature_black_body(phi_shadow),
                n_grid,
                shadow_backend,
                max_workers,
            ]
        )

    table_rows = [
        [
            row[0],
            fmt_float(row[1], 3),
            fmt_float(row[2], 3),
            f"{row[3]:+.3e}%",
            fmt_float(row[4], 3),
            f"{row[5]:+.3f}%",
            fmt_float(row[6], 3),
            f"{row[7]:+.3f}%",
            fmt_float(row[8], 3),
            fmt_float(row[9], 3),
        ]
        for row in data
    ]

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
    return headers, data


def run_lsdplus_diagnosis(
    n_grid: int, g_search_n: int, max_workers: int
) -> tuple[List[str], List[List[Any]]]:
    """
    Engenharia reversa do lsdplus_flux: estágios naive → sombra k=1 → pipeline completo.
    """
    print("\n=== TESTE 4a — diagnóstico LSDPlus (engenharia reversa) ===")
    print(
        "Compara Φ do paper com integral sem sombra, sombra k=1 e pipeline "
        "(zoom por convergência). Ajuda a ver onde o erro entra.\n"
    )
    print(f"n_grid={n_grid}, g_search_n={g_search_n}, max_workers={max_workers}\n")

    headers = [
        "anchor",
        "phi_paper_W_m2",
        "phi_naive_W_m2",
        "phi_naive_err_pct",
        "phi_shadow_k1_W_m2",
        "phi_shadow_k1_err_pct",
        "phi_pipeline_no_zoom_W_m2",
        "phi_pipeline_zoom_W_m2",
        "phi_pipeline_zoom_err_pct",
        "g_paper",
        "g_max_found",
        "g_err_pct",
        "zoom_k",
        "r_used",
        "n_grid",
        "g_search_n",
        "max_workers",
    ]
    data: List[List[Any]] = []
    for label, _r, _a, phi_ref, g_ref in PAPER_HZ_ANCHORS:
        d = diagnose_paper_anchor(
            label,
            n_grid=n_grid,
            g_search_n=g_search_n,
            max_workers=max_workers,
        )
        data.append(
            [
                d["anchor"],
                d["phi_paper"],
                d["phi_naive"],
                d["err_naive_pct"],
                d["phi_shadow_k1"],
                d["err_shadow_k1_pct"],
                d["phi_pipeline_no_zoom"],
                d["phi_pipeline_zoom"],
                d["err_pipeline_zoom_pct"],
                d["g_paper"],
                d["g_max_found"],
                d["err_g_pct"],
                d["zoom_k"],
                d["r_used"],
                n_grid,
                g_search_n,
                max_workers,
            ]
        )

    table_rows = [
        [
            row[0],
            fmt_float(row[1], 3),
            fmt_float(row[2], 3),
            f"{row[3]:+.1f}%",
            fmt_float(row[4], 3),
            f"{row[5]:+.1f}%",
            fmt_float(row[7], 3),
            f"{row[8]:+.1f}%",
            fmt_float(row[10], 0),
            fmt_float(row[11], 0),
            f"{row[12]:+.1f}%",
            fmt_float(row[13], 2),
        ]
        for row in data
    ]

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
        table_rows,
    )
    print(
        "\nLeitura rápida:\n"
        "- Φ naive ~ Φ shadow k=1: sombra Carter não é o único gargalo nesta grade.\n"
        "- Φ pipeline ≈ 2 W/m² com k zoom = 256: bug antigo de zoom (já corrigido; reexecute o teste).\n"
        "- Com código atual, pipeline ≈ shadow k=1 (~−50% vs paper em n=128).\n"
        "- Para se aproximar do paper na integral: tente --n-grid 500 (naive ~−14% na Terra).\n"
        "- g_found << g_paper: definição/discretização de g_max difere do LSDPlus."
    )
    return headers, data


def run_paper_zoom_tests(
    n_grid: int, g_search_n: int, max_workers: int, use_free_r: bool
) -> tuple[List[str], List[List[Any]]]:
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

    headers = [
        "anchor",
        "orbit_mode",
        "r_used",
        "r_isco",
        "isco_offset",
        "phi_paper_W_m2",
        "phi_computed_W_m2",
        "phi_err_pct",
        "g_paper",
        "g_computed",
        "g_err_pct",
        "t_eq_K",
        "t_eq_C",
        "zoom_k",
        "n_grid",
        "g_search_n",
        "max_workers",
        "use_free_r",
    ]
    data: List[List[Any]] = []

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

        data.append(
            [
                label,
                "FREE" if use_free_r else "ISCO",
                result.r_orb_used,
                result.r_isco,
                result.isco_offset,
                comp["phi_paper"],
                comp["phi_computed"],
                pct(comp["phi_rel_err"]),
                comp["g_paper"],
                comp["g_computed"],
                pct(comp["g_rel_err"]),
                result.temperature_K,
                temperature_C,
                result.zoom_k_opt,
                n_grid,
                g_search_n,
                max_workers,
                int(use_free_r),
            ]
        )

    table_rows = [
        [
            row[0],
            row[1],
            fmt_float(row[2], 8),
            fmt_float(row[3], 8),
            fmt_float(row[4], 3),
            fmt_float(row[5], 3),
            fmt_float(row[6], 3),
            f"{row[7]:+.3f}%",
            fmt_float(row[8], 3),
            fmt_float(row[9], 3),
            f"{row[10]:+.3f}%",
            fmt_float(row[11], 3),
            fmt_float(row[12], 3),
            fmt_float(row[13], 3),
        ]
        for row in data
    ]

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
            "T ºC",
            "zoom k",
        ],
        table_rows,
    )
    return headers, data


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
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="pasta base para CSV (padrão: blackhole_env/test_runs). Cria subpasta AAAA-MM-DD_HH-MM-SS.",
    )

    args = parser.parse_args()

    base_out = Path(args.output_dir) if args.output_dir else None
    out_dir = create_run_output_dir(base_out)
    save_run_metadata(out_dir, args)

    print("\nTeste de referência — cenário Black Sun / Bakala et al. (2020)")
    print(f"Massa usada para teste de marés: {args.mass_solar:.3e} M_sun")
    print(f"Resultados CSV: {out_dir}")

    rows = run_calibrated_tests(args.mass_solar)
    h1, d1 = calibrated_rows_to_csv(rows)
    save_csv(out_dir / "teste1_calibrado.csv", h1, d1)
    print_calibrated_report(rows)

    h2, d2 = run_scenario_consistency(args.mass_solar)
    save_csv(out_dir / "teste2_black_sun_scenario.csv", h2, d2)
    print_scenario_consistency(args.mass_solar, data=d2)

    print("\nInterpretação rápida:")
    print("- Se TESTE 1 e TESTE 2 batem em Φ, o caminho calibrado está funcionando.")
    print("- TESTE 3 (numeric) costuma divergir do paper: naive superestima; shadow depende da grade.")
    print("- R/Rt perto de 1 significa que a massa mínima do paper está no limite de marés.")
    print("- O paper_zoom (TESTE 4) é o mais ambicioso, mas também o mais sensível.")

    if args.run_numeric:
        h3, d3 = run_numeric_integral_tests(
            mass_solar=args.mass_solar,
            n_grid=args.numeric_n_grid,
            shadow_backend=args.shadow_backend,
            max_workers=args.max_workers,
        )
        save_csv(out_dir / "teste3_numeric.csv", h3, d3)
    else:
        print("\nTeste 3 (numeric) não foi executado. Para rodar:")
        print("python test_bakala_reference.py --run-numeric --numeric-n-grid 24")

    if args.run_lsdplus_diagnose:
        h4a, d4a = run_lsdplus_diagnosis(
            n_grid=args.n_grid,
            g_search_n=args.g_search_n,
            max_workers=args.max_workers,
        )
        save_csv(out_dir / "teste4a_lsdplus_diagnose.csv", h4a, d4a)

    if args.run_paper_zoom:
        h4, d4 = run_paper_zoom_tests(
            n_grid=args.n_grid,
            g_search_n=args.g_search_n,
            max_workers=args.max_workers,
            use_free_r=args.use_free_r,
        )
        save_csv(out_dir / "teste4_paper_zoom.csv", h4, d4)
    elif not args.run_lsdplus_diagnose:
        print("\nTestes LSDPlus não executados. Sugestão:")
        print("python test_bakala_reference.py --run-lsdplus-diagnose --n-grid 128 --g-search-n 512")
        print("python test_bakala_reference.py --run-paper-zoom --n-grid 500 --g-search-n 512")

    print(f"\nArquivos gravados em: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
