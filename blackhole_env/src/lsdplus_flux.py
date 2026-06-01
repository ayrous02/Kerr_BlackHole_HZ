from __future__ import annotations

"""
Integração do fluxo CMB de acordo com Bakala et al. (2020) e os cálculos
do LSDPlus ( do mesmo autor, Bakala et al. 2015).

Implementa as eqs. (20)–(25) com exclusão da sombra (Apêndice A, máscara Carter)
e zoom adaptativo centrado no máximo de g(χ,ψ), com otimização de Brent na escala
de zoom k (Sec. 3 do paper).

Dois modos de órbita (ENUM):
  • ``free`` — (r_orb, spin_a) arbitrários (desde r ≥ r_ISCO corrotante);
  • ``isco`` — força r = r_ISCO(a) para reproduzir o protocolo da HZ do artigo.

``shadow_backend='carter'`` é o padrão (rápido, alinhado ao Apêndice A).
``einsteinpy`` opcional para testes de ray tracing. Não usar, não é a melhor opção
"""

import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple

from bakala_cmb_flux import (
    XY_kepler_coefficients,
    equilibrium_temperature_black_body,
    g_frequency_shift,
)
from constants import SIGMA_SB, T_CMB
from kerr_orbit import KerrOrbit
from kerr_shadow_ray import mollweide_xy_to_chi_psi, shadow_pixel_flux_weight


class OrbitMode(str, Enum):
    """Restrição orbital para a integração"""

    FREE = "free"
    ISCO = "isco"


# Bakala et al. (2020) Sec. 4 —> pontos A, B, C na ISCO (g_max do texto).
PAPER_HZ_ANCHORS: Tuple[Tuple[str, float, float, float, float], ...] = (
    ("mars_cold", 1.00090, 1.0 - 1.64e-10, 589.0, 11_450.0),
    ("earth", 1.00057, 1.0 - 4.5e-11, 1366.0, 18_260.0),
    ("venus_hot", 1.00042, 1.0 - 1.75e-11, 2611.0, 24_742.0),
)

# Corta direções onde X − Y sinχ cosψ → 0⁺ e g diverge (singularidade da eq. 20).
# 1e-7 permitia g ~ 10⁶ e Φ ~ 10¹³; 1e-3 aproxima g_max e Φ da ordem do paper.
DEFAULT_EXCLUDE_DENOMINATOR: float = 1e-3


@dataclass
class FluxIntegrationConfig:
    """Parâmetros da integração bolométrica Φ (eq. 25)."""

    r_orb_geom: float
    spin_a: float
    orbit_mode: OrbitMode = OrbitMode.FREE #ToDo: mudar para ISCO
    n_grid: int = 500
    shadow_backend: str = "carter"
    exclude_denominator_below: float = DEFAULT_EXCLUDE_DENOMINATOR
    max_workers: Optional[int] = None
    # Zoom adaptativo (Sec. 3): aumenta k até Φ estabilizar (não maximizar Φ).
    optimize_zoom: bool = True
    zoom_k_min: float = 1.0
    zoom_k_max: float = 256.0
    zoom_rel_tol: float = 0.05
    g_search_n: int = 512
    # Se False: uma passagem em céu inteiro (k=1), sem busca em k
    skip_zoom_optimization: bool = False


@dataclass
class FluxIntegrationResult:
    phi_W_m2: float
    temperature_K: float
    g_max: float
    chi_max: float
    psi_max: float
    x_center: float
    y_center: float
    zoom_k_opt: float
    r_orb_used: float
    spin_a_used: float
    orbit_mode: OrbitMode
    on_isco: bool
    r_isco: float
    isco_offset: float


def risco_corotating(spin_a: float) -> float:
    return KerrOrbit(spin=spin_a, r_orb=6.0).risco_corotating()


def resolve_orbit(
    r_orb_geom: float,
    spin_a: float,
    orbit_mode: OrbitMode,
) -> Tuple[float, float, bool, float]:
    """
    Retorna (r_usado, spin, on_isco, r_isco).
    Em modo ISCO substitui r pelo r_ISCO(a).
    """
    r_isco = risco_corotating(spin_a)
    if orbit_mode == OrbitMode.ISCO:
        r_used = r_isco
        on_isco = True
    else:
        r_used = r_orb_geom
        on_isco = abs(r_orb_geom - r_isco) <= max(1e-9, 1e-7 * r_isco)
    KerrOrbit(spin=spin_a, r_orb=r_used).check_stable()
    return r_used, spin_a, on_isco, r_isco


def chi_psi_to_mollweide_xy(chi: float, psi: float) -> Tuple[float, float]:
    """Inverte eq. (22): (χ,ψ) → (x,y) na elipse Mollweide."""
    target = math.sin(math.pi / 2.0 - chi)
    target = max(-1.0, min(1.0, target))

    lo, hi = 0.0, math.pi / 2.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        val = (2.0 * mid + math.sin(2.0 * mid)) / math.pi
        if val < target:
            lo = mid
        else:
            hi = mid
    xi = 0.5 * (lo + hi)
    y = math.sin(xi) / 2.0
    cos_xi = math.cos(xi)
    if abs(cos_xi) < 1e-14:
        x = 0.0
    else:
        x = psi * cos_xi / math.pi
    return x, y


def _in_ellipse(x: float, y: float) -> bool:
    return (x / 1.0) ** 2 + (y / 0.5) ** 2 <= 1.0 + 1e-12


def find_g_maximum(
    r: float,
    a: float,
    *,
    n: int = 128,
    shadow_backend: str = "carter",
    exclude_denominator_below: float = DEFAULT_EXCLUDE_DENOMINATOR,
) -> Tuple[float, float, float, float, float]:
    """
    Localiza o máximo de g(χ,ψ) na grade Mollweide (eq. 20), ignorando sombra
    apenas para localizar o centro do zoom (como no paper).
    Ignora pixels com denominador da eq. (20) abaixo do limiar (evita singularidade).
    Retorna (g_max, chi, psi, x, y).
    """
    X, Y, beta = XY_kepler_coefficients(r, a)
    beta_factor = math.sqrt(max(0.0, 1.0 - beta * beta))
    g_best = -1.0
    chi_b = psi_b = x_b = y_b = 0.0

    for i in range(1, n + 1):
        for j in range(1, 2 * n + 1):
            x = -1.0 + (j - 0.5) / n
            y = -0.5 + (i - 0.5) / n
            if not _in_ellipse(x, y):
                continue
            chi, psi = mollweide_xy_to_chi_psi(x, y)
            den = X - Y * math.sin(chi) * math.cos(psi)
            if den <= exclude_denominator_below:
                continue
            g_val = beta_factor / den
            if g_val > g_best:
                g_best = g_val
                chi_b, psi_b, x_b, y_b = chi, psi, x, y

    if g_best < 0.0:
        raise RuntimeError(f"não foi possível localizar g_max em r={r}, a={a}")
    return g_best, chi_b, psi_b, x_b, y_b


def _mollweide_pixel_tasks(
    r: float,
    a: float,
    n: int,
    *,
    x_center: float,
    y_center: float,
    zoom_k: float,
    shadow_backend: str,
    exclude_denominator_below: float,
) -> List[Tuple]:
    """Lista de argumentos para shadow_pixel_flux_weight."""
    X, Y, beta = XY_kepler_coefficients(r, a)
    k = max(float(zoom_k), 1.0)
    # Área sólida por pixel escala com o zoom (janela ∝ 1/k²); eq. (25) / Sec. 3.
    omega_pix = 8.0 / (n * n * k * k)
    pref = SIGMA_SB * T_CMB**4 * omega_pix / math.pi

    if k <= 1.0 + 1e-12:
        x_lo, x_hi = -1.0, 1.0
        y_lo, y_hi = -0.5, 0.5
    else:
        x_lo = x_center - 1.0 / k
        x_hi = x_center + 1.0 / k
        y_lo = y_center - 0.5 / k
        y_hi = y_center + 0.5 / k

    tasks: List[Tuple] = []
    for i in range(1, n + 1):
        for j in range(1, 2 * n + 1):
            x = x_lo + (x_hi - x_lo) * (j - 0.5) / n
            y = y_lo + (y_hi - y_lo) * (i - 0.5) / n
            if not _in_ellipse(x, y):
                continue
            chi, psi = mollweide_xy_to_chi_psi(x, y)
            tasks.append(
                (
                    chi,
                    psi,
                    r,
                    a,
                    X,
                    Y,
                    beta,
                    pref,
                    exclude_denominator_below,
                    shadow_backend,
                )
            )
    return tasks


def _resolve_max_workers(max_workers: Optional[int]) -> int:
    """No Windows (spawn), ProcessPool exige código sob ``if __name__ == '__main__'``."""
    if max_workers is not None:
        return max(1, int(max_workers))
    if sys.platform == "win32":
        return 1
    return min(8, os.cpu_count() or 4)


def _sum_pixel_weights(tasks: List[Tuple], max_workers: Optional[int]) -> float:
    if not tasks:
        return 0.0
    workers = _resolve_max_workers(max_workers)
    if workers == 1:
        return sum(shadow_pixel_flux_weight(t) for t in tasks)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return sum(pool.map(shadow_pixel_flux_weight, tasks))


def integrate_phi_mollweide(
    r: float,
    a: float,
    n: int,
    *,
    x_center: float = 0.0,
    y_center: float = 0.0,
    zoom_k: float = 1.0,
    shadow_backend: str = "carter",
    exclude_denominator_below: float = DEFAULT_EXCLUDE_DENOMINATOR,
    max_workers: Optional[int] = None,
) -> float:
    """Soma eq. (25) com sombra; ω_pix já inclui o fator 1/k² na discretização."""
    if shadow_backend not in ("carter", "einsteinpy"):
        raise ValueError("shadow_backend deve ser 'carter' ou 'einsteinpy'")
    if n < 12:
        raise ValueError("n_grid muito pequeno (paper usa n=500)")

    tasks = _mollweide_pixel_tasks(
        r,
        a,
        n,
        x_center=x_center,
        y_center=y_center,
        zoom_k=zoom_k,
        shadow_backend=shadow_backend,
        exclude_denominator_below=exclude_denominator_below,
    )
    return _sum_pixel_weights(tasks, max_workers)


def _phi_at_zoom_k(
    r: float,
    a: float,
    cfg: FluxIntegrationConfig,
    x_c: float,
    y_c: float,
    k: float,
) -> float:
    return integrate_phi_mollweide(
        r,
        a,
        cfg.n_grid,
        x_center=x_c,
        y_center=y_c,
        zoom_k=max(k, cfg.zoom_k_min),
        shadow_backend=cfg.shadow_backend,
        exclude_denominator_below=cfg.exclude_denominator_below,
        max_workers=cfg.max_workers,
    )


def _resolve_zoom_k_converged(
    r: float,
    a: float,
    cfg: FluxIntegrationConfig,
    x_c: float,
    y_c: float,
) -> Tuple[float, float]:
    """
    Varre k e mantém o maior Φ estável. Para quando Φ cai (<50% do melhor) ou
    diverge (>10×), pois janelas muito pequenas + máscara Carter zeram a soma.
    """
    ks = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, cfg.zoom_k_max]
    ks = sorted({k for k in ks if cfg.zoom_k_min <= k <= cfg.zoom_k_max})
    if not ks:
        ks = [cfg.zoom_k_min]

    k_best = ks[0]
    phi_best = _phi_at_zoom_k(r, a, cfg, x_c, y_c, k_best)

    for k in ks[1:]:
        phi_k = _phi_at_zoom_k(r, a, cfg, x_c, y_c, k)
        if phi_k > phi_best * 10.0:
            break
        if phi_k > phi_best:
            phi_best = phi_k
            k_best = k
            continue
        if phi_k < phi_best * 0.5:
            break
        rel = abs(phi_k - phi_best) / max(phi_best, 1e-30)
        if rel < cfg.zoom_rel_tol:
            break

    return k_best, phi_best


def integrate_cmb_flux(config: FluxIntegrationConfig) -> FluxIntegrationResult:
    """
    Pipeline principal: resolve órbita → g_max → (opcional) Brent em k → Φ, T.
    """
    r_used, a_used, on_isco, r_isco = resolve_orbit(
        config.r_orb_geom, config.spin_a, config.orbit_mode
    )

    g_max, chi_m, psi_m, x_m, y_m = find_g_maximum(
        r_used,
        a_used,
        n=config.g_search_n,
        shadow_backend=config.shadow_backend,
        exclude_denominator_below=config.exclude_denominator_below,
    )

    if config.skip_zoom_optimization or not config.optimize_zoom:
        k_opt = 1.0
        phi = integrate_phi_mollweide(
            r_used,
            a_used,
            config.n_grid,
            x_center=x_m,
            y_center=y_m,
            zoom_k=1.0,
            shadow_backend=config.shadow_backend,
            exclude_denominator_below=config.exclude_denominator_below,
            max_workers=config.max_workers,
        )
    else:
        k_opt, phi = _resolve_zoom_k_converged(r_used, a_used, config, x_m, y_m)

    return FluxIntegrationResult(
        phi_W_m2=phi,
        temperature_K=equilibrium_temperature_black_body(phi),
        g_max=g_max,
        chi_max=chi_m,
        psi_max=psi_m,
        x_center=x_m,
        y_center=y_m,
        zoom_k_opt=k_opt,
        r_orb_used=r_used,
        spin_a_used=a_used,
        orbit_mode=config.orbit_mode,
        on_isco=on_isco,
        r_isco=r_isco,
        isco_offset=r_used - r_isco,
    )


def config_for_paper_anchor(name: str, *, n_grid: int = 500) -> FluxIntegrationConfig:
    """Configuração nos três pontos A/B/C (órbita na ISCO)."""
    key = name.lower().strip()
    for label, r_o, a, _phi, _g in PAPER_HZ_ANCHORS:
        if label == key or key in label:
            return FluxIntegrationConfig(
                r_orb_geom=r_o,
                spin_a=a,
                orbit_mode=OrbitMode.ISCO,
                n_grid=n_grid,
            )
    raise ValueError(f"âncora HZ desconhecida: {name!r}; use mars_cold, earth, venus_hot")


def make_free_orbit_config(
    spin_a: float,
    r_orb_geom: float,
    *,
    n_grid: int = 500,
    shadow_backend: str = "carter",
    **kwargs,
) -> FluxIntegrationConfig:
    """
    Cenário exploratório: r arbitrário (ex.: ISCO ± δ).
    ``r_orb_geom`` deve ser ≥ r_ISCO(a).
    """
    return FluxIntegrationConfig(
        r_orb_geom=r_orb_geom,
        spin_a=spin_a,
        orbit_mode=OrbitMode.FREE,
        n_grid=n_grid,
        shadow_backend=shadow_backend,
        **kwargs,
    )


def diagnose_paper_anchor(
    name: str,
    *,
    n_grid: int = 128,
    g_search_n: int = 512,
    exclude_denominator_below: float = DEFAULT_EXCLUDE_DENOMINATOR,
    shadow_backend: str = "carter",
    max_workers: Optional[int] = None,
) -> dict:
    """
    Engenharia reversa: compara Φ paper com estágios do pipeline (naive, sombra, zoom).
    Útil para localizar onde a reimplementação diverge do LSDPlus.
    """
    from bakala_cmb_flux import integrate_cmb_flux_mollweide_naive

    cfg = config_for_paper_anchor(name, n_grid=n_grid)
    cfg.g_search_n = g_search_n
    cfg.exclude_denominator_below = exclude_denominator_below
    cfg.shadow_backend = shadow_backend
    cfg.max_workers = max_workers

    r_used, a_used, on_isco, r_isco = resolve_orbit(
        cfg.r_orb_geom, cfg.spin_a, cfg.orbit_mode
    )
    phi_paper = g_paper = None
    key = name.lower().strip()
    for label, _r_o, _a, phi_ref, g_ref in PAPER_HZ_ANCHORS:
        if label == key or key in label:
            phi_paper, g_paper = phi_ref, g_ref
            break
    if phi_paper is None:
        raise ValueError(f"âncora HZ desconhecida: {name!r}")

    g_max, chi_m, psi_m, x_m, y_m = find_g_maximum(
        r_used,
        a_used,
        n=g_search_n,
        exclude_denominator_below=exclude_denominator_below,
    )

    phi_naive = integrate_cmb_flux_mollweide_naive(
        r_used,
        a_used,
        n=n_grid,
        exclude_when_denominator_below=exclude_denominator_below,
    )
    phi_shadow_k1 = integrate_phi_mollweide(
        r_used,
        a_used,
        n_grid,
        x_center=x_m,
        y_center=y_m,
        zoom_k=1.0,
        shadow_backend=shadow_backend,
        exclude_denominator_below=exclude_denominator_below,
        max_workers=max_workers,
    )

    cfg_no_zoom = FluxIntegrationConfig(
        r_orb_geom=cfg.r_orb_geom,
        spin_a=cfg.spin_a,
        orbit_mode=cfg.orbit_mode,
        n_grid=n_grid,
        shadow_backend=shadow_backend,
        exclude_denominator_below=exclude_denominator_below,
        max_workers=max_workers,
        g_search_n=g_search_n,
        skip_zoom_optimization=True,
        optimize_zoom=False,
    )
    res_no_zoom = integrate_cmb_flux(cfg_no_zoom)

    cfg_zoom = FluxIntegrationConfig(
        r_orb_geom=cfg.r_orb_geom,
        spin_a=cfg.spin_a,
        orbit_mode=cfg.orbit_mode,
        n_grid=n_grid,
        shadow_backend=shadow_backend,
        exclude_denominator_below=exclude_denominator_below,
        max_workers=max_workers,
        g_search_n=g_search_n,
        skip_zoom_optimization=False,
        optimize_zoom=True,
    )
    res_zoom = integrate_cmb_flux(cfg_zoom)

    return {
        "anchor": name,
        "phi_paper": phi_paper,
        "g_paper": g_paper,
        "phi_naive": phi_naive,
        "phi_shadow_k1": phi_shadow_k1,
        "phi_pipeline_no_zoom": res_no_zoom.phi_W_m2,
        "phi_pipeline_zoom": res_zoom.phi_W_m2,
        "g_max_found": g_max,
        "zoom_k": res_zoom.zoom_k_opt,
        "r_used": r_used,
        "err_naive_pct": 100.0 * (phi_naive / phi_paper - 1.0),
        "err_shadow_k1_pct": 100.0 * (phi_shadow_k1 / phi_paper - 1.0),
        "err_pipeline_zoom_pct": 100.0 * (res_zoom.phi_W_m2 / phi_paper - 1.0),
        "err_g_pct": 100.0 * (g_max / g_paper - 1.0),
    }


def compare_to_paper_anchor(result: FluxIntegrationResult, name: str) -> dict:
    for label, r_o, a, phi_ref, g_ref in PAPER_HZ_ANCHORS:
        if label == name.lower().strip() or name.lower() in label:
            return {
                "anchor": label,
                "phi_paper": phi_ref,
                "phi_computed": result.phi_W_m2,
                "phi_rel_err": (result.phi_W_m2 - phi_ref) / phi_ref,
                "g_paper": g_ref,
                "g_computed": result.g_max,
                "g_rel_err": (result.g_max - g_ref) / g_ref,
                "T_C": result.temperature_K - 273.15,
            }
    raise ValueError(f"âncora desconhecida: {name}")
