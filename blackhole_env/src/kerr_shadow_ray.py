from __future__ import annotations

"""
Máscara da sombra para o fluxo CMB (Bakala et al. 2020, Apêndice A).

Calcula constantes de movimento (l, q) a partir do tetrad local na órbita
Kepleriana corrotante e usa o polinómio radial de Carter R(r). Raízes reais
estritamente entre o horizonte externo r₊ e o raio orbital indicam barreira radial
incompatível com fótons vindos do infinito nesta parametrização equatorial.

É uma aproximação ao traço completo + zoom do código LSDPlus dos autores.
Requsito: numpy (numpy.roots).
"""

import math
from typing import Tuple

Vec4 = Tuple[float, float, float, float]


def sigma_bl(r: float, a: float, theta: float) -> float:
    return r * r + a * a * math.cos(theta) ** 2


def delta_bl(r: float, a: float) -> float:
    return r * r - 2.0 * r + a * a


def A_metric(r: float, a: float) -> float:
    return (r * r + a * a) ** 2 - a * a * delta_bl(r, a)


def omega_zamo_t_components(r: float, a: float, theta: float = math.pi / 2) -> Tuple[float, float, float]:
    Sig = sigma_bl(r, a, theta)
    D = delta_bl(r, a)
    A = A_metric(r, a)
    if A <= 0.0 or Sig <= 0.0 or D <= 0.0:
        raise ValueError("Parâmetros métricos inválidos para ω ZAMO")
    w_tt = math.sqrt(D * Sig / A)
    w_phit = -2.0 * a * r / math.sqrt(A * Sig)
    w_phiphi = math.sqrt(A / Sig)
    return w_tt, w_phit, w_phiphi


def covariant_p_from_local_sky(chi: float, psi: float, r_obs: float, a: float) -> Vec4:
    theta_obs = math.pi / 2.0
    Sig = sigma_bl(r_obs, a, theta_obs)
    D = delta_bl(r_obs, a)
    A = A_metric(r_obs, a)
    if D <= 0.0 or Sig <= 0.0 or A <= 0.0:
        raise ValueError("órbita inválida")

    w_tt, w_phit, w_phiphi = omega_zamo_t_components(r_obs, a, theta_obs)
    num = r_obs * r_obs + a * a - 2.0 * a * math.sqrt(r_obs)
    den = math.sqrt(D) * (r_obs ** 1.5 + a)
    beta = num / den
    gamma_l = 1.0 / math.sqrt(max(1e-30, 1.0 - beta * beta))

    wlat_t = gamma_l * (w_tt - beta * w_phit)
    wlat_ph = -gamma_l * beta * w_phiphi
    wlaphi_t = gamma_l * (w_phit - beta * w_tt)
    wlaphi_ph = gamma_l * w_phiphi

    p_lt, p_lr, p_lth, p_lph = -1.0, -math.sin(chi) * math.sin(psi), -math.cos(chi), -math.sin(chi) * math.cos(
        psi
    )

    p_t = wlat_t * p_lt + wlaphi_t * p_lph
    p_r = math.sqrt(Sig / D) * p_lr
    p_th = math.sqrt(Sig) * p_lth
    p_ph = wlat_ph * p_lt + wlaphi_ph * p_lph

    return (p_t, p_r, p_th, p_ph)


def carter_R_value(r: float, l: float, q: float, a: float) -> float:
    dlt = delta_bl(r, a)
    blk = r * r + a * a - a * l
    return blk * blk - dlt * (q + (l - a) ** 2)


def pixel_receives_cmb_from_infinity_carter(
    chi: float,
    psi: float,
    r_obs: float,
    a: float,
    *,
    outward_scan_points: int = 320,
    outward_r_factor: float = 120.0,
) -> bool:
    """
    Aproximação equatorial: exige R(r) ≥ 0 numa malha geométrica com r ≥ r_obs
    até grandes distâncias (troço exterior na retropropagação).

    Não equivale ao ray tracing + máscara de sombra completos do LSDPlus; Φ
    integrado depende fortemente da resolução Mollweide `n`. Para valores HZ da
    Fig. 3 de Bakala et al., use `hz_calibrated_flux_W_m2`.
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "numpy é necessário para a máscara Carter (pip install numpy)."
        ) from exc

    try:
        p_cov = covariant_p_from_local_sky(chi, psi, r_obs, a)
    except ValueError:
        return False

    pt, _pr, pth, pph = p_cov
    if abs(pt) < 1e-28:
        return False

    l = -pph / pt
    q = (pth / pt) ** 2

    r_plus = 1.0 + math.sqrt(max(0.0, 1.0 - a * a))
    r_obs_eff = max(r_obs, r_plus + 1e-10)

    if carter_R_value(r_obs_eff, l, q, a) <= 0.0:
        return False

    r_far = max(r_obs_eff * outward_r_factor, r_obs_eff + 400.0)
    rs = np.geomspace(r_obs_eff + 1e-10, r_far, outward_scan_points)
    vals = np.array([carter_R_value(float(rr), l, q, a) for rr in rs])
    return float(vals.min()) >= -1e-9


def mollweide_xy_to_chi_psi(x: float, y: float) -> Tuple[float, float]:
    xi = math.asin(max(-1.0, min(1.0, 2.0 * y)))
    cos_xi = math.cos(xi)
    if abs(cos_xi) < 1e-14:
        psi_w = 0.0
    else:
        psi_w = math.pi * x / cos_xi
    inner = (2.0 * xi + math.sin(2.0 * xi)) / math.pi
    inner = max(-1.0, min(1.0, inner))
    chi_w = math.pi / 2.0 - math.asin(inner)
    return chi_w, psi_w


def shadow_pixel_flux_weight(args: Tuple) -> float:
    """(chi, psi, r, a, X, Y, beta, pref, exclude_den) ou …, shadow_backend)."""
    chi, psi, r, a, X, Y, beta, pref, exclude_den = args[:9]
    backend = args[9] if len(args) > 9 else "einsteinpy"
    den = X - Y * math.sin(chi) * math.cos(psi)
    if den <= exclude_den:
        return 0.0
    if backend == "einsteinpy":
        try:
            from einsteinpy_visibility import pixel_receives_cmb_from_infinity_einsteinpy as _vis_ep

            visible = _vis_ep(chi, psi, r, a)
        except ImportError:
            visible = pixel_receives_cmb_from_infinity_carter(chi, psi, r, a)
    else:
        visible = pixel_receives_cmb_from_infinity_carter(chi, psi, r, a)
    if not visible:
        return 0.0
    g_val = math.sqrt(max(0.0, 1.0 - beta * beta)) / den
    return pref * g_val**4
