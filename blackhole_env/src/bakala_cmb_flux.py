from __future__ import annotations

"""
Geometria Kerr em coordenadas de Boyer–Lindquist (G=c=M=1) e fluxo CMB
calibrado na zona habitável de Bakala

A integral bolométrica rigorosa do paper usa exclusão da sombra + zoom no
céu local (LSDPlus). Aqui foi feito:

  • coeficientes para o desvio de frequência g(χ,ψ) na órbita circular
    corotante Kepleriana (eqs. 14–21);
  • `integrate_cmb_flux_mollweide_with_shadow`: grade Mollweide + sombra por omissão
    via geodésicas nulas EinsteinPy (FANTASY), com fallback Carter - precisa do numpy
    e opcionalmente einsteinpy;
  • interpolação HZ (`hz_calibrated_flux_W_m2`) nos três pontos da Fig. 3 do Bakala.

Para valores “oficiais” da HZ, o melhor é usar a interpolação calibrada; `numeric_shadow`
serve para explorar (r, a) fora da tabela.
"""

import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from typing import Callable, Optional, Sequence, Tuple

from constants import SIGMA_SB, T_CMB


# Bakala et al. (2020), Sec. 4: pontos da HZ em órbita na ISCO para Φ tipo Kopparapu
# r em GM/c²; a adimensional; Φ em W m^-2 (Fig. 3).
# Cada item da lista de tuplas (sequencia imutavel dos valores) é (r, a, Φ). raio, spin, fluxo bolométrico
HZ_CALIBRATION: Tuple[Tuple[float, float, float], ...] = (
    (1.00090, 1.0 - 1.64e-10, 589.0),   # “Marte” limite frio
    (1.00057, 1.0 - 4.5e-11, 1366.0),   # “Terra”
    (1.00042, 1.0 - 1.75e-11, 2611.0),  # “Vênus” limite quente
)


# Σ=r^2+a^2cos^2(θ)
def sigma_bl(r: float, a: float, theta: float = math.pi / 2) -> float:
    return r * r + a * a * math.cos(theta) ** 2

# Δ=r^2−2r+a^2
def delta_bl(r: float, a: float) -> float:
    return r * r - 2.0 * r + a * a


# A=(r^2+a^2)^2−a^2Δ   puxa o delta_bl
# No plano equatorial, θ=π/2, então sin^2(θ)=1 e por isso aparece só: A=(r^2+a^2)^2−a^2Δ
def A_metric(r: float, a: float) -> float:
    return (r * r + a * a) ** 2 - a * a * delta_bl(r, a)

# eq 14 Bakala et al.
def omega_zamo_t_components(r: float, a: float, theta: float = math.pi / 2) -> Tuple[float, float, float]:
    """Componentes ω^(t)_t, ω^(φ)_t, ω^(φ)_φ na eq. (14), plano equatorial."""
    Sig = sigma_bl(r, a, theta)
    D = delta_bl(r, a)
    A = A_metric(r, a)
    if A <= 0.0 or Sig <= 0.0 or D <= 0.0:
        raise ValueError("Parâmetros métricos inválidos para ω ZAMO")
    w_tt = math.sqrt(D * Sig / A)
    w_phit = -2.0 * a * r / math.sqrt(A * Sig)
    w_phiphi = math.sqrt(A / Sig)
    w_rr = math.sqrt(Sig / D)
    w_theta_theta = math.sqrt(Sig)
    return w_tt, w_phit, w_phiphi



def beta_kepler_corotating(r: float, a: float) -> float:
    """Velocidade orbital Kepleriana corrotante medida pelo ZAMO (eq. 15)"""
    d = delta_bl(r, a)
    if d <= 0.0:
        raise ValueError(f"Δ não positivo em r={r}, a={a}")
    num = r * r + a * a - 2.0 * a * math.sqrt(r)
    den = math.sqrt(d) * (r ** 1.5 + a)
    return num / den

def XY_kepler_coefficients(r: float, a: float) -> Tuple[float, float, float]:
    """Retorna (X, Y, β) das eqs. (15), (21)."""
    beta = beta_kepler_corotating(r, a)
    w_tt, w_phit, _ = omega_zamo_t_components(r, a)
    X = w_tt - beta * w_phit
    Y = w_phit - beta * w_tt
    return X, Y, beta


# eq 20 e eq 4
# χ,ψ são angulos no ceu vistos pelo planeta
# a função g(χ,ψ) calcula esse blueshift para cada pedacinho do céu
def g_frequency_shift(chi: float, psi: float, r: float, a: float) -> float:
    """Fator de desvio de frequência g(χ,ψ) na eq. (20). χ e ψ em radianos."""
    X, Y, beta = XY_kepler_coefficients(r, a)
    den = X - Y * math.sin(chi) * math.cos(psi)
    if den <= 0.0:
        raise ValueError("Denominador não positivo em g(χ,ψ): direção não física ou singularidade")
    return math.sqrt(max(0.0, 1.0 - beta * beta)) / den


def hz_calibrated_flux_W_m2(r_orb_geom: float) -> float:
    """
    Fluxo bolométrico Φ interpolado entre os três pontos HZ do Bakala et al.
    (2020), função apenas de r_orb em GM/c² — válido nas órbitas que os autores
    associam à ISCO para cada Φ-alvo (Sol tipo Kopparapu).
    """
    pts = sorted(HZ_CALIBRATION, key=lambda x: x[0])
    if r_orb_geom <= pts[0][0]:
        return pts[0][2]
    if r_orb_geom >= pts[-1][0]: #-1: ultimo elemento da lista
        return pts[-1][2]

    for (r1, _a1, f1), (r2, _a2, f2) in zip(pts[:-1], pts[1:]): #[:-1]: todos menos o ultimo, [1:]: todos a partir do segundo
        if r1 <= r_orb_geom <= r2:
            t = (r_orb_geom - r1) / (r2 - r1)
            ln_f = math.log(f1) * (1.0 - t) + math.log(f2) * t
            return math.exp(ln_f)
    raise RuntimeError("interpolação HZ falhou")


# para um valor de r_orb, estima qual spin a estaria “na curva” da zona habitável
# do paper, interpolando entre os três pontos Marte–Terra–Vênus
# Apenas uma aproximação útil!
def suggested_spin_for_hz_orbit(r_orb_geom: float) -> float:
    """
    Spin a interpolado linearmente em r entre os pontos calibrados do paper,
    para órbitas na ISCO em cada extremo da HZ.
    """

    #ordem crescente de r_orb
    pts = sorted(HZ_CALIBRATION, key=lambda x: x[0])
    if r_orb_geom <= pts[0][0]:
        return pts[0][1]
    if r_orb_geom >= pts[-1][0]:
        return pts[-1][1]

    for (r1, a1, _), (r2, a2, _) in zip(pts[:-1], pts[1:]):
        if r1 <= r_orb_geom <= r2:
            # pega dois pontos vizinhos e calcula um spin intermediário
            t = (r_orb_geom - r1) / (r2 - r1)
            return a1 * (1.0 - t) + a2 * t
    raise RuntimeError("interpolação de spin falhou")


def equilibrium_temperature_black_body(phi_W_m2: float) -> float:
    """Eq. (3) Bakala et al.: T = (Φ/(4σ))^1/4 — usando planeta como corpo negro."""
    return (phi_W_m2 / (4.0 * SIGMA_SB)) ** 0.25


def time_dilation_gamma_kepler(r: float, a: float, theta: float = math.pi / 2) -> float:
    """Fator Γ na eq. (27): γ √(A/(ΔΣ))."""
    Sig = sigma_bl(r, a, theta)
    D = delta_bl(r, a)
    A = A_metric(r, a)
    _, _, beta = XY_kepler_coefficients(r, a)
    gamma_lorentz = 1.0 / math.sqrt(max(1e-30, 1.0 - beta * beta))
    return gamma_lorentz * math.sqrt(A / (D * Sig))


def _mollweide_chi_psi(x: float, y: float) -> Tuple[float, float]:
    """Projeção de Mollweide (eq. 22), x,y ∈ elipse [-1,1]×[-1/2,1/2]."""
    xi = math.asin(max(-1.0, min(1.0, 2.0 * y)))
    cos_xi = math.cos(xi)
    if abs(cos_xi) < 1e-14:
        psi = 0.0
    else:
        psi = math.pi * x / cos_xi
    inner = (2.0 * xi + math.sin(2.0 * xi)) / math.pi
    inner = max(-1.0, min(1.0, inner))
    chi = math.pi / 2.0 - math.asin(inner)
    return chi, psi


#integra só o céu que não está bloqueado pelo buraco negro
# pixels na sombra entram como peso 0 (descartados no somatorio, não contribui)
def integrate_cmb_flux_mollweide_with_shadow(
    r: float, #raio orbital
    a: float, #spin
    n: int = 40, #numero de pixels na grade/mapa dop céu
    *, # depois daqui: keyword-only arguments
    exclude_when_denominator_below: float = 1e-7,
    max_workers: Optional[int] = None, #q2uantos processos paralelos para calcular a sombra
    shadow_backend: str = "carter", #mudei pra carter, pois achei mais estável
) -> float:
    """
    Eq. (25) com exclusão de pixels na sombra.

    ``shadow_backend='einsteinpy'`` (predefinido): geodésicas nulas com EinsteinPy
    (FANTASY); perto do horizonte o traço em BL pode degradar e faz-se fallback
    para a máscara Carter. ``shadow_backend='carter'``: apenas R(r) Carter.

    Sem zoom adaptativo dos autores — Φ depende fortemente de n e pode
    desviar dos valores da Fig. 3. Para HZ calibrada use `hz_calibrated_flux_W_m2`.

    `max_workers`: paralelização por pixel (1 = sequencial).
    """
    warnings.warn(
        "integrate_cmb_flux_mollweide_with_shadow é exploratório: Φ varia com a "
        "grade e não reproduz o LSDPlus. Use hz_calibrated_flux_W_m2 para HZ "
        "conforme Bakala et al. (2020).",
        UserWarning,
        stacklevel=2,
    )

    if shadow_backend not in ("einsteinpy", "carter"):
        raise ValueError("shadow_backend deve ser 'einsteinpy' ou 'carter'")

# importa a função de calculo da sombra
    from kerr_shadow_ray import mollweide_xy_to_chi_psi, shadow_pixel_flux_weight

    if n < 12: #numero minimo de pixels para calcular a sombra
        raise ValueError("n muito pequeno para sombra")

    omega_pix = 8.0 / (n * n) #peso de cada pixel
    pref = SIGMA_SB * T_CMB**4 * omega_pix / math.pi

    X, Y, beta = XY_kepler_coefficients(r, a)

    tasks = []
    # monta a lista de pixels do céu (tasks)
    for i in range(1, n + 1): # coordenada y
        for j in range(1, 2 * n + 1): #coordenada x
            # cada parzinho x,y é o centro de um pixel da grade
            x = -1.0 + (j - 0.5) / n
            y = -0.5 + (i - 0.5) / n
            # verifica se o pixel está dentro da elipse do céu
            if (x / 1.0) ** 2 + (y / 0.5) ** 2 > 1.0 + 1e-12:
                continue
            # converte o pixel para coordenadas angulares chi e psi (direção nno céu)
            chi, psi = mollweide_xy_to_chi_psi(x, y)
            # adiciona o pixel à lista de tasks
            tasks.append(
                (chi, psi, r, a, X, Y, beta, pref, exclude_when_denominator_below, shadow_backend)
            )

    # escolhe quantos processos paralelos
    # se none e windows, usa 1 processo
    if max_workers is not None:
        workers = max(1, int(max_workers))
    elif sys.platform == "win32":
        workers = 1
    else:
        workers = min(8, os.cpu_count() or 4)

    if workers == 1:
        # soma de todas as contribuições para o fluxo no céu
        return sum(shadow_pixel_flux_weight(t) for t in tasks)

    with ProcessPoolExecutor(max_workers=workers) as pool:
        return sum(pool.map(shadow_pixel_flux_weight, tasks))


# integra “todo o céu brilhante pelo redshift” usando a integral da eq. (25)
# não remove a sombra
# superestima o fluxo, pois trata o céu como se não tivesse mascara de sombra
def integrate_cmb_flux_mollweide_naive(
    r: float, #raio orbital
    a: float, #spin
    n: int = 128, #numero de pixels na grade/mapa dop céu
    *, # depois daqui: keyword-only arguments
    exclude_when_denominator_below: float = 1e-6,
) -> float:
    """
    Estimativa da eq. (25) em grade Mollweide sem exclusão da sombra nem zoom
    adaptativo — pode divergir do Φ publicado (subestima ou superestima conforme
    grade). Útil para sanity-check qualitativo do mapa g⁴, não para HZ fina.

    Retorna Φ em W m^-2. Pixels com X − Y sinχ cosψ ≤ limiar são ignorados.
    """
    if n < 8:
        raise ValueError("n muito pequeno")

    omega_pix = 8.0 / (n * n)
    pref = SIGMA_SB * T_CMB**4 * omega_pix / math.pi

    total = 0.0
    X, Y, _beta = XY_kepler_coefficients(r, a)

    for i in range(1, n + 1):
        for j in range(1, 2 * n + 1):
            x = -1.0 + (j - 0.5) / n
            y = -0.5 + (i - 0.5) / n
            if (x / 1.0) ** 2 + (y / 0.5) ** 2 > 1.0 + 1e-12:
                continue
            chi, psi = _mollweide_chi_psi(x, y)
            den = X - Y * math.sin(chi) * math.cos(psi)
            if den <= exclude_when_denominator_below:
                continue
            g_val = math.sqrt(max(0.0, 1.0 - _beta * _beta)) / den
            total += pref * g_val**4

    return total


def make_constant_flux_callable(phi_W_m2: float) -> Callable[[float], float]:
    """Para acoplar ao climate0d.integrate_temperature (órbita circular)."""
    return lambda _t: phi_W_m2
