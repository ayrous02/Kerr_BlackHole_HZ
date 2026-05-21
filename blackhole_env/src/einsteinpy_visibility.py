from __future__ import annotations

"""
Geodésicas nulas Kerr via EinsteinPy (integrador FANTASY — Christian & Chan 2021 ApJ).

O tetrad local em ``kerr_shadow_ray`` fornece ``p_μ`` covariante em Boyer–Lindquist;
elevamos com ``g^{μν}`` de EinsteinPy e integramos com ``GeodesicIntegrator``, passando o
4‑momento completo (evita a escolha de ramo da raiz em ``Nulllike`` / ``_P``).

Convenção numérica: usa-se ``p0 = -(g^{μν} p_ν)`` como direção afim na integração,
com melhor estabilidade em ``r`` BL para ``r_obs`` moderados. Perto do horizonte
(``r_obs ~ 1``) as coordenadas BL degradam; quando a fração de passos com ``r`` não
positivo/não finito é alta, faz-se fallback para a máscara Carter já existente.
"""


# antes do criterio de confiança, max(r) explodia para valores absurdos 
# não usar esse metodo. Pode fazer rodar para sempre.
import math
from typing import Tuple, Union

import numpy as np

from einsteinpy.geodesic.utils import _kerr
from einsteinpy.integrators import GeodesicIntegrator

from kerr_shadow_ray import covariant_p_from_local_sky, pixel_receives_cmb_from_infinity_carter


def _dn_val(x: Union[float, object]) -> float:
    return float(getattr(x, "val", x))


def contravariant_null_from_tetrad(
    chi: float,
    psi: float,
    r_obs: float,
    a: float,
    *,
    theta_obs: float,
) -> Tuple[np.ndarray, np.ndarray]:
    q0 = np.array([0.0, float(r_obs), float(theta_obs), 0.0])
    G = _kerr(q0.tolist(), a)
    guu = np.array([[_dn_val(G[i, j]) for j in range(4)] for i in range(4)], dtype=float)
    pc = np.array(covariant_p_from_local_sky(chi, psi, r_obs, a), dtype=float)
    p_up = guu @ pc
    return q0, p_up


def _outer_horizon_radius(a: float) -> float:
    return 1.0 + math.sqrt(max(0.0, 1.0 - a * a))


def _null_ray_stats(
    q0: np.ndarray,
    p0: np.ndarray,
    a: float,
    *,
    steps: int,
    delta: float,
    omega: float,
    order: int,
    r_target: float,
    r_plus: float,
    horizon_tol: float,
) -> Tuple[float, float, float]:
    """
    Retorna (max_r válido, min_r válido, fração de passos inválidos).
    Ignora apenas passos com r não finito ou r <= 0 (singularidade BL).
    Sai cedo se ``max_r >= r_target`` (escape) ou ``min_r < r_plus - horizon_tol`` (captura).
    """
    gi = GeodesicIntegrator(
        _kerr,
        (float(a),),
        q0,
        p0.astype(float),
        False,
        steps=max(int(steps), 2),
        delta=float(delta),
        omega=float(omega),
        suppress_warnings=True,
        order=int(order),
    )
    max_seen = float(q0[1])
    min_seen = float(q0[1])
    bad = 0
    n_steps = int(steps)
    for k in range(n_steps):
        gi.step()
        rr = _dn_val(gi.res_list[0][1])
        if not math.isfinite(rr) or rr <= 0.0:
            bad += 1
            continue
        max_seen = max(max_seen, rr)
        min_seen = min(min_seen, rr)
        denom_k = float(k + 1)
        if min_seen < r_plus - horizon_tol:
            return max_seen, min_seen, bad / denom_k
        if max_seen >= r_target:
            return max_seen, min_seen, bad / denom_k
    return max_seen, min_seen, bad / float(n_steps)


def pixel_receives_cmb_from_infinity_einsteinpy(
    chi: float,
    psi: float,
    r_obs: float,
    a: float,
    *,
    theta_obs: float = math.pi / 2,
    steps: int = 1200,
    delta: float = 0.2,
    omega: float = 1.0,
    order: int = 2,
    r_escape_factor: float = 120.0,
    min_escape_r_add: float = 60.0,
    max_bad_frac: float = 0.35,
    max_trusted_escape_r: float = 1.0e10,
    horizon_tol: float = 5e-5,
) -> bool:
    """
    True se o raio retropropagado atinge distância grande em BL sem mergulhar dentro
    do horizonte externo. Se a integração for numericamente pouco confiável, usa
    ``pixel_receives_cmb_from_infinity_carter``.
    """
    try:
        q0, p_up = contravariant_null_from_tetrad(chi, psi, r_obs, a, theta_obs=theta_obs)
    except ValueError:
        return False

    r_plus = _outer_horizon_radius(a)
    r_target = max(float(r_obs) * float(r_escape_factor), float(r_obs) + float(min_escape_r_add))

    # Ramo afim mais estável nos testes numéricos anteriores (retropropagação).
    p0 = -p_up

    max_r, min_r, bad_frac = _null_ray_stats(
        q0,
        p0,
        a,
        steps=steps,
        delta=delta,
        omega=omega,
        order=order,
        r_target=r_target,
        r_plus=r_plus,
        horizon_tol=horizon_tol,
    )

    trusted = (
        bad_frac <= max_bad_frac
        and max_r <= max_trusted_escape_r
        and min_r >= r_plus - horizon_tol
    )
    if not trusted:
        return pixel_receives_cmb_from_infinity_carter(chi, psi, r_obs, a)

    escaped = max_r >= r_target and min_r >= r_plus - horizon_tol
    return bool(escaped)
