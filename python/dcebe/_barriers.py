"""Smooth barrier functions for the GCV search domain.

Translates the inline ``barrier(x, b)`` defined in DCEBE_GCVscore.m. The
barrier is ``-log(sin((x-b) * pi/2))``, which is ``+inf`` at ``x = b``,
``0`` at ``x = b + 1``, and smoothly bridges the two. Outside the band
``[b, b+1]`` the function is clamped to ``+inf`` (below) or ``0`` (above)
so it can be added to any score without leaking into the feasible
region.

The score in DCEBE_GCVscore.m adds three barriers along three different
broadcast axes: ``beta`` along the j-axis, ``t`` along the i-axis (both
lower and upper bounds). The reshape semantics are preserved by
:func:`broadcast_barriers`.
"""
from __future__ import annotations

import numpy as np


def barrier(x, b: float) -> np.ndarray:
    """Smooth barrier at the lower bound ``b``.

    Returns ``-log(sin((x-b) * pi/2))`` clamped to ``+inf`` for ``x < b``
    and ``0`` for ``x > b + 1``.
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = -np.log(np.sin((x - b) * np.pi / 2.0))
    y = np.where(x < b, np.inf, y)
    y = np.where(x > b + 1, 0.0, y)
    return y


def broadcast_barriers(
    score: np.ndarray,
    t_arr: np.ndarray,
    beta_arr: np.ndarray,
    N: int,
    k: int,
    min_beta: float,
) -> np.ndarray:
    """Add the three barriers (beta lower, t lower, t upper) to a
    ``(M, T, A)`` score array, with the broadcast axes that match
    DCEBE_GCVscore.m.

    Parameters
    ----------
    score : ndarray, shape ``(M, T, A)``
        GCV scores before barrier addition.
    t_arr : ndarray, shape ``(T,)``
        Coarse-grid BAT candidates (1-based).
    beta_arr : ndarray, shape ``(A,)``
        Stiffness candidates.
    N : int
        Signal length.
    k : int
        Spline order.
    min_beta : float
        Lower bound for the beta barrier (typically ``beta_arr[0]``).

    Returns
    -------
    score : ndarray, shape ``(M, T, A)``
        Score plus the three barriers, broadcast in place.
    """
    beta_b = barrier(beta_arr, min_beta - 1.0)        # shape (A,)
    t_lo = barrier(t_arr, 1.0)                        # shape (T,)
    t_hi = barrier(-t_arr, -(N - k))                  # shape (T,)
    return (
        score
        + beta_b.reshape(1, 1, -1)
        + t_lo.reshape(1, -1, 1)
        + t_hi.reshape(1, -1, 1)
    )
