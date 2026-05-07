"""Top-level driver for bolus arrival time estimation.

Skeleton for day 2: argument validation, defaults resolution, and the
:class:`EstimateResult` container. The coarse-then-fine optimisation
loop lands in day 3.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Sequence

import numpy as np

from ._search import default_interval

_DEFAULT_BETA_COARSE = np.linspace(1.0, 25.0, 50)
_DEFAULT_ORDERS = (3, 4, 5, 6)


@dataclass
class EstimateResult:
    """Container for :func:`estimate_bat` output.

    Attributes
    ----------
    bat : ndarray, shape ``(M,)``
        Estimated bolus arrival time per signal, in 1-based sample units
        (matching the MATLAB reference). Convert to seconds via
        ``(bat - 1) * delta_t``.
    cp_int : ndarray, shape ``(M,)``
        ``floor(bat)``: the integer changepoint index.
    beta_opt : ndarray, shape ``(M,)``
        Optimal stiffness parameter per signal.
    k_opt : ndarray, shape ``(M,)``
        Optimal spline order per signal.
    score_opt : ndarray, shape ``(M,)``
        Final log-GCV score at the optimum.
    y_hat : ndarray, shape ``(N, M)``
        Reconstructed signal at the optimum.
    y_hat_x : ndarray, shape ``(N, M)``
        x-axis sample positions, with the row of ``floor(bat)`` replaced
        by the fractional ``bat`` so plots show the true changepoint.
    """

    bat: np.ndarray
    cp_int: np.ndarray
    beta_opt: np.ndarray
    k_opt: np.ndarray
    score_opt: np.ndarray
    y_hat: np.ndarray
    y_hat_x: np.ndarray


def _validate_inputs(
    y: np.ndarray,
    search_interval: tuple[float, float] | None,
    coarse_res: float,
    beta_coarse_search: np.ndarray | None,
    orders: Sequence[int],
    common_bat: bool,
    solver: str,
    verbosity: int,
) -> tuple[np.ndarray, tuple[float, float], np.ndarray, tuple[int, ...]]:
    """Resolve defaults, validate, and crop the search interval to the
    feasible region. Returns ``(y_2d, (lo, hi), beta_arr, orders)``.
    """
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    elif y.ndim != 2:
        raise ValueError(f"y must be 1-D or 2-D, got shape {y.shape}")
    N, M = y.shape

    if coarse_res <= 0:
        raise ValueError(f"coarse_res must be positive, got {coarse_res}")
    if beta_coarse_search is None:
        beta_arr = _DEFAULT_BETA_COARSE.copy()
    else:
        beta_arr = np.asarray(beta_coarse_search, dtype=float).ravel()
        if beta_arr.size == 0:
            raise ValueError("beta_coarse_search must be non-empty")
    orders = tuple(int(k) for k in orders)
    if not orders:
        raise ValueError("orders must be non-empty")
    if any(k < 1 for k in orders):
        raise ValueError(f"orders must be >= 1, got {orders}")

    if solver not in ("L-BFGS-B", "Nelder-Mead"):
        raise ValueError(
            f"solver must be 'L-BFGS-B' or 'Nelder-Mead', got {solver!r}"
        )
    if verbosity not in (0, 1, 2):
        raise ValueError(f"verbosity must be 0, 1, or 2, got {verbosity}")

    max_interval = N - max(orders) - 1
    if max_interval < 1:
        raise ValueError(
            f"signal too short: N={N} cannot accommodate orders={orders}"
        )

    if search_interval is None:
        lo, hi = default_interval(y)
    else:
        if len(search_interval) != 2:
            raise ValueError(f"search_interval must be (lo, hi), got {search_interval}")
        lo, hi = float(search_interval[0]), float(search_interval[1])

    hi = min(hi, max_interval)
    if lo > hi:
        lo, hi = 1.0, float(max_interval)

    if not (1.0 <= lo <= hi <= max_interval):
        raise ValueError(
            f"resolved search_interval=({lo}, {hi}) out of feasible "
            f"range [1, {max_interval}]"
        )

    return y, (lo, hi), beta_arr, orders


def estimate_bat(
    y,
    *,
    search_interval: tuple[float, float] | None = None,
    coarse_res: float = 0.25,
    beta_coarse_search=None,
    orders: Sequence[int] = _DEFAULT_ORDERS,
    common_bat: bool = False,
    solver: Literal["L-BFGS-B", "Nelder-Mead"] = "L-BFGS-B",
    verbosity: int = 1,
) -> EstimateResult:
    """Estimate bolus arrival time for one or more DCE-MRI signals.

    Spline-based estimator with parameters selected by generalized
    cross-validation. Mirrors ``DCEBE_estimateBAT.m`` from the MATLAB
    reference; returns BAT in 1-based sample units.

    Parameters
    ----------
    y : array_like, shape ``(N,)`` or ``(N, M)``
        One or several DCE signals as columns of length ``N``.
    search_interval : (lo, hi), optional
        Coarse-search bounds in 1-based sample units. If ``None``, use
        the per-signal envelope from
        :func:`dcebe._search.default_interval`.
    coarse_res : float, default 0.25
        Coarse-grid spacing in sample units.
    beta_coarse_search : array_like, optional
        Stiffness grid; ``alpha = beta ** (2 * k)``. Defaults to
        ``np.linspace(1, 25, 50)``.
    orders : sequence of int, default ``(3, 4, 5, 6)``
        Spline orders to evaluate.
    common_bat : bool, default False
        If True, estimate a single BAT shared by all signals.
    solver : {"L-BFGS-B", "Nelder-Mead"}, default "L-BFGS-B"
        Fine-search optimiser.
    verbosity : int, default 1
        ``0`` for silent, ``1`` for per-order progress, ``2`` for
        per-signal progress.

    Returns
    -------
    EstimateResult
        See :class:`EstimateResult` for field definitions.
    """
    y, (lo, hi), beta_arr, orders = _validate_inputs(
        y,
        search_interval,
        coarse_res,
        beta_coarse_search,
        orders,
        common_bat,
        solver,
        verbosity,
    )
    raise NotImplementedError(
        "estimate_bat main loop is in port-day 3; "
        "validated inputs would be: "
        f"shape={y.shape}, search=({lo}, {hi}), "
        f"|beta|={beta_arr.size}, orders={orders}"
    )
