"""Top-level driver for bolus arrival time estimation.

Orchestrates the per-order coarse-then-fine search documented in
``DCEBE_estimateBAT.m``: for each spline order ``k``, evaluate the
log-GCV score on a coarse ``(t, beta)`` grid, then run a local optimiser
from the best-coarse cell to refine the parameters. The best-over-orders
``(BAT, beta, k)`` is kept per signal (per-signal mode) or jointly
(common-BAT mode).

BAT values returned in :class:`EstimateResult` are 1-based sample units,
matching the MATLAB reference.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from scipy.optimize import minimize

from ._gcv import gcv_score_qr, hat_fun
from ._search import default_interval
from ._spline import make_deriv_pattern, make_matrix

_DEFAULT_BETA_COARSE = np.linspace(1.0, 25.0, 50)
_DEFAULT_ORDERS = (3, 4, 5, 6)
_FMINUNC_OPTIONS = {"maxiter": 1000, "gtol": 1e-6}
_FMINSEARCH_OPTIONS = {"maxiter": 100000, "maxfev": 100000, "xatol": 1e-4, "fatol": 1e-4}


@dataclass
class EstimateResult:
    """Container for :func:`estimate_bat` output.

    Attributes
    ----------
    bat : ndarray, shape ``(M,)``
        Estimated bolus arrival time per signal, in 1-based sample units
        (matching the MATLAB reference). Convert to seconds via
        ``(bat - 1) * delta_t``.
    cp_int : ndarray, shape ``(M,)`` of int
        ``floor(bat)``: the integer changepoint index (1-based).
    beta_opt : ndarray, shape ``(M,)``
        Optimal stiffness parameter per signal.
    k_opt : ndarray, shape ``(M,)`` of int
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
    y,
    search_interval,
    coarse_res,
    beta_coarse_search,
    orders,
    common_bat,
    solver,
    verbosity,
):
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

    return y, (lo, hi), beta_arr, orders


def _fine_search(f, x0, solver_name):
    """Mirror of MATLAB's ``fminunc`` ('quasi-newton', central FD) /
    ``fminsearch`` (Nelder-Mead) options. Returns ``(x_opt, f_opt)``."""
    if solver_name == "L-BFGS-B":
        result = minimize(f, x0, method="L-BFGS-B", options=_FMINUNC_OPTIONS)
    else:  # Nelder-Mead
        result = minimize(f, x0, method="Nelder-Mead", options=_FMINSEARCH_OPTIONS)
    return result.x, float(result.fun)


def estimate_bat(
    y,
    *,
    search_interval=None,
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
        Fine-search optimiser. ``L-BFGS-B`` mirrors MATLAB's ``fminunc``
        with the ``quasi-newton`` algorithm; ``Nelder-Mead`` mirrors
        ``fminsearch``. Note: SciPy uses forward finite differences by
        default whereas MATLAB's ``fminunc`` defaults to central — the
        fine-search optima may differ at the 1e-3 level.
    verbosity : int, default 1
        ``0`` for silent, ``1`` for per-order progress, ``2`` for
        per-signal progress.
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
    N, M = y.shape
    coarse_t = np.arange(lo, hi + 0.5 * coarse_res, coarse_res)
    T = len(coarse_t)
    A = len(beta_arr)
    min_beta = float(beta_arr[0])

    score_opt = np.full(M, np.inf)
    bat = np.full(M, np.nan)
    beta_opt = np.full(M, np.nan)
    k_opt = np.full(M, -1, dtype=int)

    L = len(orders)
    for li, k in enumerate(orders):
        if verbosity > 0:
            print(f"Order {k} ({li + 1} of {L}):\n  Coarse search...")

        score_coarse = np.log(
            gcv_score_qr(y, coarse_t, beta_arr, k, min_beta=min_beta)
        )

        if verbosity > 0:
            print("  Fine search...")

        if common_bat:
            score_summary = np.mean(score_coarse, axis=0)  # (T, A)
            min_idx_flat = int(np.argmin(score_summary))
            idx_t, idx_b = np.unravel_index(min_idx_flat, (T, A))
            x0 = np.array([coarse_t[idx_t], beta_arr[idx_b]])

            def f(par, _k=k):
                t_, b_ = float(par[0]), float(par[1])
                s = gcv_score_qr(y, [t_], [b_], _k, min_beta=min_beta)
                return float(np.mean(np.log(s)))

            x_opt, f_opt = _fine_search(f, x0, solver)
            if f_opt > f(x0):
                warnings.warn("Fine search not successful (common_bat)", UserWarning)

            if f_opt < score_opt[0]:
                score_opt[:] = f_opt
                bat[:] = x_opt[0]
                beta_opt[:] = x_opt[1]
                k_opt[:] = k

        else:
            for m in range(M):
                min_idx_flat = int(np.argmin(score_coarse[m]))
                idx_t, idx_b = np.unravel_index(min_idx_flat, (T, A))
                x0 = np.array([coarse_t[idx_t], beta_arr[idx_b]])
                y_m = y[:, m : m + 1]

                def f(par, _k=k, _y=y_m):
                    t_, b_ = float(par[0]), float(par[1])
                    s = gcv_score_qr(_y, [t_], [b_], _k, min_beta=min_beta)
                    return float(np.log(s[0, 0, 0]))

                x_opt, f_opt = _fine_search(f, x0, solver)
                if f_opt > f(x0):
                    warnings.warn(
                        f"Fine search not successful (signal {m})", UserWarning
                    )

                if f_opt < score_opt[m]:
                    score_opt[m] = f_opt
                    bat[m] = x_opt[0]
                    beta_opt[m] = x_opt[1]
                    k_opt[m] = k

                if verbosity > 1 or (verbosity > 0 and (m + 1) % 50 == 0):
                    print(
                        f"  {m + 1} of {M}, LogGCV: {score_opt[m]:.6f}, "
                        f"beta: {beta_opt[m]:.6f}, k_opt {k_opt[m]}, "
                        f"CP: {bat[m]:.3f}"
                    )

    cp_int = np.floor(bat).astype(int)
    y_hat = np.zeros((N, M))
    y_hat_x = np.zeros((N, M))
    for m in range(M):
        k_m = int(k_opt[m])
        deriv = make_deriv_pattern(k_m)
        X, nablaK = make_matrix(N, float(bat[m]), k_m, deriv)
        y_hat[:, m] = hat_fun(beta_opt[m] ** (2 * k_m), y[:, m], X, nablaK)
        y_hat_x[:, m] = np.arange(1, N + 1)
        # cp_int is 1-based; replace that row with the fractional bat
        y_hat_x[cp_int[m] - 1, m] = bat[m]

    return EstimateResult(
        bat=bat,
        cp_int=cp_int,
        beta_opt=beta_opt,
        k_opt=k_opt,
        score_opt=score_opt,
        y_hat=y_hat,
        y_hat_x=y_hat_x,
    )
