"""Generalized cross-validation score and stable hat-matrix evaluation.

QR-based path of ``DCEBE_GCVscore.m``. The factorisation-based path
(``DCEBE_factorize.m``) is in :mod:`dcebe._factorize` and is opt-in via
``factorize=True``; the QR path is the default and is the only one
exercised by the parity tests, since it is consistently better
conditioned (per the comment block in the MATLAB source).

The hat matrix construction
::

    H = X * ([X; sqrt(alpha) * nablaK] \\ [I_N; 0])

solves a tall least-squares system for the columns of ``H``. We use
``numpy.linalg.lstsq`` on the dense augmented matrix; ``s`` and ``N``
are typically below 100 so the dense solve is fast and conditioning is
bounded.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

from ._barriers import broadcast_barriers
from ._spline import make_deriv_pattern, make_matrix

_H_TOL = 1.0e-6


def hat_fun(
    alpha: float,
    y: np.ndarray,
    X: sparse.spmatrix,
    nablaK: sparse.spmatrix,
) -> np.ndarray:
    """Spline reconstruction ``y_hat = H y`` via the stable QR-based
    formulation that mirrors ``DCEBE_hat_fun.m``.

    Parameters
    ----------
    alpha : float
        Smoothing parameter ``alpha = beta ** (2 * k)``.
    y : ndarray, shape ``(N,)`` or ``(N, M)``
        Signal(s).
    X, nablaK : scipy.sparse.spmatrix
        From :func:`dcebe._spline.make_matrix`.
    """
    y = np.asarray(y, dtype=float)
    one_d = y.ndim == 1
    if one_d:
        y = y[:, None]
    A_aug = sparse.vstack([X, np.sqrt(alpha) * nablaK]).toarray()
    rhs = np.vstack([y, np.zeros((nablaK.shape[0], y.shape[1]))])
    coeffs, *_ = np.linalg.lstsq(A_aug, rhs, rcond=None)
    y_hat = X.toarray() @ coeffs
    return y_hat[:, 0] if one_d else y_hat


def gcv_score_qr(
    y: np.ndarray,
    t_arr: np.ndarray,
    beta_arr: np.ndarray,
    k: int,
    *,
    min_beta: float | None = None,
) -> np.ndarray:
    """Generalized cross-validation score over a ``(t, beta)`` grid.

    Returns a ``(M, T, A)`` array of scores including the three barrier
    contributions. The score formula is

    .. math::

        \\text{score}_{m,i,j} = \\frac{1}{N} \\frac{\\|y_m - H y_m\\|^2}
                                            {(1 - \\operatorname{tr}(H) / N)^2}

    where ``H`` depends on ``(t_i, beta_j, k)``.

    Parameters
    ----------
    y : ndarray, shape ``(N,)`` or ``(N, M)``
        DCE signals as columns.
    t_arr : array_like, shape ``(T,)``
        BAT candidates (1-based).
    beta_arr : array_like, shape ``(A,)``
        Stiffness candidates; ``alpha = beta ** (2 * k)``.
    k : int
        Spline order.
    min_beta : float, optional
        Lower bound for the beta barrier. Defaults to ``beta_arr[0]``.
    """
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    N, M = y.shape
    t_arr = np.atleast_1d(np.asarray(t_arr, dtype=float))
    beta_arr = np.atleast_1d(np.asarray(beta_arr, dtype=float))
    T, A = len(t_arr), len(beta_arr)
    alpha = beta_arr ** (2 * k)
    if min_beta is None:
        min_beta = float(beta_arr[0])

    deriv_pattern = make_deriv_pattern(k)
    score = np.zeros((M, T, A))

    for i in range(T):
        t = float(t_arr[i])
        # Mirror DCEBE_GCVscore.m lines 22-27: keep h floored at tol so
        # the fractional row of nablaK does not collapse to zero scale.
        t_frac = t - np.floor(t)
        h = 1.0 - t_frac
        if abs(h) < _H_TOL:
            t = float(np.ceil(t)) - _H_TOL
        # MATLAB clamps to [1, N - k - 1] before constructing the matrix;
        # this lets the barriers (added later) be the only thing that
        # penalises the score outside the feasible region.
        t_inbound = max(min(t, N - k - 1), 1.0)
        X, nablaK = make_matrix(N, t_inbound, k, deriv_pattern)
        Xd = X.toarray()
        nablaKd = nablaK.toarray()
        nrows_n = nablaKd.shape[0]

        for j in range(A):
            A_aug = np.vstack([Xd, np.sqrt(alpha[j]) * nablaKd])
            n_aug = A_aug.shape[0]
            Z = np.eye(n_aug, N)
            coeffs, *_ = np.linalg.lstsq(A_aug, Z, rcond=None)
            H = Xd @ coeffs
            df = float(np.trace(H))
            y_hat = H @ y
            rss = np.sum((y_hat - y) ** 2, axis=0)
            denom = (1.0 - df / N) ** 2
            # If denom <= 0 (hat matrix near full rank), the formula
            # blows up; barriers handle infeasible regions, so we
            # accept the +inf here.
            with np.errstate(divide="ignore", invalid="ignore"):
                score[:, i, j] = rss / (N * denom)

    score = broadcast_barriers(score, t_arr, beta_arr, N=N, k=k, min_beta=min_beta)
    return score
