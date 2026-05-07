"""Spline matrix construction and finite-difference weights.

Translation of ``DCEBE_make_matrix.m``, ``DCEBE_make_deriv_pattern.m`` and
``external/fdweights/fdweights.m`` to Python. Given a signal of length
``N`` and a fractional bolus arrival time ``t in [1, N-l]``, build the
data fidelity matrix ``X`` (a constant baseline plateau before ``t``,
identity after) and the k-th order derivative operator ``nablaK`` with
a special boundary row that accounts for ``t`` falling between integer
samples.

BAT indices in this module follow the MATLAB 1-based convention so that
parity tests against the MATLAB reference are direct comparisons.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse


def fdweights(xi: float, x, m: int) -> np.ndarray:
    """Finite-difference weights at the evaluation point ``xi`` using
    nodes ``x`` for the m-th derivative.

    Re-implementation of the algorithm by Toby Driscoll (Copyright 2008,
    BSD-2-Clause; see ``external/fdweights/license.txt`` in this
    repository). Reference: Fornberg, *A Practical Guide to Pseudospectral
    Methods*, Cambridge University Press.

    Parameters
    ----------
    xi : float
        Point at which to evaluate the derivative.
    x : array_like, 1-D
        Nodes (need not be uniformly spaced).
    m : int, non-negative
        Order of the derivative (0 for interpolation).

    Returns
    -------
    w : ndarray, shape ``(len(x),)``
        Weights such that ``w @ f(x) ~ f^(m)(xi)`` exactly for
        polynomials of degree at most ``len(x) - 1``.
    """
    x_shifted = np.asarray(x, dtype=float).ravel() - xi
    p = len(x_shifted) - 1
    memo: dict[tuple[int, int, int], float] = {}

    def weight(m_: int, j: int, k_: int) -> float:
        key = (m_, j, k_)
        cached = memo.get(key)
        if cached is not None:
            return cached
        if m_ < 0 or m_ > j:
            r = 0.0
        elif m_ == 0 and j == 0:
            r = 1.0
        elif k_ < j:
            r = (
                x_shifted[j] * weight(m_, j - 1, k_)
                - m_ * weight(m_ - 1, j - 1, k_)
            ) / (x_shifted[j] - x_shifted[k_])
        else:  # k_ == j
            num = float(np.prod(x_shifted[j - 1] - x_shifted[: j - 1])) if j >= 1 else 1.0
            den = float(np.prod(x_shifted[j] - x_shifted[:j])) if j >= 1 else 1.0
            beta_ = num / den
            r = beta_ * (
                m_ * weight(m_ - 1, j - 1, j - 1)
                - x_shifted[j - 1] * weight(m_, j - 1, j - 1)
            )
        memo[key] = r
        return r

    return np.array([weight(m, p, k_) for k_ in range(p + 1)])


def make_deriv_pattern(k: int, l: int | None = None, h: float = 1.0) -> np.ndarray:
    """k-th derivative weights at 0 with nodes ``[0, h, h+1, ..., h+l-2]``.

    With ``h = 1`` (default) this gives the standard k-th derivative
    weights on the integer grid ``[0, 1, ..., l-1]``. The free ``h``
    parameter is what allows DCEBE to support a fractional changepoint:
    the first node sits at 0 and the remaining ``l-1`` nodes start at
    ``h`` and continue at unit spacing.
    """
    if l is None:
        l = k + 1
    nodes = np.concatenate([[0.0], h + np.arange(l - 1, dtype=float)])
    return fdweights(0.0, nodes, k)


def make_matrix(
    N: int,
    t: float,
    k: int,
    deriv_pattern: np.ndarray | None = None,
):
    """Spline matrices ``X`` and ``nablaK`` for the DCEBE model.

    Builds the data fidelity matrix ``X`` (shape ``N x s``) and the k-th
    order derivative operator ``nablaK`` (shape ``(s-k) x s``) for a
    signal of length ``N`` with bolus arrival at fractional sample ``t``
    (1-based, matching the MATLAB reference).

    Parameters
    ----------
    N : int
        Signal length.
    t : float
        Bolus arrival time in 1-based sample units. Must satisfy
        ``1 <= t <= N - (k+1)``.
    k : int
        Spline order.
    deriv_pattern : ndarray, optional
        Pre-computed integer-grid k-th-derivative weights. Computed via
        :func:`make_deriv_pattern` if not supplied.

    Returns
    -------
    X : scipy.sparse.csr_matrix, shape ``(N, s)``
        Rows ``0..t_int-2`` are constant baseline (only the first column
        is 1); rows ``t_int-1..N-1`` form an identity over the post-bolus
        samples. Here ``s = N - t_int + 1`` and ``t_int = floor(t)``.
    nablaK : scipy.sparse.csr_matrix, shape ``(s-k, s)``
        k-th-order derivative operator. The first row is replaced by
        ``sqrt(h) * fd_weights_at_0_with_h``, where ``h = 1 - (t - t_int)``,
        accounting for the fractional bolus position. Remaining rows form
        a standard banded difference operator.
    """
    l = k + 1
    if t < 1 or t > N - l:
        raise ValueError(f"t={t} out of bounds [1, {N - l}]")
    if deriv_pattern is None:
        deriv_pattern = make_deriv_pattern(k)
    t_int = int(np.floor(t))
    t_frac = t - t_int
    h = 1.0 - t_frac
    s = N - t_int + 1

    if t_int > 1:
        rows = np.arange(t_int - 1)
        cols = np.zeros(t_int - 1, dtype=int)
        XU = sparse.coo_matrix(
            (np.ones(t_int - 1), (rows, cols)), shape=(t_int - 1, s)
        )
        X = sparse.vstack([XU, sparse.eye(s, format="csr")], format="csr")
    else:
        X = sparse.eye(s, format="csr")

    nrows = s - l + 1
    diagonals = [deriv_pattern[j] * np.ones(nrows) for j in range(l)]
    offsets = list(range(l))
    nablaK = sparse.diags(diagonals, offsets, shape=(nrows, s), format="lil")

    frac_pattern = make_deriv_pattern(k, l=l, h=h) * np.sqrt(h)
    nablaK[0, :l] = frac_pattern

    return X.tocsr(), nablaK.tocsr()
