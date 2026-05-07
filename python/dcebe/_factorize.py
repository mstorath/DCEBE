"""Matrix factorisation for the SVD/eigen-based GCV path.

Translation of ``DCEBE_factorize.m`` and the ``factorize=true`` branch
of ``DCEBE_GCVscore.m``. This path computes the factorisation once per
``(t, k)`` and then evaluates the GCV score for many ``alpha`` values
cheaply by reusing the factorised representation.

The MATLAB source flags this path as less well-conditioned than the
QR-based variant; we ship it for completeness and future profiling but
the user-facing :func:`dcebe.estimate_bat` calls only the QR path.

Reference: Kent and Mohammadzadeh, "Global optimization of the
generalized cross-validation criterion", 2000.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import sqrtm

from ._barriers import broadcast_barriers
from ._spline import make_deriv_pattern, make_matrix


def factorize(X: np.ndarray, C: np.ndarray, k: int):
    """Three-step matrix factorisation enabling fast multi-alpha GCV.

    Returns ``(U, d, p, q, r)`` where ``U`` is an orthogonal ``N x N``
    matrix, ``d`` is a length-``q`` vector of positive eigenvalues, and
    ``(p, q, r)`` are the dimensions used in the GCV score formula.
    """
    X = np.asarray(X, dtype=float)
    C = np.asarray(C, dtype=float)
    p = C.shape[0]
    N = X.shape[0]

    # Step 1: eigendecomposition of C, descending order
    delta_full, Gamma_full = np.linalg.eigh(C)
    Gamma = Gamma_full[:, ::-1]
    delta = delta_full[::-1]

    sqrt_eps = float(np.sqrt(np.finfo(float).eps))
    r = min(int(np.sum(delta > sqrt_eps)), p - k)
    delta_1 = np.maximum(0.0, delta[:r])

    scale = np.concatenate([delta_1 ** (-0.5), np.ones(p - r)])
    P_1 = Gamma * scale[None, :]
    X_2 = X @ P_1
    X_2_1 = X_2[:, :r]
    X_2_2 = X_2[:, r:]

    # Step 2
    B = np.linalg.inv(X_2_2.T @ X_2_2)
    ReBSqrt = np.real(sqrtm(B))
    P_2 = np.block(
        [
            [np.eye(r), np.zeros((r, p - r))],
            [-B @ X_2_2.T @ X_2_1, ReBSqrt],
        ]
    )
    X_3 = X_2 @ P_2

    # Step 3: SVD of X_3
    U, S_diag, _ = np.linalg.svd(X_3, full_matrices=True)

    eps_ = float(np.finfo(float).eps)
    q = max(1, min(int(np.sum(S_diag > eps_)), r))

    # Reorder so the (p - q) singular values closest to 1 come last among
    # the first p columns. The MATLAB source assumes those indices form a
    # contiguous block at the bottom of the spectrum; we use a robust
    # set-difference instead so non-contiguous configurations also work.
    sort_order = np.argsort(np.abs(S_diag[:p] - 1.0))
    idx_ones = sort_order[: p - q]
    not_ones = [i for i in range(p) if i not in set(idx_ones.tolist())]
    permutation = not_ones + list(idx_ones)
    permutation_full = permutation + list(range(p, N))

    U = U[:, permutation_full]
    # MATLAB: S = S(permutation2, permutation), S_diag = diag(S)
    # The relevant entries for d are the first q diagonal entries of the
    # permuted singular-value matrix, which are the singular values at
    # positions not_ones[:q] in the original spectrum.
    S_perm = S_diag[permutation[:q]]
    d = S_perm ** 2

    return U, d, p, q, r


def gcv_score_factorize(
    y: np.ndarray,
    t_arr,
    beta_arr,
    k: int,
    *,
    min_beta: float | None = None,
) -> np.ndarray:
    """``gcv_score`` via the factorisation path. Same interface as
    :func:`dcebe._gcv.gcv_score_qr`."""
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
    tol = 1e-6

    for i in range(T):
        t = float(t_arr[i])
        t_frac = t - np.floor(t)
        h = 1.0 - t_frac
        if abs(h) < tol:
            t = float(np.ceil(t)) - tol
        t_inbound = max(min(t, N - k - 1), 1.0)
        X, nablaK = make_matrix(N, t_inbound, k, deriv_pattern)

        Xd = X.toarray()
        nablaKd = nablaK.toarray()
        NTN = nablaKd.T @ nablaKd

        U, d, p, q, r = factorize(Xd, NTN, k)
        yStar = U.T @ y
        z = yStar[:q, :]                 # (q, M)
        w = yStar[p:, :]                 # (N - p, M); yStar2 in MATLAB is empty
        w_sum = np.sum(w ** 2, axis=0)   # (M,)

        npqr = N - p - q + r

        v = np.sum(alpha[None, :] / (d[:, None] + alpha[None, :]), axis=0) + npqr

        for m in range(M):
            u_mat = (
                (z[:, m, None] ** 2)
                * (alpha[None, :] ** 2)
                / ((d[:, None] + alpha[None, :]) ** 2)
            )
            u_per_alpha = np.sum(u_mat, axis=0) + w_sum[m]
            score[m, i, :] = N * u_per_alpha / v ** 2

    return broadcast_barriers(score, t_arr, beta_arr, N=N, k=k, min_beta=min_beta)
