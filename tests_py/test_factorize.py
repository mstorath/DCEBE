"""Tests for ``dcebe._factorize`` (the SVD/eigen GCV path).

The factorise path is opt-in and not exposed via :func:`estimate_bat`;
these tests validate it as a private utility and check internal
consistency with the QR path on small inputs.
"""
from __future__ import annotations

import numpy as np
import pytest

from dcebe._factorize import factorize, gcv_score_factorize
from dcebe._gcv import gcv_score_qr
from dcebe._spline import make_matrix


class TestFactorize:
    def test_shapes(self):
        N, t, k = 30, 7.0, 3
        X, nabla = make_matrix(N, t, k)
        Xd = X.toarray()
        NTN = nabla.toarray().T @ nabla.toarray()
        U, d, p, q, r = factorize(Xd, NTN, k)
        s = Xd.shape[1]
        assert p == s
        assert U.shape == (N, N)
        assert d.shape == (q,)
        assert 1 <= q <= r <= p
        assert np.all(d >= 0)

    def test_U_is_orthogonal(self):
        N, t, k = 30, 7.0, 3
        X, nabla = make_matrix(N, t, k)
        NTN = nabla.toarray().T @ nabla.toarray()
        U, *_ = factorize(X.toarray(), NTN, k)
        np.testing.assert_allclose(U.T @ U, np.eye(N), atol=1e-10)


class TestGcvScoreFactorizeMatchesQr:
    """Internal consistency: both GCV paths compute the same score
    formula, so they should agree to a comfortable tolerance."""

    @pytest.mark.parametrize("k", [3, 4, 5])
    def test_agrees_with_qr_on_synthetic_signal(self, k):
        rng = np.random.default_rng(0)
        N, M = 30, 2
        # Use a smooth-then-step signal so the GCV is well-behaved
        y = rng.standard_normal((N, M)) * 0.1
        y[15:] += np.linspace(0, 1, N - 15)[:, None]
        t_arr = np.array([8.0, 10.0, 12.0])
        beta_arr = np.array([5.0, 10.0, 20.0])
        s_qr = gcv_score_qr(y, t_arr, beta_arr, k, min_beta=1.0)
        s_fc = gcv_score_factorize(y, t_arr, beta_arr, k, min_beta=1.0)
        # The two paths use different numerical decompositions; the
        # MATLAB source comments that the factorise path is less
        # well-conditioned, so we use a generous relative tolerance.
        np.testing.assert_allclose(s_qr, s_fc, rtol=1e-4, atol=1e-4)
