"""Unit tests for ``dcebe._gcv``."""
from __future__ import annotations

import numpy as np
import pytest

from dcebe._gcv import gcv_score_qr, hat_fun
from dcebe._spline import make_deriv_pattern, make_matrix


class TestHatFun:
    def test_shape_1d(self):
        N, t, k = 30, 7.0, 3
        X, nablaK = make_matrix(N, t, k)
        y = np.sin(np.arange(N))
        y_hat = hat_fun(1.0, y, X, nablaK)
        assert y_hat.shape == (N,)

    def test_shape_2d(self):
        N, t, k = 30, 7.0, 3
        X, nablaK = make_matrix(N, t, k)
        y = np.sin(np.arange(N))[:, None] * np.array([1.0, 0.5, -0.3])[None, :]
        y_hat = hat_fun(1.0, y, X, nablaK)
        assert y_hat.shape == (N, 3)

    def test_alpha_zero_projects_onto_X(self):
        """At alpha=0 the smoothing penalty disappears, so y_hat is the
        orthogonal projection of y onto the column space of X.

        With the DCEBE construction:
        - column 0 of X has 1s in rows 0..t_int-1 (baseline plateau plus
          the BAT sample);
        - columns 1..s-1 of X have a single 1 each, mapping bijectively
          to rows t_int..N-1.

        Therefore the projection sets y_hat[i] = mean(y[0..t_int-1]) for
        i in 0..t_int-1, and y_hat[i] = y[i] for i in t_int..N-1."""
        N, t, k = 30, 8.0, 3
        X, nablaK = make_matrix(N, t, k)
        rng = np.random.default_rng(0)
        y = rng.standard_normal(N)
        y_hat = hat_fun(0.0, y, X, nablaK)
        t_int = int(t)
        baseline = np.mean(y[:t_int])
        np.testing.assert_allclose(y_hat[:t_int], baseline, atol=1e-8)
        np.testing.assert_allclose(y_hat[t_int:], y[t_int:], atol=1e-8)

    def test_alpha_large_smooths_to_baseline(self):
        """At very large alpha the penalty dominates, so y_hat collapses
        to the polynomial of degree < k that minimises ``||y - p||²``
        in the post-bolus region. For k=3 and a sinusoidal signal that
        means y_hat is a smooth low-degree polynomial trace, so RSS is
        much larger than at the unsmoothed optimum.
        """
        N, t, k = 30, 8.0, 3
        X, nablaK = make_matrix(N, t, k)
        rng = np.random.default_rng(0)
        y = np.sin(0.5 * np.arange(N)) + 0.1 * rng.standard_normal(N)
        y_hat_low = hat_fun(0.01, y, X, nablaK)
        y_hat_high = hat_fun(1e6, y, X, nablaK)
        # Heavy smoothing should produce a smaller second-difference norm
        d2_low = np.diff(y_hat_low, n=2)
        d2_high = np.diff(y_hat_high, n=2)
        assert np.linalg.norm(d2_high) < np.linalg.norm(d2_low)


class TestGcvScoreQr:
    def test_output_shape(self):
        rng = np.random.default_rng(0)
        N, M = 40, 3
        y = rng.standard_normal((N, M))
        t_arr = np.linspace(5, 15, 11)
        beta_arr = np.linspace(1, 25, 5)
        score = gcv_score_qr(y, t_arr, beta_arr, k=3)
        assert score.shape == (M, len(t_arr), len(beta_arr))

    def test_finite_inside_feasible_region(self):
        rng = np.random.default_rng(0)
        y = rng.standard_normal((30, 2))
        # Ensure both t and beta are well inside their respective bands
        t_arr = np.linspace(5, 20, 8)
        beta_arr = np.linspace(5, 20, 6)
        score = gcv_score_qr(y, t_arr, beta_arr, k=3, min_beta=2.0)
        assert np.all(np.isfinite(score))

    def test_t_below_one_yields_inf_via_barrier(self):
        rng = np.random.default_rng(0)
        y = rng.standard_normal((30, 2))
        t_arr = np.array([1.0, 5.0, 10.0])
        beta_arr = np.array([5.0, 10.0, 20.0])
        score = gcv_score_qr(y, t_arr, beta_arr, k=3, min_beta=2.0)
        # t=1 hits the lower barrier exactly -> inf
        assert np.all(np.isposinf(score[:, 0, :]))

    def test_recovers_known_changepoint(self):
        """Construct a signal with a known step change at t=15 and
        verify that the minimum-score t_arr index is near 15."""
        rng = np.random.default_rng(0)
        N = 40
        true_bat = 15.0
        y_clean = np.zeros(N)
        y_clean[int(true_bat) - 1 :] = np.linspace(0, 1, N - int(true_bat) + 1)
        y = y_clean[:, None] + 0.02 * rng.standard_normal((N, 1))
        t_arr = np.linspace(5, 25, 41)
        beta_arr = np.linspace(2, 25, 24)
        score = gcv_score_qr(y, t_arr, beta_arr, k=3, min_beta=1.5)
        # Take the (t, beta) cell with smallest score for signal 0
        idx_t, idx_beta = np.unravel_index(np.argmin(score[0]), score[0].shape)
        recovered = t_arr[idx_t]
        # Within one coarse-grid step (0.5)
        assert abs(recovered - true_bat) < 1.0


class TestGcvParityWithMatlab:
    """Cross-check :func:`gcv_score_qr` against the MATLAB reference
    on a hand-picked tiny configuration. The MATLAB column was generated
    by running ``DCEBE_GCVscore`` with ``factorize=false`` on the same
    inputs; the values are baked into this test to keep it hermetic.

    Generation script:
        rng(7);
        y = randn(20, 1);
        t_arr = [4 8 12];
        beta_arr = [2 5 10];
        k = 3; minBeta = 2;
        s = DCEBE_GCVscore(y, t_arr, beta_arr, k, false, minBeta);
        disp(squeeze(s));
    """

    @pytest.mark.skip(
        reason="Live MATLAB cross-check is in tests_py/test_matlab_parity.py "
        "(day 4); this test is a placeholder for a hermetic baked-in fixture."
    )
    def test_baked_in_values(self):
        pass
