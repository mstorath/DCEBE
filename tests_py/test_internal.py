"""Unit tests for ``dcebe._spline``.

Three groups: ``fdweights`` correctness against closed-form polynomials,
``make_deriv_pattern`` consistency under ``h=1``, and ``make_matrix``
shape/structure properties at integer and fractional ``t``.
"""
from __future__ import annotations

from math import factorial

import numpy as np
import pytest

from dcebe._spline import fdweights, make_deriv_pattern, make_matrix


class TestFdweights:
    def test_first_derivative_centred(self):
        """First-derivative central FD on ``[0, 1, 2]`` at the middle
        node gives ``[-1/2, 0, 1/2]``."""
        w = fdweights(1.0, [0.0, 1.0, 2.0], 1)
        np.testing.assert_allclose(w, [-0.5, 0.0, 0.5], atol=1e-12)

    def test_zero_derivative_is_lagrange_interpolation(self):
        """0-th 'derivative' weights interpolate the function at the
        evaluation point. Weights sum to 1, and applying to sin(x)
        recovers sin(0.5) to within polynomial-fit accuracy."""
        x = np.array([0.0, 0.25, 1.0, 1.5, 2.5])
        w = fdweights(0.5, x, 0)
        np.testing.assert_allclose(np.sum(w), 1.0, atol=1e-12)
        np.testing.assert_allclose(w @ np.sin(x), np.sin(0.5), atol=1e-3)

    def test_second_derivative_polynomial(self):
        """Applying 2nd-derivative weights to ``f(x) = x^3`` at ``xi=1``
        yields ``f''(1) = 6``."""
        x = np.array([0.0, 1.0, 2.0, 3.0])
        w = fdweights(1.0, x, 2)
        np.testing.assert_allclose(w @ x**3, 6.0, atol=1e-10)

    def test_third_derivative_constant(self):
        """``f(x) = x^3`` has constant 3rd derivative ``= 6``."""
        x = np.array([0.0, 1.0, 2.0, 3.0])
        w = fdweights(0.0, x, 3)
        np.testing.assert_allclose(w @ x**3, 6.0, atol=1e-10)

    def test_first_derivative_one_sided_known_coefficients(self):
        """One-sided 4-point first-derivative weights at 0 with nodes
        ``[0, h, 2h, 3h]`` are ``[-11/6, 3, -3/2, 1/3] / h``."""
        h = 0.4
        x = h * np.arange(4)
        w = fdweights(0.0, x, 1)
        expected = np.array([-11.0 / 6, 3.0, -1.5, 1.0 / 3]) / h
        np.testing.assert_allclose(w, expected, atol=1e-10)

    def test_polynomial_exactness_up_to_degree_p(self):
        """Weights are exact for polynomials of degree at most ``p =
        len(x)-1``. With 5 nodes, exact through degree 4."""
        x = np.linspace(-1.0, 2.0, 5)
        rng = np.random.default_rng(0)
        coeffs = rng.standard_normal(5)
        f = np.polynomial.polynomial.polyval(x, coeffs)
        for m in range(5):
            w = fdweights(0.5, x, m)
            true_deriv = np.polynomial.polynomial.polyval(
                0.5, np.polynomial.polynomial.polyder(coeffs, m)
            )
            np.testing.assert_allclose(w @ f, true_deriv, atol=1e-9)


class TestMakeDerivPattern:
    @pytest.mark.parametrize("k", [3, 4, 5, 6])
    def test_default_pattern_recovers_kth_derivative_of_xk(self, k):
        """Integer-grid k-th-derivative weights applied to ``f(x)=x^k``
        evaluated at ``x=0`` give ``k!``."""
        p = make_deriv_pattern(k)
        x = np.arange(k + 1, dtype=float)
        np.testing.assert_allclose(p @ x**k, factorial(k), atol=1e-9)

    @pytest.mark.parametrize("k", [3, 4, 5, 6])
    def test_h_eq_1_matches_default(self, k):
        """``make_deriv_pattern(k, l=k+1, h=1.0) == make_deriv_pattern(k)``."""
        p_default = make_deriv_pattern(k)
        p_h1 = make_deriv_pattern(k, l=k + 1, h=1.0)
        np.testing.assert_allclose(p_default, p_h1, atol=1e-12)

    def test_h_continuity(self):
        """Pattern is continuous in ``h``: a tiny perturbation to ``h``
        produces a small change in weights (no discontinuity)."""
        p1 = make_deriv_pattern(3, h=1.0)
        p2 = make_deriv_pattern(3, h=0.99)
        np.testing.assert_allclose(p1, p2, atol=0.5)


class TestMakeMatrix:
    def test_shapes_at_integer_t(self):
        N, t, k = 50, 10.0, 3
        X, nabla = make_matrix(N, t, k)
        s = N - int(t) + 1
        assert X.shape == (N, s)
        assert nabla.shape == (s - k, s)

    def test_X_baseline_then_identity(self):
        """``X`` has constant first column for the first ``t_int-1``
        rows and an identity over the remaining rows."""
        N, t, k = 20, 5.0, 3
        X, _ = make_matrix(N, t, k)
        Xd = X.toarray()
        for i in range(int(t) - 1):
            assert Xd[i, 0] == 1.0
            assert np.all(Xd[i, 1:] == 0.0)
        np.testing.assert_array_equal(Xd[int(t) - 1 :], np.eye(N - int(t) + 1))

    def test_first_row_collapses_to_base_at_integer_t(self):
        """At integer ``t`` (``h=1``), ``nablaK[0, :k+1]`` equals the
        integer-grid pattern and ``sqrt(h)=1`` is a no-op."""
        N, t, k = 30, 7.0, 3
        _, nabla = make_matrix(N, t, k)
        first_row = nabla.toarray()[0, : k + 1]
        np.testing.assert_allclose(first_row, make_deriv_pattern(k), atol=1e-12)

    def test_first_row_at_half_integer_t(self):
        """At ``t = t_int + 0.5`` the first row equals
        ``sqrt(0.5) * make_deriv_pattern(k, l, h=0.5)``."""
        N, t, k = 30, 7.5, 3
        _, nabla = make_matrix(N, t, k)
        first_row = nabla.toarray()[0, : k + 1]
        h = 0.5
        expected = make_deriv_pattern(k, l=k + 1, h=h) * np.sqrt(h)
        np.testing.assert_allclose(first_row, expected, atol=1e-12)

    def test_band_below_first_row_is_unchanged(self):
        """Rows ``1..nrows-1`` of ``nablaK`` carry the unchanged
        integer-grid k-th-derivative band."""
        N, t, k = 30, 7.5, 3
        _, nabla = make_matrix(N, t, k)
        nablad = nabla.toarray()
        base = make_deriv_pattern(k)
        for r in range(1, nablad.shape[0]):
            np.testing.assert_allclose(nablad[r, r : r + k + 1], base, atol=1e-12)

    @pytest.mark.parametrize("bad_t", [0.5, 0.0, -1.0])
    def test_t_below_one_raises(self, bad_t):
        with pytest.raises(ValueError, match="out of bounds"):
            make_matrix(20, bad_t, 3)

    def test_t_above_max_raises(self):
        # max is N - (k+1) = 20 - 4 = 16
        with pytest.raises(ValueError, match="out of bounds"):
            make_matrix(20, 17.0, 3)

    def test_t_at_lower_bound(self):
        """``t = 1`` is allowed; the upper block ``XU`` is empty so
        ``X = eye(N)``."""
        N, t, k = 10, 1.0, 3
        X, _ = make_matrix(N, t, k)
        assert X.shape == (N, N)
        np.testing.assert_array_equal(X.toarray(), np.eye(N))

    def test_t_at_upper_bound(self):
        """``t = N - l`` is allowed (boundary-inclusive)."""
        N, k = 10, 3
        l = k + 1
        t = N - l
        X, nabla = make_matrix(N, t, k)
        s = N - t + 1
        assert X.shape == (N, s)
        assert nabla.shape == (s - k, s)

    @pytest.mark.parametrize("k", [3, 4, 5, 6])
    def test_kernel_of_nabla_contains_polynomials_below_degree_k(self, k):
        """The k-th derivative of any polynomial of degree < k is zero,
        so ``nablaK @ p(x_grid) ~ 0`` for such ``p``. We check this on
        the homogeneous-band rows (skipping the first, which has the
        sqrt(h) scaling and the fractional pattern)."""
        N, t = 40, 5.0  # integer t so first row is also integer-grid
        X, nabla = make_matrix(N, t, k)
        s = X.shape[1]
        x_grid = np.arange(s, dtype=float)
        for d in range(k):
            poly = x_grid**d
            residual = nabla @ poly
            np.testing.assert_allclose(residual, 0.0, atol=1e-8)
