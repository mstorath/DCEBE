"""Unit tests for ``dcebe._barriers``."""
from __future__ import annotations

import numpy as np
import pytest

from dcebe._barriers import barrier, broadcast_barriers


class TestBarrier:
    def test_at_lower_bound_is_inf(self):
        # x = b: sin(0) = 0, -log(0) = +inf
        assert np.isposinf(barrier(2.0, 2.0))

    def test_at_upper_bound_is_zero(self):
        # x = b + 1: sin(pi/2) = 1, -log(1) = 0
        assert barrier(3.0, 2.0) == 0.0

    def test_below_lower_is_inf(self):
        assert np.isposinf(barrier(1.0, 2.0))

    def test_above_upper_is_zero(self):
        assert barrier(5.0, 2.0) == 0.0

    def test_in_band_is_finite_positive(self):
        v = barrier(2.5, 2.0)
        assert 0.0 < v < np.inf

    def test_vector_input(self):
        x = np.array([1.0, 2.0, 2.5, 3.0, 5.0])
        out = barrier(x, 2.0)
        assert np.isposinf(out[0])
        assert np.isposinf(out[1])
        assert 0 < out[2] < np.inf
        assert out[3] == 0.0
        assert out[4] == 0.0

    def test_monotone_decreasing_in_band(self):
        """As x moves from b toward b+1, barrier monotonically decreases."""
        b = 0.0
        xs = np.linspace(b + 1e-3, b + 1.0 - 1e-3, 50)
        vals = barrier(xs, b)
        assert np.all(np.diff(vals) < 0)


class TestBroadcastBarriers:
    def test_output_shape(self):
        M, T, A = 3, 5, 7
        score = np.zeros((M, T, A))
        t_arr = np.linspace(2, 10, T)
        beta_arr = np.linspace(1, 25, A)
        out = broadcast_barriers(score, t_arr, beta_arr, N=20, k=3, min_beta=1.0)
        assert out.shape == (M, T, A)

    def test_beta_at_lower_bound_drives_inf(self):
        """When beta_arr starts below min_beta, the corresponding
        ``score[:, :, 0]`` slice is ``+inf``."""
        M, T, A = 2, 3, 5
        score = np.zeros((M, T, A))
        t_arr = np.full(T, 5.0)
        beta_arr = np.linspace(1.0, 25.0, A)
        out = broadcast_barriers(score, t_arr, beta_arr, N=20, k=3, min_beta=2.0)
        # min_beta - 1 = 1.0, so beta = 1.0 hits the barrier exactly
        assert np.all(np.isposinf(out[:, :, 0]))
        # beta = 25 is well above min_beta; barrier contributes 0
        assert np.all(np.isfinite(out[:, :, -1]))

    def test_t_at_lower_edge_drives_inf(self):
        """When t = 1 (lower bound for t_lo barrier), the slice is
        ``+inf``."""
        M, T, A = 2, 4, 3
        score = np.zeros((M, T, A))
        t_arr = np.array([1.0, 1.5, 5.0, 10.0])
        beta_arr = np.full(A, 25.0)
        out = broadcast_barriers(score, t_arr, beta_arr, N=20, k=3, min_beta=1.0)
        # t=1 hits the lower barrier
        assert np.all(np.isposinf(out[:, 0, :]))
        # t=10 is well inside the feasible region
        assert np.all(np.isfinite(out[:, 3, :]))

    def test_t_at_upper_edge_drives_inf(self):
        """When t = N - k (upper bound), the upper barrier is +inf."""
        M, T, A = 2, 3, 3
        score = np.zeros((M, T, A))
        N, k = 20, 3
        t_arr = np.array([5.0, 10.0, N - k])
        beta_arr = np.full(A, 25.0)
        out = broadcast_barriers(score, t_arr, beta_arr, N=N, k=k, min_beta=1.0)
        assert np.all(np.isposinf(out[:, 2, :]))
        assert np.all(np.isfinite(out[:, 0, :]))
