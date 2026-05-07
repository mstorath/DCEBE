"""End-to-end tests for ``dcebe.estimate_bat``.

Synthetic-recovery: build piecewise signals with known BATs, estimate,
and check we recover them within tolerance. These are not parity tests;
that role belongs to ``test_matlab_parity.py`` (day 4).
"""
from __future__ import annotations

import numpy as np
import pytest

from dcebe import EstimateResult, estimate_bat


def _bolus_signal(N: int, bat: float, peak: float = 1.0, slope: float = 0.05):
    """A piecewise signal: zero baseline, linear ramp from sample
    floor(bat) onwards, capped at ``peak``."""
    out = np.zeros(N)
    bat_int = int(np.floor(bat))
    n_post = N - bat_int + 1
    ramp_len = min(n_post, int(peak / slope) + 1)
    out[bat_int - 1 : bat_int - 1 + ramp_len] = np.linspace(0, peak, ramp_len)
    out[bat_int - 1 + ramp_len :] = peak
    return out


@pytest.fixture
def rng():
    return np.random.default_rng(42)


class TestEstimateBatSyntheticRecovery:
    def test_single_signal_clean(self, rng):
        """No noise: BAT should be recovered to within ~1 sample."""
        N = 40
        true_bat = 15.0
        y = _bolus_signal(N, true_bat)
        result = estimate_bat(y, search_interval=(8, 25), verbosity=0)
        assert isinstance(result, EstimateResult)
        assert result.bat.shape == (1,)
        assert abs(result.bat[0] - true_bat) < 1.0

    def test_single_signal_noisy(self, rng):
        """With moderate noise, BAT recovered within 2 samples."""
        N = 50
        true_bat = 20.0
        y = _bolus_signal(N, true_bat) + 0.03 * rng.standard_normal(N)
        result = estimate_bat(y, search_interval=(10, 35), verbosity=0)
        assert abs(result.bat[0] - true_bat) < 2.0

    def test_multi_signal_independent(self, rng):
        """Multiple signals each with their own BAT."""
        N = 50
        true_bats = np.array([15.0, 20.0, 25.0])
        y = np.column_stack([_bolus_signal(N, b) for b in true_bats])
        y += 0.02 * rng.standard_normal(y.shape)
        result = estimate_bat(y, search_interval=(8, 35), verbosity=0)
        assert result.bat.shape == (3,)
        for est, true in zip(result.bat, true_bats):
            assert abs(est - true) < 2.0

    def test_common_bat(self, rng):
        """common_bat=True returns a single shared BAT for all signals."""
        N = 50
        true_bat = 20.0
        y = np.column_stack([_bolus_signal(N, true_bat) for _ in range(5)])
        y += 0.05 * rng.standard_normal(y.shape)
        result = estimate_bat(y, search_interval=(10, 35), common_bat=True, verbosity=0)
        assert np.all(result.bat == result.bat[0])
        assert abs(result.bat[0] - true_bat) < 1.5

    def test_default_search_interval(self, rng):
        """Without an explicit search_interval, the heuristic should
        still produce a reasonable estimate."""
        N = 50
        true_bat = 18.0
        y = _bolus_signal(N, true_bat) + 0.02 * rng.standard_normal(N)
        result = estimate_bat(y, verbosity=0)
        assert abs(result.bat[0] - true_bat) < 3.0

    def test_solver_nelder_mead_runs_on_clean_signal(self):
        """Nelder-Mead path runs without crashing and on a noiseless
        bolus matches L-BFGS-B to within the coarse-grid resolution.
        On noisy objectives the two solvers can diverge (Nelder-Mead is
        known to be fragile on noisy log-GCV landscapes); we test only
        the clean case here."""
        N = 40
        true_bat = 15.0
        y = _bolus_signal(N, true_bat)
        r_lbfgs = estimate_bat(
            y, search_interval=(8, 25), solver="L-BFGS-B", verbosity=0
        )
        r_nm = estimate_bat(
            y, search_interval=(8, 25), solver="Nelder-Mead", verbosity=0
        )
        # Both within ~1 sample of truth, and within ~1 sample of each
        # other (the coarse_res = 0.25 default means ~4× that is the
        # functional resolution after fine search).
        assert abs(r_lbfgs.bat[0] - true_bat) < 1.0
        assert abs(r_nm.bat[0] - true_bat) < 1.0
        assert abs(r_lbfgs.bat[0] - r_nm.bat[0]) < 1.0


class TestEstimateBatOutputShape:
    def test_full_result_shapes(self):
        N, M = 40, 3
        true_bat = 15.0
        y = np.column_stack([_bolus_signal(N, true_bat) for _ in range(M)])
        result = estimate_bat(y, search_interval=(8, 25), verbosity=0)
        assert result.bat.shape == (M,)
        assert result.cp_int.shape == (M,)
        assert result.beta_opt.shape == (M,)
        assert result.k_opt.shape == (M,)
        assert result.score_opt.shape == (M,)
        assert result.y_hat.shape == (N, M)
        assert result.y_hat_x.shape == (N, M)

    def test_cp_int_is_floor_of_bat(self):
        N = 40
        y = _bolus_signal(N, 15.5)
        result = estimate_bat(y, search_interval=(8, 25), verbosity=0)
        assert result.cp_int[0] == int(np.floor(result.bat[0]))

    def test_y_hat_x_replaces_floor_bat_with_bat(self):
        """``y_hat_x[cp_int - 1]`` (0-based of 1-based cp_int) holds the
        fractional BAT; the rest is sample positions."""
        N = 40
        y = _bolus_signal(N, 15.0)
        result = estimate_bat(y, search_interval=(8, 25), verbosity=0)
        cp = int(result.cp_int[0])
        assert result.y_hat_x[cp - 1, 0] == pytest.approx(result.bat[0])
        # other rows are the sample axis 1..N
        before = result.y_hat_x[: cp - 1, 0]
        after = result.y_hat_x[cp:, 0]
        np.testing.assert_array_equal(before, np.arange(1, cp))
        np.testing.assert_array_equal(after, np.arange(cp + 1, N + 1))


class TestEstimateBatValidation:
    def test_invalid_solver_raises(self):
        with pytest.raises(ValueError, match="solver"):
            estimate_bat(np.zeros(40), solver="bad", verbosity=0)

    def test_invalid_verbosity_raises(self):
        with pytest.raises(ValueError, match="verbosity"):
            estimate_bat(np.zeros(40), verbosity=5)

    def test_short_signal_raises(self):
        with pytest.raises(ValueError, match="signal too short"):
            estimate_bat(np.zeros(5), verbosity=0)

    def test_invalid_orders_raises(self):
        with pytest.raises(ValueError, match="orders"):
            estimate_bat(np.zeros(40), orders=(0, 1, 2), verbosity=0)
