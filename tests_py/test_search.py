"""Unit tests for ``dcebe._search`` (signal_intervals + default_interval).

The most important test here is ``test_default_interval_uses_envelope``,
which guards against regression of the linear-indexing bug fixed in
commit 39ff980 (reports/06-dcebe-bug-fixes.md).
"""
from __future__ import annotations

import numpy as np
import pytest

from dcebe._search import default_interval, signal_intervals


def _toy_signal(N: int, lb: int, rb: int, baseline: float = 0.0, peak: float = 1.0):
    """Construct a 1-D signal with min at index ``lb`` (1-based) and
    max at index ``rb`` (1-based)."""
    s = np.full(N, (baseline + peak) / 2.0)
    s[lb - 1] = baseline
    s[rb - 1] = peak
    return s


def _bolus_like(N: int, lb: int, rb: int, baseline: float = 0.0, peak: float = 1.0):
    """Monotone-ish DCE-shaped signal: flat at ``baseline`` until ``lb``,
    then ramping up to ``peak`` at ``rb``, flat thereafter."""
    out = np.full(N, baseline)
    if rb > lb:
        ramp = np.linspace(baseline, peak, rb - lb + 1)
        out[lb - 1 : rb] = ramp
    out[rb - 1 :] = peak
    out[lb - 1] = baseline
    return out


class TestSignalIntervals:
    def test_shape(self):
        # Use bolus-shaped signals (argmin < argmax) so the refinement
        # step is well-defined; arbitrary Gaussian noise can have argmin
        # > argmax and trigger the degenerate-fallback path.
        N, M = 40, 5
        rng = np.random.default_rng(0)
        y = np.column_stack(
            [_bolus_like(N, lb=8 + i, rb=20 - i) for i in range(M)]
        )
        y += 0.01 * rng.standard_normal(y.shape)
        out = signal_intervals(y)
        assert out.shape == (2, M)

    def test_1d_input_treated_as_one_signal(self):
        y = _bolus_like(20, lb=5, rb=12)
        out = signal_intervals(y)
        assert out.shape == (2, 1)

    def test_degenerate_signal_falls_back_to_min_max_bracket(self):
        """If argmin > argmax (e.g. monotone-decreasing signal), the
        refinement step is skipped and we return the un-refined bracket
        (with the lb -= 4 padding still applied)."""
        y = np.linspace(1.0, 0.0, 30)  # monotone-decreasing: argmin=29, argmax=0
        out = signal_intervals(y, method="min-max-refined")
        # Should not raise; rb stays at original argmax, lb shifted by -4
        assert out.shape == (2, 1)

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="unknown method"):
            signal_intervals(np.zeros(10), method="foo")

    def test_min_max_brackets(self):
        """``min-max`` method returns argmin/argmax (clamped to
        [2, N/2])."""
        N = 40
        y = _toy_signal(N, lb=10, rb=18)
        out = signal_intervals(y, method="min-max")
        np.testing.assert_array_equal(out[:, 0], [10, 18])

    def test_min_max_clamps_lower_bound_to_2(self):
        N = 40
        y = _toy_signal(N, lb=1, rb=18)  # min at index 1 (would be lb=1)
        out = signal_intervals(y, method="min-max")
        assert out[0, 0] == 2

    def test_min_max_clamps_upper_bound_to_half_N(self):
        N = 40
        y = _toy_signal(N, lb=5, rb=39)  # max at index 39 (> N/2 = 20)
        out = signal_intervals(y, method="min-max")
        assert out[1, 0] == 20.0


class TestDefaultInterval:
    def test_default_interval_shape(self):
        N, M = 40, 5
        y = np.column_stack(
            [_bolus_like(N, lb=8 + i, rb=20 - i) for i in range(M)]
        )
        lo, hi = default_interval(y)
        assert isinstance(lo, float)
        assert isinstance(hi, float)
        assert lo <= hi

    def test_default_interval_uses_envelope(self):
        """The bug-fix: with two signals whose individual intervals
        differ, the envelope must take the *minimum* lower bound and
        the *maximum* upper bound (not just signal 0's).
        """
        N = 40
        sig0 = _toy_signal(N, lb=8, rb=14)
        sig1 = _toy_signal(N, lb=12, rb=18)
        y = np.column_stack([sig0, sig1])
        per_signal = signal_intervals(y, method="min-max")
        # signal 0 has interval [8, 14], signal 1 has interval [12, 18]
        np.testing.assert_array_equal(per_signal[:, 0], [8, 14])
        np.testing.assert_array_equal(per_signal[:, 1], [12, 18])

        lo, hi = default_interval(y, method="min-max")
        # Bug-fixed envelope: lo = min(8, 12) = 8, hi = max(14, 18) = 18
        assert lo == 8.0
        assert hi == 18.0
        # The bug would have produced lo = 8, hi = 14 (signal 0 only).

    def test_default_interval_single_signal(self):
        """For a single signal, the envelope equals the per-signal
        interval (no aggregation surprise)."""
        N = 40
        y = _toy_signal(N, lb=8, rb=18)
        lo, hi = default_interval(y, method="min-max")
        assert (lo, hi) == (8.0, 18.0)
