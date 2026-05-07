"""Search-interval heuristics for the bolus arrival time.

Translation of ``DCEBE_searchIntervals.m`` plus the call-site logic in
``DCEBE_estimateBAT.m`` that took the envelope across signals. The
linear-indexing bug at that call site (``interval_cand(1)`` instead of
``interval_cand(1,:)``) is fixed in the canonical
:func:`default_interval` here — see reports/06-dcebe-bug-fixes.md.

All bounds are returned in MATLAB-1-based sample units, matching the
rest of the API.
"""
from __future__ import annotations

import numpy as np

_METHODS = ("min-max", "min-max-refined")


def signal_intervals(y, method: str = "min-max-refined") -> np.ndarray:
    """Per-signal search intervals.

    Parameters
    ----------
    y : array_like, shape ``(N,)`` or ``(N, M)``
        DCE signals as columns. A 1-D input is treated as ``M = 1``.
    method : {"min-max", "min-max-refined"}
        Heuristic used. ``"min-max"`` brackets between argmin and argmax;
        ``"min-max-refined"`` (default) tightens the upper bound to the
        first sample whose value is closest to the midpoint between
        argmin and argmax, and pads the lower bound by 4 samples.

    Returns
    -------
    intervals : ndarray, shape ``(2, M)``
        Row 0: lower bounds (1-based). Row 1: upper bounds (1-based).
    """
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2:
        raise ValueError(f"y must be 1-D or 2-D, got shape {y.shape}")
    if method not in _METHODS:
        raise ValueError(f"unknown method {method!r}; choose from {_METHODS}")
    N, M = y.shape

    minb = 2.0
    maxb = N / 2.0

    intervals = np.zeros((2, M))
    for m in range(M):
        col = y[:, m]
        lb = int(np.argmin(col)) + 1
        rb = int(np.argmax(col)) + 1
        if method == "min-max-refined":
            # Refinement: tighten rb to the first sample within [lb, rb]
            # whose value is closest to the midpoint between argmin and
            # argmax. Defined only when lb < rb; for degenerate signals
            # (lb >= rb, e.g. random noise), skip refinement and fall
            # back to the un-refined argmin/argmax bracket.
            if lb < rb:
                tr = (col[lb - 1] + col[rb - 1]) / 2.0
                sub = col[lb - 1 : rb]
                rb_aux = int(np.argmin(np.abs(sub - tr))) + 1
                rb = rb_aux + lb - 1
            lb = lb - 4
        intervals[0, m] = max(minb, lb)
        intervals[1, m] = min(maxb, rb)

    return intervals


def default_interval(y, method: str = "min-max-refined") -> tuple[float, float]:
    """Envelope of :func:`signal_intervals` across all signals.

    This is the bug-fixed analog of the call site in DCEBE_estimateBAT.m
    lines 49–51 (post bug-fix commit ``39ff980``): the union of
    per-signal lower and upper bounds, not just signal 1.
    """
    intervals = signal_intervals(y, method=method)
    return float(np.min(intervals[0, :])), float(np.max(intervals[1, :]))
