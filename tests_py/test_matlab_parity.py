"""MATLAB parity tests for ``dcebe.estimate_bat``.

Loads pre-generated ``.mat`` fixtures from ``matlab_fixtures/`` (regen
via ``scripts/gen_fixtures.m``) and asserts that the Python output
matches the MATLAB reference within tolerance.

Tolerance philosophy: SciPy's L-BFGS-B uses forward finite differences
by default whereas MATLAB's ``fminunc`` uses central, so the fine
optima can diverge at the ~1e-3 level. We therefore use:

  - ``atol=1e-2`` on BAT (a hundredth of a sample)
  - ``rtol=0.2`` on ``beta_opt`` (the GCV objective is flat in beta
    near the optimum; both solvers find valid local minima within the
    same basin but at beta values that can differ by up to ~15-20%)
  - exact match on ``k_opt`` (integer choice)
  - ``atol=1e-3`` on ``score_opt`` (this is the load-bearing check —
    both solvers converged to the same minimum value)

These tolerances are deliberately looser than what a Rust/MEX port
would need; they validate "same algorithm, modulo solver micro-
differences" rather than bitwise equivalence.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dcebe import estimate_bat

FIXTURE_DIR = Path(__file__).parent.parent / "matlab_fixtures"


def _load(name: str) -> dict:
    if not FIXTURE_DIR.exists() or not (FIXTURE_DIR / f"{name}.mat").exists():
        pytest.skip(
            f"fixture {name} missing; regenerate via "
            "`matlab -nodisplay -batch \"run('scripts/gen_fixtures.m')\"`"
        )
    raw = scipy.io.loadmat(FIXTURE_DIR / f"{name}.mat")
    return {
        "y": raw["y"],
        "bat": raw["bat_matlab"].ravel(),
        "beta_opt": raw["beta_opt"].ravel(),
        "k_opt": raw["k_opt"].ravel().astype(int),
        "score_opt": raw["score_opt"].ravel(),
        "y_hat": raw["y_hat"],
        "search_interval": (
            tuple(np.asarray(raw["search_interval"]).ravel().astype(float))
            if "search_interval" in raw
            else None
        ),
    }


def _compare(result, ref, *, atol_bat=1e-2, atol_score=1e-3, rtol_beta=0.2):
    np.testing.assert_allclose(result.bat, ref["bat"], atol=atol_bat)
    np.testing.assert_array_equal(result.k_opt, ref["k_opt"])
    np.testing.assert_allclose(result.beta_opt, ref["beta_opt"], rtol=rtol_beta)
    np.testing.assert_allclose(result.score_opt, ref["score_opt"], atol=atol_score)


class TestParity:
    def test_fixture_small(self):
        ref = _load("fixture_small")
        result = estimate_bat(ref["y"], verbosity=0)
        _compare(result, ref)

    def test_fixture_explicit_si(self):
        ref = _load("fixture_explicit_si")
        result = estimate_bat(
            ref["y"], search_interval=ref["search_interval"], verbosity=0
        )
        _compare(result, ref)

    def test_fixture_common(self):
        ref = _load("fixture_common")
        result = estimate_bat(
            ref["y"],
            search_interval=ref["search_interval"],
            common_bat=True,
            verbosity=0,
        )
        _compare(result, ref)

    def test_fixture_small_y_hat(self):
        """Reconstructed signal y_hat should also match (looser
        tolerance because it depends on solver-found beta and bat)."""
        ref = _load("fixture_small")
        result = estimate_bat(ref["y"], verbosity=0)
        np.testing.assert_allclose(result.y_hat, ref["y_hat"], atol=5e-2)
