# DCEBE — Python package

`dcebe` is a Python re-implementation of the MATLAB DCEBE reference for bolus arrival time (BAT) estimation in DCE-MRI signals. It provides one user-facing function — `estimate_bat` — that mirrors `DCEBE_estimateBAT.m` and returns the same outputs in 1-based MATLAB-compatible sample units.

## Install

```bash
pip install dcebe
```

The package is pure Python (NumPy + SciPy); no compiled wheels, no Rust, no MATLAB runtime. Wheels are available for Python 3.9–3.13 on Linux, macOS, and Windows.

For development:

```bash
git clone https://github.com/mstorath/DCEBE.git
cd DCEBE
pip install -e .[test]
pytest tests_py/
```

## Quickstart

```python
import numpy as np
from dcebe import estimate_bat

# y: (N,) or (N, M) array of DCE-MRI signals
result = estimate_bat(y, search_interval=(10, 30))

# 1-based BAT in sample units; convert to seconds via:
bat_sec = (result.bat - 1) * delta_t
```

The full output is on `result`:

| Field           | Shape    | Meaning                                                                |
|-----------------|----------|------------------------------------------------------------------------|
| `result.bat`    | `(M,)`   | Estimated BAT, 1-based sample units                                    |
| `result.cp_int` | `(M,)`   | `floor(bat)` — integer changepoint                                     |
| `result.beta_opt` | `(M,)` | Optimal stiffness parameter                                            |
| `result.k_opt`  | `(M,)`   | Optimal spline order                                                   |
| `result.score_opt` | `(M,)` | Final log-GCV score                                                  |
| `result.y_hat`  | `(N, M)` | Reconstructed signal                                                   |
| `result.y_hat_x`| `(N, M)` | Sample x-axis with `floor(bat)` row replaced by fractional `bat`       |

## API reference

```python
estimate_bat(
    y,                                # (N,) or (N, M) array of signals
    *,
    search_interval=None,             # (lo, hi) in 1-based samples; None = auto
    coarse_res=0.25,                  # coarse-grid spacing
    beta_coarse_search=None,          # default: linspace(1, 25, 50)
    orders=(3, 4, 5, 6),              # spline orders to evaluate
    common_bat=False,                 # joint estimation across signals
    solver="L-BFGS-B",                # or "Nelder-Mead"
    verbosity=1,                      # 0 silent, 1 per-order, 2 per-signal
) -> EstimateResult
```

## Index conventions

`bat`, `cp_int`, and the `search_interval` argument are all in **1-based** MATLAB-compatible sample units. This makes parity tests against the MATLAB reference `np.allclose`-direct, at the cost of being slightly off-idiom for Python users. Convert to seconds with `(bat - 1) * delta_t`, exactly as the MATLAB demos do.

## Demos

Three runnable demonstrations live in [`demos_py/`](demos_py/):

| Script | Mirrors |
|---|---|
| [`demos_py/demo.py`](demos_py/demo.py) | `demos/DCEBE_demo.m` (per-signal mode, simulated rat data, SNR=25) |
| [`demos_py/demo_common.py`](demos_py/demo_common.py) | `demos/DCEBE_demo_commonBAT.m` (common-BAT mode, SNR=10) |
| [`demos_py/demo_in_vivo.py`](demos_py/demo_in_vivo.py) | `demos/DCEBE_demo_in_vivo.m` (94 voxels of real rat data, restricted search interval) |

Run with `python demos_py/demo.py`. Each prints estimated BATs and shows a Matplotlib plot if `matplotlib` is installed (`pip install dcebe[demos]`).

## When to use this vs. alternatives

`dcebe` is **focused**: bolus arrival time only, spline-based, GCV-tuned, designed for signals **without a fast upslope** (typical of small-animal DCE-MRI, where the bolus shape is smooth rather than abrupt). It is intentionally narrower than the broad DCE-MRI pipelines listed below.

### Related work

- **[OSIPI Python package (`osipy`)](https://github.com/OSIPI/osipy) and [DCE-DSC-MRI_CodeCollection](https://github.com/OSIPI/DCE-DSC-MRI_CodeCollection)** — the ISMRM Open Science Initiative for Perfusion Imaging maintains a comprehensive Python collection (86 implementations across BAT estimation, T1 mapping, AIF fitting, and pharmacokinetic models like Tofts / Extended Tofts / Patlak / 2CXM / 2CUM). If you need a full DCE-MRI pipeline, start there. `dcebe` complements OSIPI for the specific case of small-animal-style smooth-upslope BAT estimation.
- **[`pydcemri`](https://github.com/welcheb/pydcemri) (and the [`aydindemircioglu` fork](https://github.com/aydindemircioglu/pydcemri))** — older Python module producing K^trans, v_e, v_p maps from a T1 + dynamic + AIF stack. Different scope (full pharmacokinetic modeling vs. BAT-only) but a useful upstream reference.
- **[`DCEMRI.jl`](https://pmc.ncbi.nlm.nih.gov/articles/PMC4411523/)** — Julia toolkit for DCE-MRI analysis, not Python; included for completeness as a reference implementation in the perfusion-imaging space.
- **[MITK-ModelFit](https://link.springer.com/article/10.1186/s12859-018-2588-1)** — C++ framework with Python bindings, focused on parameter fitting in medical imaging more generally.

If your data has a sharp upslope (typical of clinical human DCE-MRI), the OSIPI BAT estimators may fit better than `dcebe`'s spline-smoothed approach. The DCEBE paper specifically targets the *small-animal* regime where the bolus profile is gradual — see the citation below.

## Citation

If you use `dcebe` in scientific work, please cite the underlying paper:

> A. Bendinger, C. Debus, C. Glowa, C. Karger, J. Peter, M. Storath, *Bolus arrival time estimation in dynamic contrast-enhanced MRI of small animals based on spline models*, Physics in Medicine & Biology 64(4), 2019. [DOI: 10.1088/1361-6560/aafce7](https://doi.org/10.1088/1361-6560/aafce7).

The `CITATION.cff` at the repository root provides machine-readable metadata.

## License

`dcebe` is released under the Apache License 2.0. The vendored `fdweights` algorithm by Toby Driscoll is BSD-2-Clause; its license is reproduced in `external/fdweights/license.txt` and the `fdweights` docstring.

## How the Python implementation differs from MATLAB

- **Numerical kernel**: same algorithm; SciPy LAPACK calls (QR, SVD, eigendecomposition) replace MATLAB's `mldivide` / `eig` / `svd`.
- **Bug fixes baked in**: the linear-indexing bug at `DCEBE_estimateBAT.m:50` (default search interval used only signal 1's bounds) is fixed in the canonical Python `default_interval`.
- **Default optimiser**: `scipy.optimize.minimize(method="L-BFGS-B")`. The MATLAB reference uses `fminunc` quasi-Newton with central finite differences; the SciPy version uses forward FD with a slightly enlarged step (`eps=1e-4`) and a sub-grid jitter on integer-`t` starting points to bridge the kink the spline-matrix structure has there.
- **`factorize=True` path**: ported as `dcebe._factorize.gcv_score_factorize` but **not exposed** through `estimate_bat`; the MATLAB source itself notes the QR path is better-conditioned.

## Status

`0.1.0.dev0`. Parity tested against MATLAB on three pinned-seed configurations:

- `atol(BAT) = 1e-2` (a hundredth of a sample)
- `atol(score_opt) = 1e-3` (load-bearing — both solvers find the same minimum value)
- exact match on `k_opt`
- `rtol(beta_opt) = 0.2` (the GCV objective is flat in beta near the optimum)

See `tests_py/test_matlab_parity.py` for the parity assertions.
