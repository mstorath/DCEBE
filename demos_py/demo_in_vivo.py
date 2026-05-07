"""BAT estimation on real in-vivo rat DCE-MRI voxel signals.

Mirrors ``demos/DCEBE_demo_in_vivo.m`` with a restricted search interval
``[40, 60]`` so the runtime stays under a minute. Estimates BAT for 94
voxels within one tumour ROI.

Run from the repository root::

    python demos_py/demo_in_vivo.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

import scipy.io

from dcebe import estimate_bat

REPO_ROOT = Path(__file__).resolve().parent.parent
IN_VIVO_MAT = REPO_ROOT / "demos" / "DCEBE_Rat_in_vivo.mat"


def main():
    data = scipy.io.loadmat(IN_VIVO_MAT)
    y_full = data["input_data"]
    # Subsample voxels to keep runtime reasonable on pure-Python.
    # MATLAB demo uses all 94 voxels; we pick the first 10.
    n_voxels_demo = 10
    y = y_full[:, :n_voxels_demo]
    N, M = y.shape
    delta_t = 0.75

    print(f"Estimating BAT for {M} in-vivo voxel signals (N={N}) ...")
    print("Using restricted search_interval=(40, 60) to keep runtime reasonable.")

    import time

    t0 = time.perf_counter()
    result = estimate_bat(y, search_interval=(40, 60), verbosity=1)
    elapsed = time.perf_counter() - t0
    bat_sec = (result.bat - 1.0) * delta_t

    print(
        f"\nElapsed: {elapsed:.1f} s. "
        f"BAT range: [{bat_sec.min():.2f}, {bat_sec.max():.2f}] s, "
        f"median {np.median(bat_sec):.2f} s, std {bat_sec.std():.2f} s."
    )

    if plt is not None:
        t = np.arange(N) * delta_t
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(t, y, ".-", alpha=0.4, linewidth=0.7)
        ax.plot(bat_sec, np.zeros_like(bat_sec), "xg", markersize=8, label="Estimated BATs")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Concentration")
        ax.legend()
        out_path = REPO_ROOT / "demos_py" / "demo_in_vivo_output.png"
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        print(f"\nPlot saved to {out_path}")


if __name__ == "__main__":
    main()
