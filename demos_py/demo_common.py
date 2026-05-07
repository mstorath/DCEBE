"""Common-BAT estimation: a single shared BAT for all signals.

Mirrors ``demos/DCEBE_demo_commonBAT.m``. SNR=10 with 20 noisy
realisations of one DCE-MRI curve, jointly estimating one BAT.

Run from the repository root::

    python demos_py/demo_common.py
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
ETM_MAT = REPO_ROOT / "demos" / "DCEBE_Rat_ETM.mat"


def _add_gaussian_noise_snr(y, snr, rng):
    M = y.shape[1]
    out = np.zeros_like(y)
    for m in range(M):
        sig = y[:, m]
        noise = rng.standard_normal(sig.shape)
        sigma = float(np.max(sig)) / (snr * float(np.std(noise)))
        out[:, m] = sig + sigma * noise
    return out


def main():
    data = scipy.io.loadmat(ETM_MAT)
    input_data = data["input_data"]
    curve_type = 1  # MATLAB index 2 -> Python index 1
    concentration_hr = input_data[:, curve_type : curve_type + 1]
    delta_t_hr = 0.25

    diffs = np.abs(np.diff(concentration_hr.ravel()))
    bat_true_sec = int(np.argmax(diffs > 0)) * delta_t_hr

    delta_t_lr = 2.0
    ratio = int(delta_t_lr / delta_t_hr)
    concentration_lr = concentration_hr[::ratio, :]
    N_lr = concentration_lr.shape[0]
    t_lr = np.arange(N_lr) * delta_t_lr

    rng = np.random.default_rng(123)
    SNR = 10.0
    M = 8  # MATLAB demo uses 20; smaller M keeps demo runtime reasonable
    y = _add_gaussian_noise_snr(np.tile(concentration_lr, (1, M)), SNR, rng)

    bat_true_idx = int(np.argmax(np.abs(np.diff(concentration_hr.ravel())) > 0)) + 1
    bat_true_idx_lr = bat_true_idx // ratio + 1
    search_interval = (
        max(2.0, bat_true_idx_lr - 8.0),
        min(N_lr - 7.0, bat_true_idx_lr + 8.0),
    )

    print(f"Estimating common BAT for {M} signals of length {N_lr}, search={search_interval} ...")
    result = estimate_bat(y, search_interval=search_interval, common_bat=True, verbosity=1)
    bat_sec = (result.bat[0] - 1.0) * delta_t_lr  # all entries equal

    print(f"\nTrue BAT:       {bat_true_sec:.2f} s")
    print(f"Estimated BAT:  {bat_sec:.2f} s")

    if plt is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t_lr, y, ".-", alpha=0.6)
        ax.plot([bat_sec], [0.0], "xg", markersize=12, label="Estimated BAT")
        ax.plot([bat_true_sec], [0.0], "ok", markersize=12, label="True BAT")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Concentration")
        ax.legend()
        out_path = REPO_ROOT / "demos_py" / "demo_common_output.png"
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        print(f"\nPlot saved to {out_path}")


if __name__ == "__main__":
    main()
