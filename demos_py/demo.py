"""Per-signal BAT estimation on simulated rat DCE-MRI signals.

Mirrors ``demos/DCEBE_demo.m`` with reduced ``M = 5`` and a restricted
search interval so the Python version finishes in ~1 minute. (The pure-
Python implementation is correct but per-call overhead is higher than
MATLAB's; performance tuning is a post-v0.1.0 concern. Increase ``M``
or remove ``search_interval`` for a heavier benchmark.)

Run from the repository root::

    python demos_py/demo.py
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


def _add_gaussian_noise_snr(y: np.ndarray, snr: float, rng: np.random.Generator):
    """Mirror of ``demos/DCEBE_add_gaussian_noise_SNR.m``: adds Gaussian
    noise with sigma = max(signal) / (snr * std(noise)) per signal."""
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
    curve_type = 0  # MATLAB index 1 -> Python index 0
    concentration_hr = input_data[:, curve_type : curve_type + 1]
    delta_t_hr = 0.25
    N_hr = concentration_hr.shape[0]

    # True BAT: last sample of high-resolution signal still at baseline
    diffs = np.abs(np.diff(concentration_hr.ravel()))
    first_change = int(np.argmax(diffs > 0))  # 0-based index of first diff
    bat_true_sec = first_change * delta_t_hr

    # Simulate lower temporal resolution
    delta_t_lr = 2.0
    ratio = int(delta_t_lr / delta_t_hr)
    concentration_lr = concentration_hr[::ratio, :]
    N_lr = concentration_lr.shape[0]
    t_lr = np.arange(N_lr) * delta_t_lr

    rng = np.random.default_rng(123)
    SNR = 25.0
    M = 5  # MATLAB demo uses 20; we stay small to keep demo runtime under a minute
    y = _add_gaussian_noise_snr(np.tile(concentration_lr, (1, M)), SNR, rng)

    # Restrict search to the rough vicinity of the true BAT to keep
    # runtime sane on the pure-Python implementation. Remove this kwarg
    # to exercise the default-interval heuristic.
    bat_true_idx_lr = bat_true_sec / delta_t_lr + 1.0  # 1-based low-res index
    search_interval = (
        max(2.0, bat_true_idx_lr - 6.0),
        min(N_lr - 7.0, bat_true_idx_lr + 6.0),
    )

    print(f"Estimating BAT for {M} signals of length {N_lr}, search={search_interval} ...")
    result = estimate_bat(y, search_interval=search_interval, verbosity=1)
    bat_sec = (result.bat - 1.0) * delta_t_lr

    print(f"\nTrue BAT:        {bat_true_sec:.2f} s")
    print(
        f"Estimated BATs:  mean={bat_sec.mean():.2f} s, "
        f"std={bat_sec.std():.2f} s, range=[{bat_sec.min():.2f}, {bat_sec.max():.2f}] s"
    )

    if plt is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t_lr, y, ".-", alpha=0.6)
        ax.plot(bat_sec, np.zeros_like(bat_sec), "xg", markersize=10, label="Estimated BATs")
        ax.plot(bat_true_sec, 0, "ok", markersize=10, label="True BAT")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Concentration")
        ax.legend()
        out_path = REPO_ROOT / "demos_py" / "demo_output.png"
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        print(f"\nPlot saved to {out_path}")


if __name__ == "__main__":
    main()
