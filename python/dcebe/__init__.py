"""DCEBE: bolus arrival time estimation for DCE-MRI signals.

Spline-based model for estimating the bolus arrival time of dynamic
contrast-enhanced MRI signals, particularly intended for signals without
a fast upslope (typical for small-animal data). Parameters are selected
via generalized cross-validation.

Reference:
    A. Bendinger, C. Debus, C. Glowa, C. Karger, J. Peter, M. Storath,
    "Bolus arrival time estimation in dynamic contrast-enhanced MRI of
    small animals based on spline models", Physics in Medicine & Biology
    64(4), 2019. DOI: 10.1088/1361-6560/aafce7

The user-facing API is :func:`estimate_bat`. The current build is
``0.1.0.dev0`` and exposes only the spline-construction primitives;
:func:`estimate_bat` lands in a later port day.
"""
from __future__ import annotations

__version__ = "0.1.0.dev0"
