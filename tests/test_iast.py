from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fit_isotherm import fit_isotherm
from iast import iast


def test_fit_isotherm_for_sample_data():
    parameters = fit_isotherm(
        ROOT / "src" / "Results.dat-ZIF-8-433K-methane",
        ROOT / "src" / "Results.dat-ZIF-8-433K-ethane",
    )

    assert len(parameters) == 2
    assert all(p.shape == (2,) for p in parameters)
    assert np.all(np.isfinite(parameters[0]))
    assert np.all(np.isfinite(parameters[1]))


def test_iast_calculation_returns_positive_loadings():
    pressure = np.linspace(1e5, 2e6, 20)
    parameters = fit_isotherm(
        ROOT / "src" / "Results.dat-ZIF-8-433K-methane",
        ROOT / "src" / "Results.dat-ZIF-8-433K-ethane",
    )

    class SampleSystem:
        y = np.array([0.5, 0.5])
        n_components = 2
        n_pressure_intervals = len(pressure)
        isotherm_parameters = np.asarray(parameters, dtype=float)

    output = iast(pressure, SampleSystem())

    assert output.shape == (2, pressure.size)
    assert np.all(output > 0)
    assert np.all(np.isfinite(output))
