from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fit_isotherm import fit_isotherm
from iast import iast
from main import build_parameters


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
    parameters = fit_isotherm(
        ROOT / "src" / "Results.dat-ZIF-8-433K-methane",
        ROOT / "src" / "Results.dat-ZIF-8-433K-ethane",
    )
    pressure, system = build_parameters(1e5, 2e6, 20, 0.5, parameters)
    assert isinstance(system, dict)
    output = iast(pressure, system)

    assert output.shape == (2, pressure.size)
    assert np.all(output > 0)
    assert np.all(np.isfinite(output))
