from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

try:
    from .isotherms import langmuir_fit
except ImportError:  # pragma: no cover - direct script execution
    from isotherms import langmuir_fit

PRESSURE_COLUMN = 0
LOADING_COLUMN = 7
HEADER_LINES_TO_SKIP = 22


def fit_isotherm(*file_names):
    """Fit a Langmuir isotherm to each RASPA output file."""
    coefficients = []
    for file_name in file_names:
        path = Path(file_name)
        if not path.is_file():
            raise FileNotFoundError(f"File '{file_name}' does not exist")

        pressure_loading = _read_pressure_loading_data(path)
        coefficients.append(_fit_single_isotherm(pressure_loading))
    return coefficients


def fitIsotherm(*file_names):
    return fit_isotherm(*file_names)


def _read_pressure_loading_data(path: Path):
    data = np.genfromtxt(
        path,
        comments="#",
        skip_header=HEADER_LINES_TO_SKIP,
        usecols=(PRESSURE_COLUMN, LOADING_COLUMN),
        invalid_raise=False,
    )
    if data.size == 0:
        raise ValueError(f"No valid adsorption data found in '{path}'")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _fit_single_isotherm(pressure_loading):
    pressure = pressure_loading[:, 0]
    loading = pressure_loading[:, 1]
    if np.any(pressure <= 0) or np.any(loading < 0):
        raise ValueError("Pressure and loading values must be positive and finite")

    initial_guess = [loading[-1], loading[1] / pressure[1]]
    parameters, _ = curve_fit(langmuir_fit, pressure, loading, p0=initial_guess)
    return parameters
