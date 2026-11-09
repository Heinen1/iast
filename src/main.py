from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    from .fit_isotherm import fit_isotherm
    from .iast import iast
except ImportError:  # pragma: no cover - direct script execution
    from fit_isotherm import fit_isotherm
    from iast import iast


DEFAULT_DATA_FILES = (
    "src/Results.dat-ZIF-8-433K-methane",
    "src/Results.dat-ZIF-8-433K-ethane",
)


def _resolve_data_file(name: str | None, fallback: str) -> Path:
    if name:
        return Path(name).expanduser().resolve()
    return (Path(__file__).resolve().parent.parent / fallback).resolve()


def build_parameters(
    pressure_start: float,
    pressure_end: float,
    points: int,
    y1: float,
    isotherm_parameters,
) -> tuple[np.ndarray, dict]:
    y2 = 1.0 - y1
    pressure = np.linspace(pressure_start, pressure_end, points)
    parameters = {
        "pressure_initial": pressure_start,
        "isotherm_parameters": np.asarray(isotherm_parameters, dtype=float),
        "y": np.asarray([y1, y2], dtype=float),
        "n_components": 2,
        "n_pressure_intervals": points,
    }
    return pressure, parameters


def _write_csv(path: Path | None, pressure: np.ndarray, loading: np.ndarray) -> None:
    if path is None:
        return
    with path.open("w", encoding="utf-8") as handle:
        handle.write("pressure_pa,component_1_loading,component_2_loading\n")
        for p, q1, q2 in zip(pressure, loading[0], loading[1]):
            handle.write(f"{p:.12e},{q1:.12e},{q2:.12e}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute binary IAST adsorption using Langmuir fits.")
    parser.add_argument("--file1", type=str, help="Path to the first adsorption-isotherm file.")
    parser.add_argument("--file2", type=str, help="Path to the second adsorption-isotherm file.")
    parser.add_argument("--y1", type=float, default=0.5, help="Molar fraction of component 1 in the gas phase.")
    parser.add_argument("--pressure-start", type=float, default=1e5, help="Starting pressure in Pa.")
    parser.add_argument("--pressure-end", type=float, default=2e6, help="Final pressure in Pa.")
    parser.add_argument("--points", type=int, default=20, help="Number of pressure points to evaluate.")
    parser.add_argument("--output", type=str, help="Optional CSV output path.")
    args = parser.parse_args()

    file1 = _resolve_data_file(args.file1, DEFAULT_DATA_FILES[0])
    file2 = _resolve_data_file(args.file2, DEFAULT_DATA_FILES[1])
    isotherm_parameters = fit_isotherm(file1, file2)
    pressure, parameters = build_parameters(
        args.pressure_start,
        args.pressure_end,
        args.points,
        args.y1,
        isotherm_parameters,
    )
    loading = iast(pressure, parameters)

    _write_csv(Path(args.output).expanduser().resolve() if args.output else None, pressure, loading)

    print("Pressure [Pa], Component 1 [mol/kg], Component 2 [mol/kg]")
    for p, q1, q2 in zip(pressure, loading[0], loading[1]):
        print(f"{p:.12e}, {q1:.12e}, {q2:.12e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
