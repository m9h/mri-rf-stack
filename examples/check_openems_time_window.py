"""Assess port-spectrum sensitivity to the length of the recorded time window.

This checks S11 and accepted fraction at the benchmark frequency. It does not
establish convergence of spatial fields or the discarded pulse tail.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    report = json.loads((args.folder / "report.json").read_text())
    frequency = report["frequency_hz"]
    voltage = np.loadtxt(args.folder / "port_ut_1", comments="%")
    current = np.loadtxt(args.folder / "port_it_1", comments="%")
    windows = []
    for fraction in (0.5, 0.75, 1.0):
        spectra = []
        for data in (voltage, current):
            part = data[: int(len(data) * fraction)]
            dt = data[1, 0] - data[0, 0]
            spectra.append(
                2 * dt * np.sum(part[:, 1] * np.exp(-2j * np.pi * frequency * part[:, 0]))
            )
        u, i = spectra
        s11 = (u - 50 * i) / (u + 50 * i)
        windows.append(
            {
                "window_fraction": fraction,
                "s11": [float(s11.real), float(s11.imag)],
                "accepted_fraction": float(1 - abs(s11) ** 2),
            }
        )
    if not np.isclose(
        complex(*windows[-1]["s11"]), complex(*report["s11"]), atol=1e-10, rtol=1e-10
    ):
        raise ValueError("time-domain probes do not reproduce the reported full-window S11")
    for window in windows:
        window["absolute_s11_change_from_full_window"] = abs(
            complex(*window["s11"]) - complex(*windows[-1]["s11"])
        )
        window["relative_accepted_fraction_change"] = (
            window["accepted_fraction"] / windows[-1]["accepted_fraction"] - 1
        )
    result = {
        "frequency_hz": frequency,
        "windows": windows,
        "scope": "port spectra only; no spatial-field convergence claim",
    }
    (args.folder / "time-window.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
