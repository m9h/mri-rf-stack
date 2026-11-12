"""Plot shared-grid B1 and pointwise SAR; no spatial or amplitude fitting."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mri_rf.fields import load_npz


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument(
        "--incident", action="store_true", help="show 1 W incident instead of accepted"
    )
    args = parser.parse_args()
    data = load_npz(args.folder / "common-grid.npz")
    report = json.loads((args.folder / "comparison.json").read_text())
    index = len(report["comparisons"]) - 1
    if args.incident:
        factors = {
            "hertzon": np.sqrt(report["hertzon_report"]["unscaled_accepted_w"]),
            f"openems_{index}": np.sqrt(
                report["comparisons"][index]["openems_report"]["unscaled_accepted_w"]
            ),
        }
        for prefix, scale in factors.items():
            for field in ("E_V_m", "H_A_m", "Bplus_T"):
                data[f"{prefix}_{field}"] *= scale
    plane = np.isclose(data["points_m"][:, 2], 0)
    points = data["points_m"][plane]
    axes = [np.unique(points[:, i]) for i in (0, 1)]
    ix, iy = [np.searchsorted(axis, points[:, i]) for i, axis in enumerate(axes)]

    def raster(values):
        grid = np.full((len(axes[1]), len(axes[0])), np.nan)
        grid[iy, ix] = values[plane]
        return grid

    b_fem = raster(abs(data["hertzon_Bplus_T"]) * 1e6)
    b_fdtd = raster(abs(data[f"openems_{index}_Bplus_T"]) * 1e6)
    # The benchmark uses sigma=.5 S/m and density=1000 kg/m^3. Recover these
    # from the supplied config rather than silently hardcoding phantom values.
    case = report["comparisons"][index]
    config = report["physical_config"]
    factor = config["sigma_s_m"] / (2 * config["density_kg_m3"])
    sar_fem = raster(factor * np.sum(abs(data["hertzon_E_V_m"]) ** 2, axis=1))
    sar_fdtd = raster(factor * np.sum(abs(data[f"openems_{index}_E_V_m"]) ** 2, axis=1))
    fig, axs = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    half_pixel = report["common_interior"]["spacing_m"] / 2
    extent = [
        (axes[0][0] - half_pixel) * 1e3,
        (axes[0][-1] + half_pixel) * 1e3,
        (axes[1][0] - half_pixel) * 1e3,
        (axes[1][-1] + half_pixel) * 1e3,
    ]
    for row, (a, b, label) in enumerate(
        ((b_fem, b_fdtd, "|B1+| (µT)"), (sar_fem, sar_fdtd, "Pointwise SAR (W/kg)"))
    ):
        maximum = max(np.nanmax(a), np.nanmax(b))
        for col, values in enumerate((a, b)):
            im = axs[row, col].imshow(
                values,
                origin="lower",
                extent=extent,
                vmin=0,
                vmax=maximum,
                cmap="magma",
                interpolation="none",
            )
        fig.colorbar(im, ax=axs[row, :2], shrink=0.85, label=label)
        difference = b - a
        limit = np.nanmax(abs(difference))
        im = axs[row, 2].imshow(
            difference,
            origin="lower",
            extent=extent,
            vmin=-limit,
            vmax=limit,
            cmap="RdBu_r",
            interpolation="none",
        )
        fig.colorbar(im, ax=axs[row, 2], shrink=0.85, label=f"Δ {label}")
    for ax, title in zip(
        axs[0],
        (
            "Hertzon",
            f"openEMS ({case['openems_report']['cell_m'] * 1000:g} mm)",
            "openEMS − Hertzon",
        ),
        strict=True,
    ):
        ax.set_title(title)
    for ax in axs.flat:
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_aspect("equal")
    inset = report["common_interior"]["inset_m"] * 1e3
    normalization = "incident" if args.incident else "accepted"
    fig.suptitle(
        f"{config['b0_t']:g} T birdcage · 1 W {normalization} · z = 0\n"
        f"Shared interior grid, {inset:g} mm inside phantom boundary"
    )
    filename = "fields-incident.png" if args.incident else "fields-comparison.png"
    fig.savefig(args.folder / filename, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
