"""Plot a central tissue slab from the output of mri-rf analyze."""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(folder):
    folder = Path(folder)
    with np.load(folder / "fields.npz") as d:
        points = d["points_m"] * 1000
        selected = np.abs(points[:, 2]) < 10
        fields = (np.abs(d["Bplus_T"]) * 1e6, np.abs(d["Bminus_T"]) * 1e6, d["sar_w_kg"])
        titles = ("|B+| (µT)", "|B−| (µT)", "Element-mean SAR (W/kg)")
        fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), layout="constrained")
        for ax, field, title in zip(axs, fields, titles, strict=True):
            dots = ax.scatter(
                points[selected, 0], points[selected, 1], c=field[selected], s=18, cmap="viridis"
            )
            ax.set(title=title, xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
            fig.colorbar(dots, ax=ax, shrink=0.8)
        fig.suptitle("Homogeneous phantom · |z| < 10 mm · 1 W accepted · research benchmark")
        fig.savefig(folder / "fields.png", dpi=170)
        plt.close(fig)


if __name__ == "__main__":
    plot(sys.argv[1])
