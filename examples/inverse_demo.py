"""Deterministic toy encoding example, not an EM simulation of a head coil.

Run: python examples/inverse_demo.py outputs/inverse-demo
This models one unresolved spatial line after the other acquisition axes have
been encoded. Replace A and C with an acquisition forward operator and measured
channel noise covariance for scientific use.
"""

import json
import sys
from pathlib import Path

import numpy as np

from mri_rf.inverse import minimum_norm_operator, point_spread, reconstruct


def run(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-1, 1, 128)
    centers = np.linspace(-1.2, 1.2, 32)
    a = np.exp(-(((x[None, :] - centers[:, None]) / 0.4) ** 2))
    a = a * np.exp(1j * centers[:, None] * x[None, :] * 2)
    c = 0.85 * np.eye(32) + 0.15 * np.ones((32, 32))
    truth = np.zeros(128)
    truth[45], truth[82] = 1, 0.7
    rng = np.random.default_rng(1234)
    noise = reconstruct(
        np.linalg.cholesky(c), rng.normal(size=32) + 1j * rng.normal(size=32)
    ) / np.sqrt(2)
    y = reconstruct(a, truth) + 0.01 * noise
    w, info = minimum_norm_operator(a, c, alpha=0.03)
    np.savez_compressed(
        folder / "input.npz", forward=a, covariance=c, data=y, positions=x, truth=truth
    )
    np.savez_compressed(
        folder / "result.npz",
        reconstruction=reconstruct(w, y),
        point_spread=point_spread(w, a, 45),
        **info,
    )
    (folder / "report.json").write_text(
        json.dumps(
            {
                "kind": "synthetic coil-like encoding; not Maxwell-derived",
                "channels": 32,
                "source_locations": 128,
                "effective_dof": float(info["effective_dof"]),
                "seed": 1234,
                "alpha": 0.03,
            },
            indent=2,
        )
        + "\n"
    )
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4), layout="constrained")
    axs[0].plot(x, a.real.T, alpha=0.45)
    axs[0].set(title="Toy channel sensitivity profiles", xlabel="Source position (arbitrary)")
    axs[1].plot(x, truth, label="Truth")
    axs[1].plot(x, np.abs(reconstruct(w, y)), label="Minimum norm")
    axs[1].legend()
    axs[1].set(title="Regularized reconstruction", xlabel="Source position (arbitrary)")
    axs[2].semilogy(info["singular_values"], "o-")
    axs[2].set(
        title=f"Whitened spectrum; effective DOF {info['effective_dof']:.1f}",
        xlabel="Mode",
        ylabel="Singular value",
    )
    fig.savefig(folder / "inverse-demo.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "outputs/inverse-demo")
