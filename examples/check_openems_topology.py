"""Audit the actual PEC edge graph exported by openEMS (requires PyVista)."""

import argparse
import json
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    mesh = pv.read(args.folder / "PEC_dump.vtp")
    lines = mesh.lines.reshape(-1, 3)
    if not np.all(lines[:, 0] == 2):
        raise ValueError("expected two-point PEC edges")
    points, inverse = np.unique(np.round(mesh.points, 10), axis=0, return_inverse=True)
    edges = inverse[lines[:, 1:]]
    graph = coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(len(points), len(points))
    )
    count, labels = connected_components(graph, directed=False)
    components = []
    for label in range(count):
        p = points[labels == label]
        components.append(
            {"nodes": len(p), "min_m": p.min(axis=0).tolist(), "max_m": p.max(axis=0).tolist()}
        )
    report = {
        "pec_edges": len(edges),
        "components": components,
        "expected_components": "shield, lower ring/rungs, upper ring/rungs, feed-to-capacitor strip",
    }
    (args.folder / "pec-topology.json").write_text(json.dumps(report, indent=2) + "\n")
    if count != 4:
        raise ValueError(f"expected four isolated PEC networks, found {count}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
