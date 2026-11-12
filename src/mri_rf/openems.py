"""Independent Cartesian FDTD discretization of the nominal birdcage.

openEMS/CSXCAD are optional native dependencies. PEC sheets are snapped to
Cartesian faces explicitly: a zero-thickness rotated box would otherwise miss
Yee edges, and rotating a lumped-element box does not rotate its circuit.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import numpy as np

from .model import Birdcage
from .physics import PROTON_HZ_PER_T, circular_b, finite, local_sar


def axis_mesh(anchors, step):
    """Subdivide anchor intervals without moving conductor/gap endpoints."""
    anchors = np.unique(np.round(anchors, 12))
    return np.concatenate(
        [
            np.linspace(a, b, max(1, int(np.ceil((b - a) / step - 1e-9))), endpoint=False)
            for a, b in zip(anchors[:-1], anchors[1:], strict=True)
        ]
        + [anchors[-1:]]
    )


def contour_faces(x, y, inside):
    """Grid-face segments separating inside and outside cell centres in xy."""
    segments = []
    for i, j in np.argwhere(inside[:-1] != inside[1:]):
        segments.append([[x[i + 1], y[j]], [x[i + 1], y[j + 1]]])
    for i, j in np.argwhere(inside[:, :-1] != inside[:, 1:]):
        segments.append([[x[i], y[j + 1]], [x[i + 1], y[j + 1]]])
    return np.asarray(segments)


def discretize(config, cell_m):
    """Return a mesh and zero-thickness sheet faces, plus circuit node weights."""
    config.validate()
    c = config
    if not np.isfinite(cell_m) or not 0 < cell_m <= min(c.strip_width_m / 2, c.gap_m):
        raise ValueError("cell size must be positive and <= strip width/2 and gap")
    # The air box matches Hertzon Auto Open - Filter (OCC's 0.1 um padding omitted).
    side = c.shield_radius_m + c.air_margin_m
    end = c.shield_length_m / 2 + c.air_margin_m
    anchors = [
        0,
        c.strip_width_m / 2,
        c.phantom_radius_m,
        c.radius_m - c.strip_width_m / 2,
        c.radius_m,
        c.radius_m + c.strip_width_m / 2,
        c.shield_radius_m,
    ]
    # Anchor every strip endpoint as well as its centre. Otherwise the 45-degree
    # strips can lose a substantial fraction of their width on a coarse grid.
    for index in range(c.rungs):
        angle = 2 * np.pi * index / c.rungs
        n = np.array([np.cos(angle), np.sin(angle)])
        t = np.array([-n[1], n[0]])
        for offset in (-c.strip_width_m / 2, 0, c.strip_width_m / 2):
            anchors.extend(abs(c.radius_m * n + offset * t))
    positive = axis_mesh(anchors, cell_m)
    outer = axis_mesh([c.shield_radius_m, side], 3 * cell_m)[1:]
    positive = np.r_[positive, outer]
    x = np.r_[-positive[:0:-1], positive]
    y = x.copy()
    zanchors = [
        -c.shield_length_m / 2,
        -c.length_m / 2,
        -c.phantom_length_m / 2,
        -c.length_m / 4 - c.gap_m / 2,
        -c.length_m / 4 + c.gap_m / 2,
        -c.gap_m / 2,
        0,
        c.gap_m / 2,
        c.phantom_length_m / 2,
        c.length_m / 2,
        c.shield_length_m / 2,
    ]
    z = np.r_[
        axis_mesh([-end, zanchors[0]], 3 * cell_m)[:-1],
        axis_mesh(zanchors, cell_m),
        axis_mesh([zanchors[-1], end], 3 * cell_m)[1:],
    ]
    xx, yy = np.meshgrid((x[:-1] + x[1:]) / 2, (y[:-1] + y[1:]) / 2, indexing="ij")
    shield = contour_faces(x, y, xx**2 + yy**2 < c.shield_radius_m**2)
    ring_mask = (xx**2 + yy**2 >= (c.radius_m - c.strip_width_m / 2) ** 2) & (
        xx**2 + yy**2 <= (c.radius_m + c.strip_width_m / 2) ** 2
    )
    rings = np.array([[[x[i], y[j]], [x[i + 1], y[j + 1]]] for i, j in np.argwhere(ring_mask)])
    rungs = []
    for index in range(c.rungs):
        angle = index * 2 * np.pi / c.rungs
        normal = np.array([np.cos(angle), np.sin(angle)])
        tangent = np.array([-normal[1], normal[0]])
        segments = contour_faces(x, y, normal[0] * xx + normal[1] * yy < c.radius_m)
        nodes = np.unique(segments.reshape(-1, 2), axis=0)
        order = np.argsort(nodes @ tangent)
        nodes = nodes[order]
        s = nodes @ tangent
        lo, hi = [int(np.argmin(abs(s - a))) for a in (-c.strip_width_m / 2, c.strip_width_m / 2)]
        nodes = nodes[lo : hi + 1]
        if len(nodes) < 2:
            raise ValueError("mesh does not resolve a rung; reduce cell size")
        selected = np.array(
            [
                seg
                for seg in segments
                if all(np.any(np.all(nodes == point, axis=1)) for point in seg)
            ]
        )
        if len(selected) != len(nodes) - 1:
            raise ValueError("disconnected rung staircase")
        # Distribute the total capacitor over UNIQUE axial edges. Adjacent 2-D
        # RLC boxes share edges and overwrite each other; these 1-D boxes do not.
        s = nodes @ tangent
        weights = np.diff(np.r_[s[0], (s[1:] + s[:-1]) / 2, s[-1]])
        weights /= weights.sum()
        rungs.append(
            {
                "nodes": nodes,
                "segments": selected,
                "weights": weights,
                "projected_width_m": float(s[-1] - s[0]),
                "max_plane_error_m": float(np.max(abs(nodes @ normal - c.radius_m))),
            }
        )
    return (x, y, z), shield, rings, rungs


def build(config, cell_m, end_criteria, max_steps):
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    c = config
    axes, shield, rings, rungs = discretize(c, cell_m)
    frequency = c.b0_t * PROTON_HZ_PER_T
    fdtd = openEMS(NrTS=max_steps, EndCriteria=end_criteria, CellConstantMaterial=True)
    fdtd.SetGaussExcite(frequency, 80e6)
    fdtd.SetBoundaryCond(["MUR"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1)
    for axis, lines in zip("xyz", axes, strict=True):
        grid.AddLine(axis, lines)
    phantom = csx.AddMaterial(
        "Phantom", epsilon=c.epsilon_r, kappa=c.sigma_s_m, density=c.density_kg_m3
    )
    phantom.AddCylinder(
        [0, 0, -c.phantom_length_m / 2],
        [0, 0, c.phantom_length_m / 2],
        c.phantom_radius_m,
        priority=0,
    )
    metal = csx.AddMetal("PEC")

    def sheet(segment, lo, hi):
        a, b = segment
        metal.AddBox([*a, lo], [*b, hi], priority=10)

    for segment in shield:
        sheet(segment, -c.shield_length_m / 2, c.shield_length_m / 2)
    for z in (-c.length_m / 2, c.length_m / 2):
        for rectangle in rings:
            sheet(rectangle, z, z)
    for i, rung in enumerate(rungs):
        intervals = [(-c.length_m / 2, -c.gap_m / 2), (c.gap_m / 2, c.length_m / 2)]
        if i == 0:
            lo, hi = -c.length_m / 4 - c.gap_m / 2, -c.length_m / 4 + c.gap_m / 2
            intervals = [(intervals[0][0], lo), (hi, intervals[0][1]), intervals[1]]
        for segment in rung["segments"]:
            for lo, hi in intervals:
                sheet(segment, lo, hi)
        for j, (node, weight) in enumerate(zip(rung["nodes"], rung["weights"], strict=True)):
            cap = csx.AddLumpedElement(
                f"C{i}_{j}", ny="z", caps=False, C=float(c.capacitance_f * weight)
            )
            cap.AddBox([*node, -c.gap_m / 2], [*node, c.gap_m / 2], priority=20)
    port = fdtd.AddLumpedPort(
        1,
        50,
        [c.radius_m, -c.strip_width_m / 2, -c.length_m / 4 - c.gap_m / 2],
        [c.radius_m, c.strip_width_m / 2, -c.length_m / 4 + c.gap_m / 2],
        "z",
        excite=1,
        priority=30,
    )
    # Include one air layer for material-interface accounting. Raw SAR supplies
    # native cell volumes, density, conductivity, and cell-interpolated E.
    start = [
        -c.phantom_radius_m - cell_m,
        -c.phantom_radius_m - cell_m,
        -c.phantom_length_m / 2 - cell_m,
    ]
    stop = [-v for v in start]
    for name, kind in (("E_raw", 29), ("H", 11)):
        dump = csx.AddDump(name, dump_type=kind, dump_mode=2, file_type=1, frequency=[frequency])
        dump.AddBox(start, stop)
    manifest = {
        "config": asdict(c),
        "cell_m": cell_m,
        "frequency_hz": frequency,
        "end_criteria": end_criteria,
        "max_steps": max_steps,
        "boundary": "first-order Mur on the same nominal outer box as Hertzon",
        "mesh_lines": [len(a) for a in axes],
        "mesh_cells": int(np.prod([len(a) - 1 for a in axes])),
        "box_m": [[float(a[0]), float(a[-1])] for a in axes],
        "rungs": [
            {
                "projected_width_m": r["projected_width_m"],
                "max_plane_error_m": r["max_plane_error_m"],
                "capacitor_nodes": len(r["nodes"]),
                "capacitance_f": float(np.sum(r["weights"]) * c.capacitance_f),
            }
            for r in rungs
        ],
        "sheet_discretization": "zero-thickness Cartesian staircase; no finite conductor thickness",
        "material_discretization": "cell-centre constant epsilon, conductivity and density",
        "python_versions": {p: version(p) for p in ("openEMS", "CSXCAD", "numpy", "h5py")},
    }
    return fdtd, port, manifest, (axes, shield, rings, rungs)


def analyze(folder, port, manifest):
    from openEMS.utilities import HDF5Dump

    folder = Path(folder)
    f = manifest["frequency_hz"]
    port.CalcPort(str(folder), np.array([f]))
    s11 = complex(port.uf_ref[0] / port.uf_inc[0])
    accepted_fraction = 1 - abs(s11) ** 2
    power = float(port.P_acc[0])
    if not np.isfinite(power) or power <= 0 or accepted_fraction <= 0:
        raise ValueError("nonpositive or invalid accepted spectral power")
    if not np.isclose(power / float(port.P_inc[0]), accepted_fraction, rtol=1e-7):
        raise ValueError("port power-wave and VI normalization disagree")
    # Both port and field DFT use 2*dt*sum(x*exp(-iwt)). Their pulse spectral
    # units cancel here; rotate to real-positive incident voltage as in FEM.
    scale = np.exp(-1j * np.angle(port.uf_inc[0])) / np.sqrt(power)
    with HDF5Dump(str(folder / "E_raw.h5")) as dump:
        if len(dump.frequencies) != 1 or not np.isclose(dump.frequencies[0], f, rtol=1e-9):
            raise ValueError("electric field dump frequency disagrees with configuration")
        axes = dump.GetMesh()["lines"]
        e = np.moveaxis(dump.GetFieldAtIndex(f_idx=0), 0, -1) * scale

        # Require the new native XYZ layout; never silently guess axis order.
        def cell_data(name):
            ds = dump.file[f"CellData/{name}"]
            order = ds.attrs.get("d_order", "")
            if isinstance(order, bytes):
                order = order.decode()
            if str(order).upper() != "XYZ":
                raise ValueError(f"expected native XYZ cell data, got {order!r}")
            return np.asarray(ds)

        volumes, density, sigma = [cell_data(k) for k in ("Volume", "Density", "Conductivity")]
    with HDF5Dump(str(folder / "H.h5")) as dump:
        if len(dump.frequencies) != 1 or not np.isclose(dump.frequencies[0], f, rtol=1e-9):
            raise ValueError("magnetic field dump frequency disagrees with configuration")
        if not all(
            np.array_equal(a, b) for a, b in zip(axes, dump.GetMesh()["lines"], strict=True)
        ):
            raise ValueError("E and H dump grids disagree")
        h = np.moveaxis(dump.GetFieldAtIndex(f_idx=0), 0, -1) * scale
    finite("E", e)
    finite("H", h)
    if e.shape != h.shape or e.shape[:-1] != density.shape or np.any(volumes <= 0):
        raise ValueError("invalid field/cell layout")
    selected = density > 0
    c = manifest["config"]
    if not np.allclose(sigma[selected], c["sigma_s_m"]) or not np.allclose(
        density[selected], c["density_kg_m3"]
    ):
        raise ValueError("native material dump disagrees with phantom configuration")
    bplus, bminus = circular_b(h)
    sar = local_sar(e[selected], sigma[selected], density[selected])
    mass = float(np.sum(density[selected] * volumes[selected]))
    absorbed = float(np.sum(sar * density[selected] * volumes[selected]))
    bmean = float(np.average(abs(bplus[selected]), weights=volumes[selected]))
    bstd = float(
        np.sqrt(np.average((abs(bplus[selected]) - bmean) ** 2, weights=volumes[selected]))
    )
    log = (folder / "solver.log").read_text()
    steps, decay_db = termination(log, manifest["end_criteria"])
    native_version = re.search(r"openEMS 64bit -- version (\S+)", log)
    native_engine = re.search(r"Create FDTD engine (.+)", log)
    report = {
        "schema_version": 1,
        "source": "openEMS / Cartesian FDTD",
        "native_version": native_version[1] if native_version else "unknown; see solver.log",
        "native_engine": native_engine[1] if native_engine else "unknown; see solver.log",
        "frequency_hz": f,
        "field_convention": "peak phasors; exp(+i omega t); B0 along +z",
        "s11": [s11.real, s11.imag],
        "reference_impedance_ohm": 50,
        "unscaled_incident_w": 1,
        "unscaled_accepted_w": accepted_fraction,
        "target_accepted_w": 1,
        "duty_cycle": 1,
        "spectral_normalization": "E,H / sqrt(P_acc), phase referenced to incident voltage; pulse units cancel",
        "tissue_mass_kg": mass,
        "tissue_absorbed_w": absorbed,
        "tissue_absorbed_fraction": absorbed,
        "mean_sar_w_kg": absorbed / mass,
        "max_cell_local_sar_w_kg": float(sar.max()),
        "bplus_mean_uT": bmean * 1e6,
        "bplus_coefficient_of_variation": bstd / bmean,
        "energy_decay_db": decay_db,
        "timesteps": steps,
        "mesh_cells": manifest["mesh_cells"],
        "cell_m": manifest["cell_m"],
        "power_balance_check": "FAIL" if absorbed > 1.05 else "within bound",
        "sar_averaging": "cell-interpolated E local SAR; no 1 g or 10 g averaging",
        "validation": "independent solver output; discretization and boundary convergence required",
    }
    np.savez_compressed(
        folder / "fields.npz",
        x_m=axes[0],
        y_m=axes[1],
        z_m=axes[2],
        E_V_m=e,
        H_A_m=h,
        Bplus_T=bplus,
        Bminus_T=bminus,
        volumes_m3=volumes,
        tissue_mask=selected,
        frequency_hz=f,
        target_accepted_w=1,
    )
    (folder / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def termination(log, end_criteria):
    match = re.search(
        r"end-criteria of -([\d.]+)dB reached after (\d+) timesteps \(-([\d.]+)dB\)", log
    )
    if match is None or -float(match[3]) > 10 * np.log10(end_criteria) + 0.011:
        raise ValueError("openEMS did not report reaching the requested energy-decay criterion")
    return int(match[2]), -float(match[3])


def run(
    folder,
    config,
    *,
    cell_m=0.006,
    end_criteria=1e-6,
    max_steps=2000000,
    engine="multithreaded",
    threads=4,
    postprocess=False,
):
    folder = Path(folder).resolve()
    if not np.isfinite(end_criteria) or not 0 < end_criteria < 1:
        raise ValueError("end criterion must be between 0 and 1")
    if max_steps <= 0 or threads <= 0:
        raise ValueError("maximum steps and threads must be positive")
    if folder.exists() and any(folder.iterdir()) and not postprocess:
        raise ValueError("output directory must be new or empty")
    if postprocess:
        previous = json.loads((folder / "manifest.json").read_text())
        config = Birdcage(**previous["config"])
        cell_m, end_criteria, max_steps = [
            previous[k] for k in ("cell_m", "end_criteria", "max_steps")
        ]
    fdtd, port, manifest, geometry = build(config, cell_m, end_criteria, max_steps)
    folder.mkdir(parents=True, exist_ok=True)
    if not postprocess:
        manifest["requested_engine"] = engine
        manifest["requested_threads"] = threads
        fdtd.Write2XML(str(folder / "birdcage.xml"))
        (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        axes, shield, rings, rungs = geometry
        np.savez_compressed(
            folder / "geometry.npz",
            x_m=axes[0],
            y_m=axes[1],
            z_m=axes[2],
            shield_segments=shield,
            ring_rectangles=rings,
            **{f"rung_{i}": r["segments"] for i, r in enumerate(rungs)},
        )
        (folder / "status.json").write_text('{"state": "running"}\n')
        command = [
            sys.executable,
            "-m",
            "mri_rf.openems",
            str(folder),
            "--worker",
            "--engine",
            engine,
            "--threads",
            str(threads),
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(str(Path(p).resolve()) for p in sys.path if p)
        with (folder / "solver.log").open("w") as log:
            result = subprocess.run(
                command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        if result.returncode:
            (folder / "status.json").write_text('{"state": "failed"}\n')
            raise ValueError(f"openEMS failed; see {folder / 'solver.log'}")
    try:
        report = analyze(folder, port, manifest)
    except Exception as exc:
        (folder / "status.json").write_text(
            json.dumps({"state": "incomplete", "reason": str(exc)}) + "\n"
        )
        raise
    (folder / "status.json").write_text('{"state": "done"}\n')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--cell-mm", type=float, default=6)
    parser.add_argument("--end-criteria", type=float, default=1e-6)
    parser.add_argument("--max-steps", type=int, default=2000000)
    parser.add_argument(
        "--engine", choices=("basic", "multithreaded", "gpu"), default="multithreaded"
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--postprocess", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.worker:
        folder = args.output.resolve()
        m = json.loads((folder / "manifest.json").read_text())
        fdtd, _, _, _ = build(
            Birdcage(**m["config"]), m["cell_m"], m["end_criteria"], m["max_steps"]
        )
        code = fdtd.Run(
            str(folder),
            engine=args.engine,
            numThreads=args.threads,
            dump_statistics=True,
            debug_PEC=True,
            exact_endcriteria=True,
        )
        if code not in (None, 0):
            raise RuntimeError(f"openEMS setup failed: {code}")
        return
    config = Birdcage(**json.loads(args.config.read_text())) if args.config else Birdcage()
    report = run(
        args.output,
        config,
        cell_m=args.cell_mm / 1000,
        end_criteria=args.end_criteria,
        max_steps=args.max_steps,
        engine=args.engine,
        threads=args.threads,
        postprocess=args.postprocess,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
