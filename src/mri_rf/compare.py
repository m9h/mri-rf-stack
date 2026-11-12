"""Compare power-normalized FEM and FDTD results on shared interior points."""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .fields import analyze_hertzon, load_npz
from .model import Birdcage
from .physics import circular_b, finite, shape_functions


def sample_tetrahedra(points, tetrahedra, samples):
    """Evaluate each tetrahedron's own interpolant; do not blend across elements.

    A point on a shared face uses the first containing element. Points outside
    all tetrahedra raise rather than being silently extrapolated or set to zero.
    """
    points, tetrahedra, samples = [
        finite(n, v)
        for n, v in (("points", points), ("tetrahedra", tetrahedra), ("samples", samples))
    ]
    if (
        points.ndim != 2
        or points.shape[1] != 3
        or tetrahedra.ndim != 3
        or tetrahedra.shape[1:] != (4, 3)
        or samples.ndim != 3
        or len(samples) != len(tetrahedra)
        or samples.shape[1] not in (4, 10)
    ):
        raise ValueError("incompatible tetrahedral interpolation arrays")
    inv = np.linalg.inv(np.swapaxes(tetrahedra[:, 1:] - tetrahedra[:, :1], 1, 2))
    output = np.empty((len(points), samples.shape[-1]), dtype=np.result_type(samples, float))
    for start in range(0, len(points), 128):
        q = points[start : start + 128]
        xyz = np.einsum("tij,qtj->qti", inv, q[:, None] - tetrahedra[None, :, 0])
        bary = np.concatenate([1 - xyz.sum(axis=-1, keepdims=True), xyz], axis=-1)
        inside = np.all((bary >= -1e-9) & (bary <= 1 + 1e-9), axis=-1)
        if not np.all(inside.any(axis=1)):
            raise ValueError("comparison points outside FEM tissue support")
        owners = inside.argmax(axis=1)
        phi = shape_functions(bary[np.arange(len(q)), owners], samples.shape[1])
        output[start : start + len(q)] = np.einsum("qi,qic->qc", phi, samples[owners])
    return output


def field_metrics(reference, candidate):
    """Magnitude error and complex similarity, without an amplitude fit."""
    a, b = finite("reference", reference).ravel(), finite("candidate", candidate).ravel()
    if a.shape != b.shape or not a.size or np.linalg.norm(a) == 0 or np.linalg.norm(b) == 0:
        raise ValueError("field comparison needs matching, nonzero arrays")
    cross = np.vdot(b, a)
    phase = cross / abs(cross) if abs(cross) else 1 + 0j
    return {
        "magnitude_nrmse": float(np.linalg.norm(abs(b) - abs(a)) / np.linalg.norm(a)),
        "complex_correlation": float(abs(cross) / (np.linalg.norm(a) * np.linalg.norm(b))),
        "global_phase_alignment_rad": float(np.angle(phase)),
        "phase_aligned_complex_nrmse": float(np.linalg.norm(phase * b - a) / np.linalg.norm(a)),
        "reference_mean_magnitude": float(np.mean(abs(a))),
        "candidate_mean_magnitude": float(np.mean(abs(b))),
    }


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(hertzon, config_path, openems, output, *, inset_m=0.015, spacing_m=0.005):
    hertzon, output = Path(hertzon), Path(output)
    config = Birdcage(**json.loads(Path(config_path).read_text()))
    config.validate()
    if not 0 < inset_m < min(config.phantom_radius_m, config.phantom_length_m / 2):
        raise ValueError("interior inset must be positive and leave a nonempty phantom")
    if not np.isfinite(spacing_m) or spacing_m <= 0:
        raise ValueError("sampling spacing must be positive")
    output.mkdir(parents=True, exist_ok=True)
    reference = analyze_hertzon(
        hertzon,
        {"Phantom": {"sigma_s_m": config.sigma_s_m, "density_kg_m3": config.density_kg_m3}},
        output / "hertzon-analysis",
    )
    mesh = load_npz(hertzon / "mesh.npz")
    fields = load_npz(hertzon / "fields_LastAdaptive.npz")
    params = load_npz(hertzon / "sparams_LastAdaptive.npz")
    if params["S"].shape != (1, 1, 1) or not np.allclose(params["z0"], 50):
        raise ValueError(
            "comparison currently requires a single-frequency, 50 ohm, single-port run"
        )
    s11 = complex(params["S"][0, 0, 0])
    reference["s11"] = [s11.real, s11.imag]
    selected = mesh["tet_obj"] == mesh["obj_names"].tolist().index("Phantom")
    tetrahedra = mesh["nodes"][mesh["tets"][selected]]
    npt = 10 if "tet10" in fields else 4
    scale = np.sqrt(1 / reference["unscaled_accepted_w"])
    es = fields["E_0_1"].reshape(-1, npt, 3)[selected] * scale
    hs = fields["H_0_1"].reshape(-1, npt, 3)[selected] * scale
    radius = config.phantom_radius_m - inset_m
    half_length = config.phantom_length_m / 2 - inset_m
    xy = (
        np.arange(-int(np.floor(radius / spacing_m)), int(np.floor(radius / spacing_m)) + 1)
        * spacing_m
    )
    z = (
        np.arange(
            -int(np.floor(half_length / spacing_m)), int(np.floor(half_length / spacing_m)) + 1
        )
        * spacing_m
    )
    points = np.stack(np.meshgrid(xy, xy, z, indexing="ij"), axis=-1).reshape(-1, 3)
    points = points[np.linalg.norm(points[:, :2], axis=1) <= radius + 1e-12]
    reference_samples = sample_tetrahedra(points, tetrahedra, np.concatenate([es, hs], axis=-1))
    e_ref, h_ref = reference_samples[:, :3], reference_samples[:, 3:]
    b_ref, _ = circular_b(h_ref)
    saved = dict(points_m=points, hertzon_E_V_m=e_ref, hertzon_H_A_m=h_ref, hertzon_Bplus_T=b_ref)
    comparisons = []
    input_hashes = {
        f"hertzon/{name}": digest(hertzon / name)
        for name in ("mesh.npz", "fields_LastAdaptive.npz", "sparams_LastAdaptive.npz")
    }
    for index, folder in enumerate(map(Path, openems)):
        manifest = json.loads((folder / "manifest.json").read_text())
        report = json.loads((folder / "report.json").read_text())
        if json.loads((folder / "status.json").read_text()).get("state") != "done":
            raise ValueError(f"incomplete openEMS run: {folder}")
        physical_keys = set(manifest["config"]) - {"mesh_m", "passes"}
        if any(manifest["config"][k] != getattr(config, k) for k in physical_keys):
            raise ValueError("FEM and FDTD physical input configurations differ")
        if not np.isclose(report["frequency_hz"], reference["frequency_hz"], rtol=1e-10):
            raise ValueError("FEM and FDTD frequencies differ")
        if report["target_accepted_w"] != 1 or report["duty_cycle"] != 1:
            raise ValueError("comparison expects 1 W accepted and duty cycle 1")
        data = load_npz(folder / "fields.npz")
        axes = [data[f"{a}_m"] for a in "xyz"]
        e = RegularGridInterpolator(axes, data["E_V_m"], bounds_error=True)(points)
        h = RegularGridInterpolator(axes, data["H_A_m"], bounds_error=True)(points)
        bplus, _ = circular_b(h)
        bulk = {}
        for key in (
            "unscaled_accepted_w",
            "tissue_mass_kg",
            "tissue_absorbed_w",
            "mean_sar_w_kg",
            "bplus_mean_uT",
        ):
            bulk[key] = {
                "hertzon": reference[key],
                "openems": report[key],
                "relative_difference": report[key] / reference[key] - 1,
            }
        comparisons.append(
            {
                "name": folder.name,
                "openems_report": report,
                "bulk": bulk,
                "absolute_complex_s11_difference": float(abs(complex(*report["s11"]) - s11)),
                "interior_bplus": field_metrics(b_ref, bplus),
                "interior_e_vector": field_metrics(e_ref, e),
                "interior_h_vector": field_metrics(h_ref, h),
            }
        )
        saved.update(
            {
                f"openems_{index}_E_V_m": e,
                f"openems_{index}_H_A_m": h,
                f"openems_{index}_Bplus_T": bplus,
            }
        )
        input_hashes.update(
            {
                f"{folder.name}/{name}": digest(folder / name)
                for name in ("fields.npz", "report.json", "manifest.json")
            }
        )
    result = {
        "schema_version": 1,
        "physical_config": asdict(config),
        "hertzon_report": reference,
        "comparisons": comparisons,
        "common_interior": {"inset_m": inset_m, "spacing_m": spacing_m, "points": len(points)},
        "interpolation": "FEM tet4/tet10 shape functions; FDTD trilinear complex cell-centre fields",
        "normalization": "1 W net accepted at 50 ohm feed; no amplitude fitting",
        "peak_sar_comparison": "omitted: FEM element means and FDTD cell-interpolated local SAR are different estimators",
        "input_sha256": input_hashes,
        "interpretation": "numerical cross-check; agreement and physical validation are not assumed",
    }
    (output / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "common-grid.npz", **saved)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hertzon", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--openems", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = compare(args.hertzon, args.config, args.openems, args.output)
    for item in result["comparisons"]:
        print(
            item["name"],
            json.dumps({"bulk": item["bulk"], "interior_bplus": item["interior_bplus"]}, indent=2),
        )


if __name__ == "__main__":
    main()
