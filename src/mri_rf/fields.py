"""Import Hertzon 0.1.25 straight-tetrahedron results without importing its GUI."""

import json
from pathlib import Path

import numpy as np

from .physics import circular_b, finite, mean_square_field, port_power, shape_functions, tet_volumes


def load_npz(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def analyze_hertzon(
    folder,
    materials,
    output,
    *,
    label="LastAdaptive",
    frequency_index=0,
    drive=None,
    target_accepted_w=1.0,
    duty_cycle=1.0,
):
    """Normalize coherent E/H to target net accepted power; report local/element SAR.

    Material keys are exact Hertzon volume names; each listed volume must have sigma
    and density. Unlisted volumes are explicitly excluded from the tissue analysis.
    """
    folder, output = Path(folder), Path(output)
    if not np.isfinite(target_accepted_w) or target_accepted_w <= 0:
        raise ValueError("target accepted power must be positive")
    if not np.isfinite(duty_cycle) or not 0 <= duty_cycle <= 1:
        raise ValueError("duty cycle must be between 0 and 1")
    if frequency_index < 0 or int(frequency_index) != frequency_index:
        raise ValueError("frequency index must be a nonnegative integer")
    if (
        not label
        or any(
            ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-"
            for ch in label
        )
        or label in (".", "..")
    ):
        raise ValueError("invalid solution label")
    status = json.loads((folder / "status.json").read_text())
    if status.get("state") != "done":
        raise ValueError("solver has not completed successfully")
    mesh = load_npz(folder / "mesh.npz")
    fields = load_npz(folder / f"fields_{label}.npz")
    params = load_npz(folder / f"sparams_{label}.npz")
    if "curv_tets" in mesh and mesh["curv_tets"].size:
        raise ValueError("curved tetrahedra are unsupported; regenerate with CurvedElements=False")
    solpath = folder / f"solution_{label}.npz"
    if solpath.exists():
        sol = load_npz(solpath)
        if "stamp" not in fields or not np.array_equal(sol.get("stamp"), fields["stamp"]):
            raise ValueError("stale field cache: recompute fields in Hertzon")
    npt = 10 if "tet10" in fields else 4 if "tet4" in fields else 0
    if npt == 0:
        raise ValueError("unrecognized Hertzon field sampling")
    freq = float(fields[f"freq_{frequency_index}"][0])
    indices = np.flatnonzero(np.isclose(params["freq"], freq, rtol=1e-9, atol=0.01))
    if len(indices) != 1:
        raise ValueError("field frequency does not uniquely match S-parameters")
    ports = json.loads(str(params["ports"]))
    if any(p.get("kind") != "Lumped Port" for p in ports):
        raise ValueError("this adapter currently verifies only 1 W real-impedance lumped ports")
    z0 = finite("reference impedance", params["z0"])
    if np.any(np.abs(z0.imag) > 1e-12) or np.any(z0.real <= 0):
        raise ValueError("complex/nonpositive port reference impedances are unsupported")
    a = np.ones(len(ports), complex) if drive is None else finite("drive", drive).astype(complex)
    s = params["S"][indices[0]]
    incident, accepted = port_power(s, a)
    scale = np.sqrt(target_accepted_w / accepted)
    nodes, tets = mesh["nodes"], mesh["tets"]
    if nodes.ndim != 2 or nodes.shape[1] != 3 or tets.ndim != 2 or tets.shape[1] != 4:
        raise ValueError("invalid mesh dimensions")
    if np.any(tets < 0) or np.any(tets >= len(nodes)):
        raise ValueError("invalid mesh connectivity")
    shape = (len(tets) * npt, 3)
    e, h = np.zeros(shape, complex), np.zeros(shape, complex)
    for weight, port in zip(a, ports, strict=True):
        for prefix, dest in (("E", e), ("H", h)):
            array = finite(prefix, fields[f"{prefix}_{frequency_index}_{port['number']}"])
            if array.shape != shape:
                raise ValueError("field and mesh dimensions disagree")
            dest += weight * scale * array
    owners, names = mesh["tet_obj"], mesh["obj_names"].tolist()
    if owners.shape != (len(tets),):
        raise ValueError("invalid mesh ownership")
    if not materials or set(materials) - set(names):
        raise ValueError("material mapping is empty or contains unknown volume names")
    sigma, density, selected = np.zeros(len(tets)), np.ones(len(tets)), np.zeros(len(tets), bool)
    for name, props in materials.items():
        mask = owners == names.index(name)
        if not mask.any():
            raise ValueError(f"no tetrahedra owned by {name}")
        cond, rho = float(props["sigma_s_m"]), float(props["density_kg_m3"])
        if not np.isfinite(cond) or not np.isfinite(rho) or cond < 0 or rho <= 0:
            raise ValueError("materials need nonnegative conductivity and positive density")
        sigma[mask], density[mask], selected[mask] = cond, rho, True
    tet_nodes = nodes[tets[selected]]
    volumes = tet_volumes(nodes, tets)[selected]
    es = e.reshape(-1, npt, 3)[selected]
    hs = h.reshape(-1, npt, 3)[selected]
    centroid_phi = shape_functions(np.full((1, 4), 0.25), npt)[0]
    ec, hc = np.einsum("i,tic->tc", centroid_phi, es), np.einsum("i,tic->tc", centroid_phi, hs)
    bplus, bminus = circular_b(hc)
    power_density = sigma[selected] * mean_square_field(es) / 2
    sar = power_density / density[selected]
    absorbed = float(np.sum(power_density * volumes))
    mass = float(np.sum(density[selected] * volumes))
    bmean = float(np.average(np.abs(bplus), weights=volumes))
    bstd = float(np.sqrt(np.average((np.abs(bplus) - bmean) ** 2, weights=volumes)))
    convergence_path = folder / "convergence.json"
    convergence = json.loads(convergence_path.read_text()) if convergence_path.exists() else None
    report = {
        "schema_version": 1,
        "source": "Hertzon 0.1.25 / EMerge 3.0.0a25",
        "frequency_hz": freq,
        "field_convention": "peak phasors; exp(+i omega t); B0 along +z",
        "drive_sqrt_w": [[float(x.real), float(x.imag)] for x in a],
        "unscaled_incident_w": incident,
        "unscaled_accepted_w": accepted,
        "target_accepted_w": target_accepted_w,
        "duty_cycle": duty_cycle,
        "selected_volumes": list(materials),
        "tissue_mass_kg": mass,
        "tissue_absorbed_w": absorbed,
        "tissue_absorbed_fraction": absorbed / target_accepted_w,
        "mean_sar_w_kg": absorbed / mass * duty_cycle,
        "max_element_mean_sar_w_kg": float(sar.max()) * duty_cycle,
        "bplus_mean_uT": bmean * 1e6,
        "bplus_coefficient_of_variation": bstd / bmean if bmean else None,
        "tetrahedra_total": len(tets),
        "tetrahedra_tissue": int(selected.sum()),
        "convergence": convergence,
        "validation": "research output; field/mesh convergence and independent solver checks required",
        "power_balance_check": "FAIL" if absorbed > 1.05 * target_accepted_w else "within bound",
        "sar_averaging": "tetrahedron mean only; no 1 g or 10 g averaging",
    }
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "fields.npz",
        points_m=tet_nodes.mean(axis=1),
        tetrahedra_m=tet_nodes,
        volumes_m3=volumes,
        E_V_m=ec,
        H_A_m=hc,
        Bplus_T=bplus,
        Bminus_T=bminus,
        sar_w_kg=sar * duty_cycle,
        frequency_hz=freq,
        target_accepted_w=target_accepted_w,
    )
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_vtk(output / "fields.vtk", tet_nodes, sar * duty_cycle, bplus, bminus)
    return report


def write_vtk(path, tetrahedra, sar, bplus, bminus):
    """Legacy VTK unstructured grid with per-element quantities; no averaging across tissues."""
    n = len(tetrahedra)
    with Path(path).open("w") as f:
        f.write(
            f"# vtk DataFile Version 3.0\nMRI RF research fields\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS {4 * n} double\n"
        )
        np.savetxt(f, tetrahedra.reshape(-1, 3))
        f.write(f"CELLS {n} {5 * n}\n")
        np.savetxt(f, np.column_stack([np.full(n, 4), np.arange(4 * n).reshape(n, 4)]), fmt="%d")
        f.write(f"CELL_TYPES {n}\n")
        np.savetxt(f, np.full(n, 10), fmt="%d")
        f.write(f"CELL_DATA {n}\n")
        for name, values in (
            ("SAR_W_kg", sar),
            ("Bplus_uT", np.abs(bplus) * 1e6),
            ("Bminus_uT", np.abs(bminus) * 1e6),
            ("Bplus_phase_rad", np.angle(bplus)),
        ):
            f.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            np.savetxt(f, values)
