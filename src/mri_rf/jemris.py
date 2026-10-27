"""Write JEMRIS ExternalCoil maps with explicit grid and amplitude conventions."""

import json
from pathlib import Path
from xml.etree.ElementTree import Element, ElementTree, SubElement, indent

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator

from .fields import load_npz
from .physics import finite


def write_coil_map(path, relative_b1_xyz, extent_mm):
    """Input array indexed [x,y,z]; HDF5 [z,y,x] gives x-fastest JEMRIS memory.

    Coordinates are (i/N - 1/2)*Extent in mm, NOT endpoint-inclusive linspace.
    Magnitude is dimensionless relative to a separately specified RF amplitude.
    """
    path = Path(path)
    b = finite("coil map", relative_b1_xyz)
    if b.ndim != 3 or len(set(b.shape)) != 1 or b.shape[0] < 2:
        raise ValueError("JEMRIS map must be cubic with at least two points per axis")
    if not np.isfinite(extent_mm) or extent_mm <= 0:
        raise ValueError("extent must be positive millimetres")
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("maps/magnitude", data=np.abs(b).transpose(2, 1, 0).astype("float64"))
        f.create_dataset("maps/phase", data=np.angle(b).transpose(2, 1, 0).astype("float64"))
        f.attrs["axis_order"] = "zyx; x fastest"
        f.attrs["extent_mm"] = extent_mm
        f.attrs["magnitude_units"] = "dimensionless relative sensitivity"
        f.attrs["phase_units"] = "radians"
    root = Element("CoilArray")
    SubElement(
        root,
        "ExternalCoil",
        Name="RFMap",
        Filename=path.name,
        Dim="3",
        Points=str(b.shape[0]),
        Extent=str(extent_mm),
        Scale="1",
        Conj="0",
    )
    indent(root)
    ElementTree(root).write(path.with_suffix(".xml"), encoding="utf-8", xml_declaration=True)


def export_jemris(
    field_file,
    output,
    *,
    points=24,
    extent_mm=200.0,
    reference_b1_t=1e-6,
    component="plus",
    conjugate=False,
):
    if points < 2 or int(points) != points:
        raise ValueError("points must be an integer >= 2")
    if not np.isfinite(reference_b1_t) or reference_b1_t <= 0:
        raise ValueError("reference B1 must be finite and positive")
    if component not in ("plus", "minus"):
        raise ValueError("component must be plus or minus")
    if not np.isfinite(extent_mm) or extent_mm <= 0:
        raise ValueError("extent must be positive")
    data = load_npz(field_file)
    values = data["Bplus_T" if component == "plus" else "Bminus_T"]
    values = values.conj() if conjugate else values
    axis = (np.arange(points) / points - 0.5) * extent_mm * 1e-3
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    interpolate = LinearNDInterpolator(data["points_m"], values, fill_value=np.nan)
    sampled = interpolate(grid)
    mask = np.isfinite(sampled)
    if not mask.any():
        raise ValueError("grid does not overlap the field sample convex hull")
    # Interpolate complex values before extracting phase, avoiding branch-cut artifacts.
    sampled = np.where(mask, sampled / reference_b1_t, 0)
    write_coil_map(output, sampled, extent_mm)
    np.savez_compressed(
        Path(output).with_suffix(".grid.npz"), axis_m=axis, valid=mask, relative_b1=sampled
    )
    metadata = {
        "reference_b1_t": reference_b1_t,
        "component": component,
        "conjugated": conjugate,
        "points": points,
        "extent_mm": extent_mm,
        "valid_fraction": float(mask.mean()),
        "sampling": "linear interpolation of tetrahedron-centre fields; zero outside convex hull",
        "warning": "convex-hull support is not an anatomical mask; validate with JEMRIS samples",
    }
    Path(output).with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
