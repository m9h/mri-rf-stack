import json

import numpy as np
import pytest

from mri_rf.fields import analyze_hertzon


def fixture(folder):
    (folder / "status.json").write_text('{"state":"done"}')
    nodes = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    np.savez(
        folder / "mesh.npz", nodes=nodes, tets=[[0, 1, 2, 3]], tet_obj=[0], obj_names=["Phantom"]
    )
    e = np.tile([1.0, 0, 0], (4, 1))
    h = np.tile([1.0, -1j, 0], (4, 1))
    np.savez(folder / "fields_LastAdaptive.npz", tet4=[1], freq_0=[128e6], E_0_1=e, H_0_1=h)
    np.savez(
        folder / "sparams_LastAdaptive.npz",
        freq=[128e6],
        S=[[[0.5]]],
        z0=[50.0],
        ports=json.dumps([{"kind": "Lumped Port", "number": 1}]),
    )
    return {"Phantom": {"sigma_s_m": 0.5, "density_kg_m3": 1000.0}}


def test_power_and_duty_scaling_preserve_units(tmp_path):
    materials = fixture(tmp_path)
    a = analyze_hertzon(tmp_path, materials, tmp_path / "a", target_accepted_w=1, duty_cycle=1)
    b = analyze_hertzon(tmp_path, materials, tmp_path / "b", target_accepted_w=4, duty_cycle=0.25)
    assert a["unscaled_accepted_w"] == 0.75
    assert a["tissue_mass_kg"] == pytest.approx(1000 / 6)
    assert a["tissue_absorbed_w"] == pytest.approx(0.5 / 2 / 6 / 0.75)
    assert b["bplus_mean_uT"] == pytest.approx(2 * a["bplus_mean_uT"])
    assert b["mean_sar_w_kg"] == pytest.approx(a["mean_sar_w_kg"])


def test_reject_incomplete_unknown_material_and_curved_mesh(tmp_path):
    materials = fixture(tmp_path)
    with pytest.raises(ValueError, match="unknown"):
        analyze_hertzon(tmp_path, {"missing": materials["Phantom"]}, tmp_path / "a")
    (tmp_path / "status.json").write_text('{"state":"failed"}')
    with pytest.raises(ValueError, match="completed"):
        analyze_hertzon(tmp_path, materials, tmp_path / "b")
    fixture(tmp_path)
    with np.load(tmp_path / "mesh.npz") as f:
        mesh = dict(f)
    np.savez(tmp_path / "mesh.npz", **mesh, curv_tets=[0])
    with pytest.raises(ValueError, match="curved"):
        analyze_hertzon(tmp_path, materials, tmp_path / "c")
