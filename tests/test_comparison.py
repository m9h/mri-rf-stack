import numpy as np
import pytest

from mri_rf.compare import field_metrics, sample_tetrahedra
from mri_rf.model import Birdcage
from mri_rf.openems import discretize, termination


@pytest.mark.parametrize("cell", [0.006, 0.004])
def test_staircase_preserves_strip_width_and_capacitance(cell):
    c = Birdcage()
    axes, shield, _, rungs = discretize(c, cell)
    assert all(np.all(np.diff(a) > 0) for a in axes)
    # A closed shield contour: each vertex must belong to two faces.
    _, counts = np.unique(shield.reshape(-1, 2), axis=0, return_counts=True)
    np.testing.assert_array_equal(counts, 2)
    for rung in rungs:
        assert rung["projected_width_m"] == pytest.approx(c.strip_width_m)
        assert rung["max_plane_error_m"] <= cell
        assert np.all(rung["weights"] > 0)
        assert (rung["weights"] * c.capacitance_f).sum() == pytest.approx(
            c.capacitance_f, abs=1e-23
        )
        assert len(np.unique(rung["nodes"], axis=0)) == len(rung["nodes"])


def test_tetrahedral_interpolation_preserves_quadratic_complex_fields():
    tet = np.array([[[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]]], float)
    from mri_rf.physics import EDGES

    nodes = np.concatenate([tet[0], [(tet[0, i] + tet[0, j]) / 2 for i, j in EDGES]])

    def field(p):
        x, y, z = p.T
        return np.column_stack([x * x + 2j * y * z, (1 + 1j) * z + x * y])

    points = np.array([[0.1, 0.2, 0.3], [0.5, 0.4, 1], [0, 0, 0]])
    actual = sample_tetrahedra(points, tet, field(nodes)[None])
    np.testing.assert_allclose(actual, field(points), atol=1e-14)
    with pytest.raises(ValueError, match="outside"):
        sample_tetrahedra(np.array([[2, 3, 4]]), tet, field(nodes)[None])


def test_comparison_does_not_hide_amplitude_disagreement():
    a = np.array([1 + 1j, 2 - 3j, 4])
    result = field_metrics(a, 2j * a)
    assert result["complex_correlation"] == pytest.approx(1)
    assert result["magnitude_nrmse"] == pytest.approx(1)
    assert result["phase_aligned_complex_nrmse"] == pytest.approx(1)
    assert result["global_phase_alignment_rad"] == pytest.approx(-np.pi / 2)


def test_time_truncated_solver_runs_are_rejected():
    good = "RunFDTD: end-criteria of -60.00dB reached after 643795 timesteps (-60.10dB)"
    assert termination(good, 1e-6) == (643795, -60.10)
    for text in (
        "RunFDTD: simulation aborted after 6000 timesteps (-2.88dB)",
        "RunFDTD: Warning: Max. number of timesteps was reached before the end-criteria",
        good.replace("60.10", "40.10"),
    ):
        with pytest.raises(ValueError, match="energy-decay"):
            termination(text, 1e-6)
