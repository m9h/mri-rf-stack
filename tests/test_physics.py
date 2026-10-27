import numpy as np
import pytest

from mri_rf.physics import (
    EDGES,
    MU0,
    circular_b,
    local_sar,
    mean_square_field,
    port_power,
    tet_volumes,
)


def test_circular_polarization_and_peak_sar():
    plus, minus = circular_b(np.array([[1, -1j, 0], [1, 1j, 0]]))
    np.testing.assert_allclose(plus, [MU0, 0])
    np.testing.assert_allclose(minus, [0, MU0])
    assert local_sar(np.array([10, 0, 0]), 0.5, 1000) == 0.025


def test_coherent_power_includes_coupling():
    s = np.array([[0.2, 0.3j], [0.3j, 0.1]])
    a = np.array([1, 1j])
    incident, accepted = port_power(s, a)
    assert incident == 2
    assert accepted == pytest.approx(2 - 0.1**2 - 0.4**2)
    _, wrong = port_power(np.diag(np.diag(s)), a)
    assert abs(accepted - wrong) > 0.1


@pytest.mark.parametrize("s,a", [([[1.2]], [1]), ([[1]], [1]), ([[0]], [0]), ([[0]], [1, 1])])
def test_reject_invalid_port_normalization(s, a):
    with pytest.raises(ValueError):
        port_power(np.array(s), np.array(a))


def test_tetrahedron_quadrature_against_analytic_integrals():
    vertices = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert tet_volumes(vertices, [[0, 1, 2, 3]])[0] == pytest.approx(1 / 6)
    # E=(x,y,z): <x^2+y^2+z^2>=3/10 in the unit simplex.
    assert mean_square_field(vertices[None])[0] == pytest.approx(3 / 10)
    midpoints = np.array([(vertices[i] + vertices[j]) / 2 for i, j in EDGES])
    points = np.concatenate([vertices, midpoints])
    # E=(x^2,y^2,z^2): exact degree-four integral is 3/35.
    assert mean_square_field((points**2)[None])[0] == pytest.approx(3 / 35)
    # Uniform field must remain uniform regardless of nodal representation.
    assert mean_square_field(np.full((1, 10, 3), 2 + 3j))[0] == pytest.approx(39)
