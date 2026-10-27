import numpy as np
import pytest

from mri_rf.inverse import minimum_norm_operator, point_spread


def test_complex_correlated_noise_against_closed_form():
    a = np.array([[1, 1j, 0.3], [0.5j, 2, 1 - 1j]])
    c = np.array([[2, 0.3j], [-0.3j, 1]])
    w, info = minimum_norm_operator(a, c, 0.2)
    expected = a.conj().T @ np.linalg.inv(a @ a.conj().T + info["regularization"] ** 2 * c)
    np.testing.assert_allclose(w, expected, atol=1e-12)
    np.testing.assert_allclose(point_spread(w, a, 1), (w @ a)[:, 1])
    assert 0 < info["effective_dof"] < 2
    np.testing.assert_allclose(info["noise_sd"] ** 2, np.diag(w @ c @ w.conj().T).real)


def test_global_phase_and_channel_basis_invariance():
    rng = np.random.default_rng(19)
    a = rng.normal(size=(4, 9)) + 1j * rng.normal(size=(4, 9))
    q, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
    c = np.diag([1, 2, 3, 4])
    w, _ = minimum_norm_operator(a, c)
    rotated, _ = minimum_norm_operator(q @ a, q @ c @ q.conj().T)
    np.testing.assert_allclose(rotated @ q, w, atol=1e-12)


@pytest.mark.parametrize("c", [np.diag([1, 0]), np.diag([1, -1]), [[1, 1j], [1j, 1]]])
def test_reject_invalid_covariance(c):
    with pytest.raises(ValueError):
        minimum_norm_operator(np.ones((2, 3)), c)


def test_small_physical_noise_units_do_not_hide_nonhermitian_covariance():
    with pytest.raises(ValueError, match="Hermitian"):
        minimum_norm_operator(np.ones((2, 3)), np.array([[1, 0.1j], [0.1j, 1]]) * 1e-24)


def test_realistic_channel_count_with_floating_point_warnings_as_errors():
    from mri_rf.inverse import reconstruct

    rng = np.random.default_rng(9)
    a = rng.normal(size=(32, 128)) + 1j * rng.normal(size=(32, 128))
    with np.errstate(all="raise"):
        w, _ = minimum_norm_operator(a, np.eye(32))
        result = reconstruct(w, np.ones(32))
    assert np.isfinite(result).all()
