"""Noise-whitened complex minimum-norm inverse and resolution diagnostics.

This is an InI-inspired linear inverse building block, not a reproduction of a
published pulse sequence or a claim to measure neural currents with RF coils.
"""

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.linalg.blas import zgemm, zgemv

from .physics import finite


def minimum_norm_operator(forward, covariance, alpha=0.1):
    """W for x_hat=W@y; lambda=alpha*largest singular value of whitened A.

    forward: (channels, sources), including any acquisition encoding/voxel weights.
    covariance: E[n n^H], positive-definite, same channel order and units as y.
    The source prior is identity. alpha must be positive.
    """
    a = finite("forward", forward).astype(complex)
    c = finite("covariance", covariance).astype(complex)
    if a.ndim != 2 or min(a.shape) == 0 or c.shape != (a.shape[0], a.shape[0]):
        raise ValueError("expected A=(channels,sources), C=(channels,channels)")
    if np.max(np.abs(c - c.conj().T)) > 1e-10 * np.max(np.abs(c)):
        raise ValueError("noise covariance must be Hermitian")
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive")
    try:
        lower = cholesky(c, lower=True)
    except np.linalg.LinAlgError as exc:
        raise ValueError("noise covariance must be positive definite") from exc
    whitening = solve_triangular(lower, np.eye(len(c)), lower=True)
    # SciPy's BLAS also avoids spurious floating-point exception flags from
    # NumPy 2.2 / Apple Accelerate complex matmul on current macOS.
    aw = zgemm(1, whitening, a)
    u, s, vh = np.linalg.svd(aw, full_matrices=False)
    if s[0] == 0:
        raise ValueError("forward model has no sensitivity")
    lam = alpha * s[0]
    filt = s / (s**2 + lam**2)
    operator = zgemm(1, zgemm(1, vh.conj().T * filt, u.conj().T), whitening)
    noise_sd = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", operator, c, operator.conj()).real, 0))
    return operator, {
        "singular_values": s,
        "regularization": lam,
        "effective_dof": np.sum(s**2 / (s**2 + lam**2)),
        "noise_sd": noise_sd,
    }


def point_spread(operator, forward, source):
    """One column of W A, avoiding allocation of the full source resolution matrix."""
    if not 0 <= source < forward.shape[1]:
        raise ValueError("source index out of range")
    return reconstruct(operator, forward[:, source])


def reconstruct(operator, data):
    """Apply a complex inverse to one observation or a channel-by-time matrix."""
    operator, data = finite("operator", operator), finite("data", data)
    if operator.ndim != 2 or data.ndim not in (1, 2) or data.shape[0] != operator.shape[1]:
        raise ValueError("data must have the same channel count as the operator")
    result = zgemv(1, operator, data) if data.ndim == 1 else zgemm(1, operator, data)
    return finite("reconstruction", result)
