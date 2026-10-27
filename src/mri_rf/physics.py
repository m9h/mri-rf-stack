"""SI units, peak complex phasors with exp(+i omega t), B0 along +z.

The circular components are named by algebra, not by a universal transmit/receive
handedness. Receiver conjugation belongs in the acquisition convention.
"""

from functools import lru_cache
from itertools import product

import numpy as np
from numpy.polynomial.legendre import leggauss

MU0 = 1.25663706212e-6
PROTON_HZ_PER_T = 42.577478518e6
EDGES = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))


def finite(name, value):
    a = np.asarray(value)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain only finite values")
    return a


def circular_b(h):
    """Return mu0*(Hx+iHy)/2 and mu0*(Hx-iHy)/2 in tesla (mu_r=1)."""
    h = finite("H", h)
    if h.shape[-1:] != (3,):
        raise ValueError("H must have final axis of length 3")
    return MU0 * (h[..., 0] + 1j * h[..., 1]) / 2, MU0 * (h[..., 0] - 1j * h[..., 1]) / 2


def port_power(s, a):
    """Incident and net accepted W for real-reference power-wave ports.

    Each a is a complex incident amplitude in sqrt(W); b=S@a includes coupling.
    This is NOT sum(|a|^2 * (1-|Sii|^2)), which drops cross terms.
    """
    s, a = finite("S", s), finite("drive", a)
    if a.ndim != 1 or s.shape != (a.size, a.size) or a.size == 0:
        raise ValueError("S must be square and match the nonempty drive vector")
    inc = float(np.vdot(a, a).real)
    if inc <= 0:
        raise ValueError("incident power must be positive")
    if np.linalg.svd(s, compute_uv=False)[0] > 1 + 1e-5:
        raise ValueError("S is not passive; check normalization or solver convergence")
    ref = float(np.vdot(s @ a, s @ a).real)
    accepted = inc - ref
    if accepted <= inc * 1e-12:
        raise ValueError("accepted power is zero or negative")
    return inc, accepted


def local_sar(e, sigma, rho):
    """Cycle-averaged local SAR in W/kg for peak E in V/m."""
    e, sigma, rho = finite("E", e), finite("sigma", sigma), finite("rho", rho)
    if e.shape[-1:] != (3,) or np.any(sigma < 0) or np.any(rho <= 0):
        raise ValueError("expected vector E, nonnegative conductivity and positive density")
    return sigma * np.sum(np.abs(e) ** 2, axis=-1) / (2 * rho)


def tet_volumes(nodes, tets):
    p = finite("nodes", nodes)[tets]
    v = np.abs(np.linalg.det(p[:, 1:] - p[:, :1])) / 6
    if np.any(v <= 0):
        raise ValueError("degenerate tetrahedron")
    return v


def shape_functions(bary, npt):
    if npt == 4:
        return bary
    if npt != 10:
        raise ValueError("only tet4 and VTK-order tet10 samples are supported")
    return np.concatenate(
        [
            bary * (2 * bary - 1),
            np.stack([4 * bary[..., i] * bary[..., j] for i, j in EDGES], axis=-1),
        ],
        axis=-1,
    )


@lru_cache(maxsize=2)
def field_mass_matrix(npt):
    """Exact integration of squared linear/quadratic interpolants on straight tets.

    4-point Gauss in each Duffy coordinate integrates the degree-4 polynomial
    including its Jacobian; weights sum to one (element mean, not integral).
    """
    x, w = leggauss(4)
    x, w = (x + 1) / 2, w / 2
    bary, weights = [], []
    for i, j, k in product(range(4), repeat=3):
        u, v, z = x[i], x[j], x[k]
        bary.append([1 - u, u * (1 - v), u * v * (1 - z), u * v * z])
        weights.append(6 * w[i] * w[j] * w[k] * u**2 * v)
    phi = shape_functions(np.asarray(bary), npt)
    return np.einsum("qi,q,qj->ij", phi, np.asarray(weights), phi)


def mean_square_field(samples):
    """Element mean |E|^2 from samples shaped (element, 4 or 10, xyz)."""
    s = finite("field samples", samples)
    if s.ndim != 3 or s.shape[-1] != 3:
        raise ValueError("expected (element, 4 or 10, 3) samples")
    q = field_mass_matrix(s.shape[1])
    result = np.einsum("tic,ij,tjc->t", s.conj(), q, s).real
    return np.maximum(result, 0)
