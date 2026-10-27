from xml.etree.ElementTree import parse

import h5py
import numpy as np
import pytest

from mri_rf.jemris import write_coil_map


def test_hdf5_matches_cpp_x_fastest_memory_and_radian_phase(tmp_path):
    x, y, z = np.indices((4, 4, 4))
    mag = 1 + x + 10 * y + 100 * z
    phase = 0.1 * x - 0.2 * y + 0.05 * z
    b = mag * np.exp(1j * phase)
    path = tmp_path / "coil.h5"
    write_coil_map(path, b, 200)
    with h5py.File(path) as f:
        flat = f["maps/magnitude"][:].ravel()
        for i, j, k in ((1, 2, 3), (3, 1, 0), (0, 0, 2)):
            # JEMRIS NDData operator()(i,j,k): i+N*j+N*N*k.
            assert flat[i + 4 * j + 16 * k] == pytest.approx(mag[i, j, k])
        np.testing.assert_allclose(f["maps/phase"][:], phase.transpose(2, 1, 0), atol=1e-14)
    element = parse(path.with_suffix(".xml")).getroot()[0]
    assert element.attrib["Filename"] == "coil.h5"
    assert element.attrib["Points"] == "4"
    assert element.attrib["Extent"] == "200"
