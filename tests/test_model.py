import pytest

from mri_rf.model import Birdcage


def test_phantom_must_fit_and_gaps_cannot_overlap():
    with pytest.raises(ValueError):
        Birdcage(phantom_radius_m=0.2).validate()
    with pytest.raises(ValueError):
        Birdcage(gap_m=0.05).validate()


def test_native_hertzon_validation(tmp_path):
    pytest.importorskip("hertzon")
    from hertzon.model.design import Project
    from hertzon.model.validate import validate

    from mri_rf.model import write_birdcage

    path = write_birdcage(tmp_path, Birdcage())
    project = Project.load(path)
    assert not [c for c in validate(project, project.designs[0]) if c.level == "error"]
