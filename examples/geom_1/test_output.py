import os
from pathlib import Path
import numpy as np
import pytest


def test_output():
    # convert from script location to corresponding location in build directory
    script_dir = Path(os.path.abspath(__file__)).parent
    script_dir_parts = list(script_dir.parts)
    script_dir_parts.insert(-2, "build")
    script_build_dir = Path(*script_dir_parts)

    assert script_build_dir.exists()
    A = np.loadtxt(script_build_dir / "A", skiprows=1)
    A_ref = np.loadtxt(script_build_dir / "A_ref", skiprows=1)
    assert A.shape == A_ref.shape
    assert A == pytest.approx(A_ref, rel=1e-6, abs=1e-6)
