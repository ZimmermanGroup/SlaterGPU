import os
from pathlib import Path
import numpy as np
import pytest

@pytest.mark.parametrize("filename", ["A", "Ciap", "pVp"])
def test_output(filename):
    # convert from script location to where the output files should be located
    script_dir = Path(os.path.abspath(__file__)).parent
    example_output_dir = script_dir.parent / "build" / "examples" / "geom_1"

    assert example_output_dir.exists()
    output = np.loadtxt(example_output_dir / filename, skiprows=1)
    reference_output = np.loadtxt(example_output_dir / f"{filename}_ref", skiprows=1)
    assert output.shape == reference_output.shape
    assert output == pytest.approx(reference_output, rel=1e-6, abs=1e-6)
