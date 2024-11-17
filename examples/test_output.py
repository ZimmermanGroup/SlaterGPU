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

    with open(example_output_dir / filename) as output, open(
        example_output_dir / f"{filename}_ref"
    ) as reference_output:
        output = list(output)
        for i, line in enumerate(reference_output):
            if line.endswith(":"):
                continue
            assert np.fromstring(output[i], sep=" ") == pytest.approx(
                np.fromstring(line, sep=" "), rel=1e-6, abs=1e-6
            )
        # matrix_boundaries = [line for line in file if line.endswith(':')] + [-1]
    # for idx in range(len(matrix_boundaries)-1):

    # output = np.loadtxt(example_output_dir / filename, skiprows=matrix_boundaries[idx])
    # reference_output = np.loadtxt(example_output_dir / f"{filename}_ref", skiprows=matrix_boundaries[idx]+1, max_rows=matrix_boundaries[idx+1]-matrix_boundaries[idx])
    # assert output.shape == reference_output.shape
    # assert output == pytest.approx(reference_output, rel=1e-6, abs=1e-6)
