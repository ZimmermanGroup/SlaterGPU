from itertools import pairwise
import os
from pathlib import Path
import numpy as np
import pytest

@pytest.mark.parametrize("filename", ["A", "Ciap", "pVp"])
def test_output(filename):
    # Convert from the script location to where the output files should be located.
    script_dir = Path(os.path.abspath(__file__)).parent
    example_output_dir = script_dir.parent / "build" / "examples" / "geom_1"

    assert example_output_dir.exists()

    with open(example_output_dir / f"{filename}_ref") as reference_output:
        # output = list(output)
        # for i, line in enumerate(reference_output):
        # if line.endswith(":"):
        #     continue
        # assert np.fromstring(output[i], sep=" ") == pytest.approx(
        #     np.fromstring(line, sep=" "), rel=1e-6, abs=1e-6
        # )
        matrix_boundaries = [
            line_num
            for line_num, line in enumerate(reference_output)
            # Get header lines for matrices and blank last line.
            if line.strip().endswith(":")
        ] + [len(reference_output)]  # JOSH - RESUME HERE
    print(f"{matrix_boundaries = }")
    for header_row, end_row in pairwise(matrix_boundaries):
        start_row = header_row + 1
        num_rows = end_row - start_row
        output = np.loadtxt(
            example_output_dir / filename,
            skiprows=start_row,
            max_rows=num_rows,
        )
        reference_output = np.loadtxt(
            example_output_dir / f"{filename}_ref",
            skiprows=start_row,
            max_rows=num_rows,
        )
        assert output.shape == reference_output.shape
        assert output + 1 == pytest.approx(reference_output, rel=1e-6, abs=1e-6)
    assert False
