from itertools import pairwise
import os
from pathlib import Path
import numpy as np
import pytest

@pytest.mark.parametrize("filename", ["A", "Ciap", "pVp", "SENT"])
def test_output(filename):
    example_output_dir = Path("lih_VK1")

    with open(example_output_dir / f"{filename}_ref") as reference_output:
        reference_output_lines = reference_output.readlines()

    matrix_boundaries = [
        line_num
        for line_num, line in enumerate(reference_output_lines)
        # Get header lines for matrices and blank last line.
        if line.strip().endswith(":")
    ] + [len(reference_output_lines)]

    # Loop over each matrix in the file and compare to the reference.
    for header_row, end_row in pairwise(matrix_boundaries):
        start_row = header_row + 1
        num_rows = end_row - start_row
        matrix = np.loadtxt(
            example_output_dir / filename,
            skiprows=start_row,
            max_rows=num_rows,
        )
        reference_matrix = np.loadtxt(
            example_output_dir / f"{filename}_ref",
            skiprows=start_row,
            max_rows=num_rows,
        )
        assert matrix == pytest.approx(reference_matrix, rel=1e-6, abs=1e-6)
