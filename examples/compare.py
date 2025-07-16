import sys
import numpy as np

def extract_blocks(filename):
    blocks = {}
    current_key = None
    current_data = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.endswith(':'):
                if current_key and current_data:
                    blocks[current_key] = np.array(current_data)
                current_key = line[:-1]
                current_data = []
            else:
                try:
                    row = [float(x) for x in line.split()]
                    current_data.append(row)
                except ValueError:
                    continue
        if current_key and current_data:
            blocks[current_key] = np.array(current_data)

    return blocks

def compare_blockwise(file1, file2, decimal_places=10):
    data1 = extract_blocks(file1)
    data2 = extract_blocks(file2)

    all_keys = set(data1.keys()).union(data2.keys())
    tolerance = 10**(-decimal_places)
    success = True

    for key in sorted(all_keys):
        if key not in data1 or key not in data2:
            print(f"Failure! Block '{key}' missing in one of the files.")
            success = False
            continue

        A, B = data1[key], data2[key]

        if A.shape != B.shape:
            print(f"Failure! Shape mismatch in block '{key}': {A.shape} vs {B.shape}")
            success = False
            continue

        diff = np.abs(A - B) > tolerance
        if np.any(diff):
            indices = np.argwhere(diff)
            print(f"Failure! Block '{key}' has {len(indices)} mismatches (>{decimal_places} decimals):")
            for i, j in indices[:5]:  # Show first 5 mismatches
                print(f"   [{i},{j}]: {A[i,j]:.12f} vs {B[i,j]:.12f}")
            if len(indices) > 5:
                print(f"   ... and {len(indices) - 5} more")
            success = False
        else:
            print(f"Success! Block '{key}' matches up to {decimal_places} decimal places.")
    
    return success

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    compare_blockwise(file1, file2, decimal_places=10)

