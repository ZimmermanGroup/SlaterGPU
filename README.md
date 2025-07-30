# SlaterGPU

Library for numerically computing Slater-type orbital integrals.
For running on GPU, OpenACC is required. A compatible compiler should be used.
Testing has been done with NVHPC SDK version 24.9. A wrapper
library for Libcint (https://github.com/sunqm/libcint) is also
provided for Gaussian integrals. It should be noted that version 5.3.0 of libcint is required at the moment.

This library uses the resolution
of the identity (RI) approximation. Thus, auxiliary and main
basis sets must be specified.

For STOs, functions up to 6h (5g for derivatives) are available 
for the auxiliary basis set. With up to 6h auxiliary functions,
it is not recommended to go beyond 4f in the main basis. The 
library requires user supplied basis sets. See the examples folder 
for formatting inputs. Atoms up to Zn are currently supported.

Note: This repository is an experimental research project, and functionality related to the prolate spheroidal coordinate system is still in progress and particularly unstable.

### Example of getting dependencies on Zimmerman group cluster

Source the script that automatically loads modules and adds dependencies to PATH:

```
source env.set.local0
```

For other users, this bash script may be a helpful example to see how to install / load necessary dependencies.

### For compiling run :
```
cmake -Bbuild
cmake --build build -j 16
```

The executable will be generated at `build/examples/sgpu.exe`

### Testing:
Go into one of the example directories inside the build directory.
```
cd build/examples/lih_VK1
```

Then execute SlaterGPU within this directory to run the test.
```
../sgpu.exe
```

After running, if successful the directory should contain the output files "A", "Ciap", "SENT", and "pVp".

### Notes:
If running on multiple GPUs, it's advised to have the number of
OpenMP threads equal to the number of GPUs. i.e. set the following
environment variable
```
export OMP_NUM_THREADS=<ngpu>
```

There are example calculations in `SlaterGPU/examples/` with integral files denoted generally as `INT_ref`. It should be noted that normally, the three center coefficients, `Ciap`, are generated, but were omitted due to Git's file size limit.

Please see LICENSE file for licensing information.

## Experimental pixi build:
- Clone this repository
- Install [pixi](https://pixi.sh/latest/installation/) if you don't already have it (`curl -fsSL https://pixi.sh/install.sh | sh`)
- Ensure NVIDIA HPC SDK is available on the path
  - On athena: `module load nvidia-sdk/24.11-multi`
  - On perlmutter:
    ```
    module use /global/cfs/cdirs/m4957/joshkamm/hpc_sdk/modulefiles
    module unload cudatoolkit
    module load nvhpc/24.11
    ```
- Install dependencies and build ZEST with `pixi -v install` (`-v` to see output from pixi and cmake build)
- If the install has succeeded, activate the pixi environment with `pixi shell` inside the folder where this repository is cloned
- Run the zest executable using `mpirun -n <NUMBER_OF_THREADS> zest` in a directory with the appropriate input files (including the integral files SENT, A, Ciap).
- Please note that MPI parallel is not available except for NHBCI and iFCI.
Pixi should automatically install whatever dependencies are necessary to run a pixi task (except NVHPC, which cannot yet be installed automatically). Run tasks with `pixi run <TASK>` where \<TASK\> is one of the following:
```
 - clean           Remove example output files
 - test            Test SlaterGPU executable on a small molecular system
 - test-in-gh-action Test SlaterGPU executable on a small molecular system in a github actions environment
```
Above list generated with `pixi task list`

Configuration of tasks and dependencies can be found in the `pixi.toml` file.

