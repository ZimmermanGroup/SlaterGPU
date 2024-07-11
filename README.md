# SlaterGPU

Library for numerically computing Slater-type orbital integrals.
For running on GPU, OpenACC is required. Code has been tested
and compiled with Nvidia HPC SDK 20.7, 21.7 and 21.9. A wrapper
library for Libcint (https://github.com/sunqm/libcint) is also
provided for Gaussian integrals. This library uses the resolution
of the identity (RI) approximation. Thus, auxiliary and main
basis sets must be specified.

For STOs, functions up to 6h (5g for derivatives) are available 
for the auxiliary basis set. With up to 6h auxiliary functions,
it is not recommended to go beyond 4f in the main basis. The 
library requires user supplied basis sets. See the examples folder 
for formatting inputs. Atoms up to Zn are currently supported.

Note: This repository is an experimental research project, and functionality related to the prolate spheroidal coordinate system is still in progress and particularly unstable.

### Example of getting dependencies on Zimmerman group cluster

Source the script that automatically loads modules on the cluster:

```
source go
```

For other users, this bash script may be a helpful example to see how to install / load necessary dependencies.

### For compiling run :
```
mkdir build
cd build
cmake ..
make
```

The executable will be generated at `build/examples/sgpu.exe`

### Testing:
Go to an example in the build directory. From the build directory, this would be:
```
cd examples/geom_ps/
```

Then run the SlaterGPU executable which will read the input files from the current working directory:
```
../sgpu.exe
```

### Notes:
If running on multiple GPUs, it's advised to have the number of
OpenMP threads equal to the number of GPUs. i.e. set the following
environment variable
```
export OMP_NUM_THREADS=<ngpu>
```

An example for computing the integrals is provided in the 
`examples` folder.

Please see LICENSE file for licensing information.
