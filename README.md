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

For compiling run 
```
mkdir build
cd build
cmake ..
make
```
By default, the integration grid is evaluated in double precision.
For additional performance, one can specify mixed precision at
compile time by setting `-DEVL64=0` in `CMake_CXX_FLAGS` i.e.
by passing `-DCMake_CXX_FLAGS="-DEVL64=0"` when configuring.

By default the GTO wrapper library is not built. To build this
wrapper include `-DDO_GTO=True` as a CMake flag, and set an 
environment variable `LIBCINT_PATH` to your Libcint install.
When using a GTO basis, a main basis file named `basis` and an
auxiliary basis file named `aux` are required. These files
should be in Gaussian format and include the specification for
all atoms of interest.

An example for computing the integrals is provided in the 
`examples` folder.

Please see LICENSE file for licensing information.
