# Installation

1. Clone this repo.
2. Download itensor to arbitrary directory.
3. Build itensor against your openblas/mkl/etc.
4. Move the resulting "itensor" and "lib" folders to the following location:

itensor -> extern/itensor/include/tensor

lib -> extern/itensor/lib

5. Change the openblas/mkl/etc location in the CMakeLists.txt
6. Run cmake

# Usage

The built binary accepts the following arguments:

* -s, --site: type of sites [std::string]
* -b, --bc: type of boundary condition [std::string]
* -n: number of sites [int]
* -d: max bond dimension [int]
* -c, --cutoff: cutoff [float]
* -t, --tol: tolerance for algorithm [float]
* --theta: parameter used in constructing Hamiltonian [float]
* -u, --penalty: coefficient for penalty terms [float]
* --nstates: number of excited states to solve for [int]
* --analysis: if set to 1, run the analysis using existing files. If set to 0, run the dmrg. [int]