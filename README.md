# Installation

1. Clone this repo.
2. Download itensor to arbitrary directory.
3. Build itensor against your openblas/mkl/etc.
4. Move the resulting "itensor" and "lib" folders to the following location (create new folders whenever needed):

itensor -> extern/itensor/include/itensor

lib -> extern/itensor/lib

5. Change the openblas/mkl/etc location in CMakeLists.txt, matching ITensor options.mk. For example, the following change is appropriate for the Harvard cluster:

-L/usr/lib64 -> -L/usr/local/opt/openblas/lib

-I/usr/include/openblas -> -I/usr/local/opt/openblas/include

Also included is a file CMakeLists_Mac.txt that works on Ying's Mac. To use it, rename it to CMakeLists.txt. Further modification may be necessary to correctly link libraries.

6. For cmake versions earlier than 2.8, you can manually modify cmake_minimum_required, but you may need to explicitly include extra libraries. If GNU version is older version, you may need to add -lstdc++fs to ITENSOR_LINK_FLAGS.
7. Make a folder called build. Change directory into it.
8. Run "cmake ..". This should generate a makefile.
9. The program by default stores data directories under the parent directory of build. To change, modify the variable kPath in src/path.h.
10. Run make.

# Usage

The built binary accepts the following arguments:

* -s, --site: type of sites [std::string]
* -b, --bc: type of boundary condition [std::string]
* -n: number of sites [int]
* -d: max bond dimension [int]
* -c, --svd_cutoff_: singular value decomposition cutoff [float]
* -t, --tol: stable_tol_ for algorithm [float]
* --couplings: comma-separated string containing the coupling coefficients of projectors [std::string], e.g., for Haagerup, -1,0,0 is pure identity projector, 0,-1,0 is pure rho projector, 0,0,-1 is pure a*rho projector.
* -u, --penalty: coefficient for the penalty term [float]
* --nstates: number of excited states to simulate/analyze [int]
* --analysis: If set to 0, run the dmrg simulation. If set to 1, run the analysis using existing files. If set to 2, run the quick analysis of energy ratios using existing files. [int]
