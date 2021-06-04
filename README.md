# Installation

1. Clone this repo.
2. Download itensor to arbitrary directory.
3. Build itensor against your openblas/mkl/etc.
4. Move the resulting "itensor" and "lib" folders to the following location(create new folders whenever needed):

itensor -> extern/itensor/include/itensor

lib -> extern/itensor/lib

5. Change the openblas/mkl/etc location in CMakeLists.txt, matching ITensor options.mk. For example,

-L/usr/lib64 -> -L/usr/local/opt/openblas/lib

-I/usr/include/openblas -> -I/usr/local/opt/openblas/include

For cmake versions ealier than 2.8, you can manually cmake_minimum_required, but you may need to explicitly include extra libraries. If GNU version is older version, you may need to Add -lstdc++fs to ITENSOR_LINK_FLAGS.

6. Make a folder called build. Change directory into it.
7. Run "cmake ..". This should generate a makefile.
8. Run make. The program by default stores data direcrtories under the parent directory of build. To change, modify the variable prefix_ in src/dmrg.h.

# Usage

The built binary accepts the following arguments:

* -s, --site: type of sites [std::string]
* -b, --bc: type of boundary condition [std::string]
* -n: number of sites [int]
* -d: max bond dimension [int]
* -c, --svd_cutoff_: svd_cutoff_ [float]
* -t, --tol: stable_tol_ for algorithm [float]
* --theta_: angle away from identity projector [float]
* --phi_: angle between rho and a*rho projectors [float]
* -u, --penalty: coefficient for penalty terms [float]
* --nstates: number of excited states to solve for [int]
* --analysis: if set to 1, Run the analysis using existing files. If set to 0, Run the dmrg. [int]
