# Installation

1. Clone this repo.
2. Download itensor to arbitrary directory.
3. Build itensor against your openblas/mkl/etc.
4. Move/copy the resulting "itensor" and "lib" folders to the following locations (create new folders whenever needed):

itensor -> extern/itensor/include/itensor

lib -> extern/itensor/lib

5. Modify the openblas/mkl/etc location in CMakeLists.txt, matching ITensor options.mk. For example, the following change is appropriate for the Harvard cluster:

-L/usr/lib64 -> -L/usr/local/opt/openblas/lib

-I/usr/include/openblas -> -I/usr/local/opt/openblas/include

Also included is a file CMakeLists_Mac.txt that works on Ying's Mac. To use it, rename it to CMakeLists.txt. Further modification may be needed to correctly link libraries.

6. For cmake versions earlier than 2.8, you can manually modify cmake_minimum_required, but you may need to explicitly include extra libraries. If GNU version is older, you may need to add -lstdc++fs to ITENSOR_LINK_FLAGS.
7. Create a folder called build. Change directory into it.
8. Simulate "cmake ..". This should generate a makefile.
9. The program by default creates and stores data to directories under the Chain folder. To change this, modify the variable kPath in src/path.h.
10. Simulate make.

# Usage

The built binary accepts the following arguments:

* -s, --site: type of sites [std::string]; golden, haagerup, haagerupq (haagerup with QN conservation).
* -b, --bc: type of boundary condition [std::string]; p(eriodic), o(pen), s(ine-squared deformed), sp (sine-squared deformed to hot start periodic).
* -l: L = length of chain = number of sites [int].
* -d: max bond dimension [int].
* -c, --cutoff: cutoff for singular value decomposition [float]; e.g., 1e-8.
* -p, --precision: precision of energy minimization is this value divided by length [float]; e.g., 1e-2.
* -j, --couplings: comma-separated string containing the coupling coefficients of projectors, positive for ferromagnetic and negative for antiferromagnetic [std::string]; e.g., for Haagerup, -1,0,0 is pure identity projector, 0,-1,0 is pure rho projector, 0,0,-1 is pure a*rho projector.
* -u, --penalty: coefficient for the penalty term [float]; e.g., 2.
* -q, --charge: QN charge (only valid for haagerupq) [int]; e.g., 0.
* -n, --nstates: number of excited states to simulate/analyze [int]; e.g., 1 for just the ground state.
* -m, --mode: currently provides 4 modes [int];
  
  0: run the dmrg simulation,
  
  1: measure energy, translation, and rho eigenvalues (measuring rho takes roughly (dim H)^(L-4) seconds), 
  
  2: measure energy and translation eigenvalues, 
  
  3: compute energy ratios.
