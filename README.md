# Installation

1. Clone this repo.
2. Download itensor to arbitrary directory.
3. Build itensor against your openblas/mkl/etc.
4. Move/copy the resulting "itensor" and "lib" folders to the following locations (create new folders whenever needed):

itensor -> extern/itensor/include/itensor

lib -> extern/itensor/lib

5. Duplicate CMakeListsTemplate.txt as CMakeLists.txt. Modify the openblas/mkl/etc location in CMakeLists.txt, matching ITensor options.mk. For example, the following change is appropriate for the Harvard cluster:

-L/usr/lib64 -> -L/usr/local/opt/openblas/lib

-I/usr/include/openblas -> -I/usr/local/opt/openblas/include

Also included is a template CMakeListsMac.txt that works on Ying's Mac. Further modification may be necessary to correctly link libraries.

6. For cmake versions earlier than 2.8, you can manually modify cmake_minimum_required, but you may need to explicitly include extra libraries. If GNU version is older, you may need to add -lstdc++fs to ITENSOR_LINK_FLAGS.
7. Create a folder called build. Change directory into it.
8. Run "cmake ..". This should generate a makefile.
9. Duplicate src/path_template.h file as src/path.h. This file specifies the directory to store data, and is by default set to the Chain folder. To change this, modify the variable kPath in src/path.h.
10. Run make.

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
* -m, --mode: currently provides 5 modes [int];
  
  0: run dmrg simulation,
  
  1: measure energy, translation and rho eigenvalues,
  
  2: measure energy and translation eigenvalues, 
  
  3: measure energies (normalization dependent),

  4: measure energy ratios (normalization independent),

  5: run dmrg simulation while performing 1 (charge=0) or 2 (otherwise) after simulation of each state.

For example, the following command simulates a length 6 periodic Haagerup chain with QN charge conservation under pure rho projector in the neutral sector, and moreover measures energy, translation and rho eigenvalues after simulation of each state:

./Chain -s haagerupq -b p -l 6 -d 1600 -c 1e-8 -p 1e-2 -j 0,-1,0 -u 2 -q 0 -n 1 -m 5

The simulation and measurement results are stored in several files.  The .an file contains the diagonalized spectrum in the form of {energy, momentum, rho}, .ee file contains the entanglement entropy curve, .en file contains the energies without diagonalization, and .m file contains the energy, translation and rho matrices without diagonalization in Mathematica-compatible format.