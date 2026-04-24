# Chain

A high-performance DMRG (density matrix renormalization group) solver for 1D **anyon chains** — quantum many-body systems whose Hilbert spaces are *not* tensor-factorized over sites, and therefore cannot be handled by standard spin-chain codes out of the box.

Written in C++ on top of [ITensor](https://itensor.org/), with a later Julia port for prototyping new fusion categories.

Research output: [*Numerical Evidence for a Haagerup Conformal Field Theory*](https://arxiv.org/abs/2110.03008), Phys. Rev. Lett. **128**, 231603 (2022).

---

## What it does

Given a fusion category (the algebraic data defining the anyon model) and a chain geometry, the solver:

- Constructs the Hamiltonian as a matrix product operator (MPO) directly from the category's F-symbols — no hard-coded constants.
- Runs DMRG to find the ground state and low-lying excited states, targeting excited states via orthogonality penalties.
- Measures translation and topological-symmetry eigenvalues alongside the energy spectrum, so that states can be organized into momentum and symmetry sectors.
- Supports multiple boundary conditions: periodic, open, Dirichlet, and sine-squared deformed (SSD), including the SSD-to-periodic "hot start" protocol.
- Handles quantum-number (QN) conservation where the category admits a grading, which dramatically reduces bond dimension at fixed accuracy.

Built-in models: **golden** (Fibonacci), **Haagerup**, **Haagerup with QN conservation**. The Julia side additionally supports user-defined fusion categories via a generic `AnyonModel` interface.

## Why it's non-trivial

Anyonic Hilbert spaces are labeled by fusion trees rather than local spin configurations. This means:

- Local operators don't decompose as tensor products; every term in the Hamiltonian involves an F-move that depends on the global fusion-tree state.
- The MPO construction has to be generated from categorical data at runtime.
- Symmetry sectors are labeled by non-invertible operators (e.g., the $\rho$ line of the Haagerup category), which the code diagonalizes jointly with energy and translation.

The solver reaches bond dimensions of $\sim 1600$ on chains up to $L \sim 16$ with SVD cutoffs of $10^{-8}$, which was sufficient to extract the central charge of the emergent CFT in the published work.

## Tech stack

| Layer | Tool |
|---|---|
| Core solver | **C++17**, [ITensor](https://itensor.org/) v3 |
| Linear algebra backend | OpenBLAS / MKL (configurable at build time) |
| Build system | CMake ≥ 2.8 |
| Data format | ITensor binary (C++), JLD2 (Julia) |
| Prototyping / generic categories | **Julia 1.x** |
| Output analysis | Mathematica-compatible matrix dumps (`.m`) |

## Repository layout

```
src/          C++ DMRG driver, MPO construction, measurement modes
julia/        Julia port with generic AnyonModel interface
extern/       Vendored ITensor headers and libraries (see INSTALL.md)
CMakeLists*   Build configuration (template + Mac variant)
```

## Install

See [INSTALL.md](INSTALL.md) for the full build procedure. In short: vendor ITensor under `extern/itensor/`, point `CMakeLists.txt` at your BLAS installation, then:

```bash
mkdir build && cd build
cmake ..
make
```

The Julia version has no build step:

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

## Quick start (C++)

Simulate a length-6 periodic Haagerup chain with QN conservation under the pure-$\rho$ projector, in the neutral sector, and measure the energy / translation / $\rho$ spectrum after each state:

```bash
./Chain -s haagerupq -b p -l 6 -d 1600 -c 1e-8 -p 1e-2 -j 0,-1,0 -u 2 -q 0 -n 1 -m 5
```

Output files:

| Extension | Contents |
|---|---|
| `.an` | Diagonalized spectrum: `{energy, momentum, ρ}` |
| `.ee` | Entanglement entropy curve |
| `.en` | Raw (undiagonalized) energies |
| `.m`  | Energy / translation / ρ matrices in Mathematica format |

Full CLI reference is in [INSTALL.md](INSTALL.md#cli-reference).

## Quick start (Julia)

```bash
julia --project=julia julia/run_chain.jl -l 6 -d 100 -s 5 -b p -j 1.0 -u 2.0 -m haagerup -o results.jld2
```

Or, using the generic interface to run a Fibonacci chain:

```julia
using Chain
energy = run_dmrg(L=6, model=Model.fibonacci_model(), U=2.0)
```

Define your own fusion category by supplying F-symbols and quantum dimensions:

```julia
using Chain.Model
m = AnyonModel(2; fsymbols=myF, qdims=[1.0, phi])
```

## Citing

If this code is useful for your work, please cite:

```bibtex
@article{HuangLinOhmoriTachikawaTezuka2022,
  title   = {Numerical Evidence for a Haagerup Conformal Field Theory},
  author  = {Huang, Tzu-Chen and Lin, Ying-Hsuan and Ohmori, Kantaro
             and Tachikawa, Yuji and Tezuka, Masaki},
  journal = {Phys. Rev. Lett.},
  volume  = {128},
  number  = {23},
  pages   = {231603},
  year    = {2022},
  doi     = {10.1103/PhysRevLett.128.231603},
}
```
