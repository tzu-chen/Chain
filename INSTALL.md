# Installation & Usage

Full build procedure and CLI reference for the C++ solver. See the [README](README.md) for a project overview.

## Prerequisites

- A C++17-capable compiler (GCC ≥ 7, Clang ≥ 5)
- CMake ≥ 2.8
- BLAS/LAPACK: OpenBLAS, MKL, or equivalent
- [ITensor](https://itensor.org/) v3 (built against the same BLAS as above)
- For older GCC: you may need to add `-lstdc++fs` to `ITENSOR_LINK_FLAGS`

## Building

### 1. Clone

```bash
git clone https://github.com/tzu-chen/Chain.git
cd Chain
```

### 2. Build ITensor

Download ITensor to any directory and build it against your chosen BLAS. Follow ITensor's own documentation for this step — make a note of the exact `options.mk` you used, since `CMakeLists.txt` must match it.

### 3. Vendor ITensor into this repo

Copy the resulting `itensor` headers and `lib` directory into `extern/`, creating directories as needed:

```
extern/itensor/include/itensor/   <- ITensor headers
extern/itensor/lib/               <- compiled libraries
```

### 4. Configure CMake

Duplicate the template:

```bash
cp CMakeListsTemplate.txt CMakeLists.txt
```

Edit `CMakeLists.txt` to point at your BLAS installation. The paths must match the `options.mk` used when building ITensor. For example, on the Harvard cluster:

```
-L/usr/lib64             ->  -L/usr/local/opt/openblas/lib
-I/usr/include/openblas  ->  -I/usr/local/opt/openblas/include
```

A working macOS variant is provided as `CMakeListsMac.txt` — use it as a starting point if you're on Mac.

For CMake older than 2.8, lower `cmake_minimum_required` manually; you may also need to declare extra libraries explicitly.

### 5. Configure the data path

```bash
cp src/path_template.h src/path.h
```

Edit `kPath` in `src/path.h` to set where simulation output is written. The default is the `Chain/` folder itself.

### 6. Build

```bash
mkdir build && cd build
cmake ..
make
```

The binary is written to `build/Chain`.

## Julia version

No build step — just install the project's dependencies:

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

Tested with Julia 1.x.

---

## CLI reference

```
./Chain [options]
```

### Model and geometry

| Flag | Long form | Type | Description |
|---|---|---|---|
| `-s` | `--site` | string | Fusion category: `golden`, `haagerup`, `haagerupq` (Haagerup with QN conservation) |
| `-b` | `--bc` | string | Boundary condition: `p` (periodic), `o` (open), `d` (Dirichlet), `s` (sine-squared deformed), `sp` (SSD → periodic hot start) |
| `-l` |  | int | Chain length $L$ (number of sites) |
| `-q` | `--charge` | int | QN charge, only meaningful for `haagerupq` (e.g., `0` for the neutral sector) |

### DMRG parameters

| Flag | Long form | Type | Description |
|---|---|---|---|
| `-d` |  | int | Maximum bond dimension |
| `-c` | `--cutoff` | float | SVD truncation cutoff (e.g., `1e-8`) |
| `-p` | `--precision` | float | Energy-minimization precision, divided by $L$ internally (e.g., `1e-2`) |
| `-n` | `--nstates` | int | Number of states to compute, including the ground state (`1` = ground state only) |

### Hamiltonian

| Flag | Long form | Type | Description |
|---|---|---|---|
| `-j` | `--couplings` | string | Comma-separated projector coefficients. Positive = ferromagnetic, negative = antiferromagnetic. For Haagerup: `-1,0,0` = pure identity projector, `0,-1,0` = pure $\rho$ projector, `0,0,-1` = pure $a\rho$ projector |
| `-u` | `--penalty` | float | Coefficient of the orthogonality penalty term used when targeting excited states (e.g., `2`) |

### Mode

| Flag | Long form | Type | Description |
|---|---|---|---|
| `-m` | `--mode` | int | See modes below |

**Modes:**

| Value | Action |
|---|---|
| `0` | Run DMRG simulation |
| `1` | Measure energy, translation, and $\rho$ eigenvalues |
| `2` | Measure energy and translation eigenvalues |
| `3` | Measure raw energies (normalization-dependent) |
| `4` | Measure energy ratios (normalization-independent) |
| `5` | Run DMRG + perform mode 1 (if `charge=0`) or mode 2 (otherwise) after each state |
| `6` | Recover already-simulated states from a corrupt vector |
| `7` | Perform direct diagonalization of the Hamiltonian |

## Output files

Results are written under the `kPath` directory configured in `src/path.h`.

| Extension | Contents |
|---|---|
| `.an` | Diagonalized spectrum: `{energy, momentum, ρ}` tuples |
| `.ee` | Entanglement entropy curve across the chain |
| `.en` | Raw (undiagonalized) energies |
| `.m`  | Energy / translation / $\rho$ matrices, Mathematica-compatible |

## Example runs

**Periodic Haagerup chain with QN conservation, pure $\rho$ projector, neutral sector, ground state with full measurement:**

```bash
./Chain -s haagerupq -b p -l 6 -d 1600 -c 1e-8 -p 1e-2 -j 0,-1,0 -u 2 -q 0 -n 1 -m 5
```

**Golden (Fibonacci) chain, periodic, lowest 3 states, measurement only (no DMRG):**

```bash
./Chain -s golden -b p -l 12 -d 400 -c 1e-10 -p 1e-2 -j -1 -u 2 -n 3 -m 1
```

## Julia usage

DMRG driver:

```bash
julia --project=julia julia/run_chain.jl \
    -l 6 -d 100 -s 5 -b p -j 1.0 -u 2.0 -m haagerup -o results.jld2
```

Flags:

- `-m` — model: `golden` (default) or `haagerup`
- `-n` — number of excited states to target
- `-w` — orthogonality weight for excited-state targeting
- `-o` — output file (JLD2)

Programmatic use with a user-defined fusion category:

```julia
using Chain.Model
m = AnyonModel(2; fsymbols=myF, qdims=[1.0, phi])
```

Or the built-in Fibonacci helper:

```julia
using Chain
energy = run_dmrg(L=6, model=Model.fibonacci_model(), U=2.0)
```

The F-move operator is constructed directly from the supplied F-symbols, so no constants (like $\varphi$) are hard-coded into the driver.

## Troubleshooting

**`ld: library not found for -lopenblas`** — Check that the `-L` path in `CMakeLists.txt` matches your actual BLAS install location and that it matches ITensor's `options.mk`.

**`std::filesystem` link errors on older GCC** — Add `-lstdc++fs` to `ITENSOR_LINK_FLAGS`.

**CMake version too old** — Lower `cmake_minimum_required` manually, and expect to declare some linked libraries explicitly that newer CMake would have figured out.

**`kPath` not writable / simulations produce no output** — Verify `src/path.h` was created from the template and that `kPath` points somewhere you have write access.
