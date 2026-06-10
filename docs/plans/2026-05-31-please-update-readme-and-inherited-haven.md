# Plan: Update README and Prepare User Documentation

## Context

OUxSBLI's existing documentation has drifted from the codebase. The build system migrated from `Makefile` to CMake, and scheme/physics configuration moved from `mod_globals.f90` to `config.fypp` — but neither README.md nor tutorials/installation.md reflects these changes. The user wants a coherent set of user-facing Markdown docs under `docs/` targeting external researchers approaching the project for the first time.

---

## Files to Create or Modify

| File | Action |
|------|--------|
| `README.md` | Update in place |
| `docs/quickstart.md` | Create new |
| `docs/installation.md` | Create new (supersedes `tutorials/installation.md`) |
| `docs/configuration.md` | Create new |
| `docs/theory.md` | Create new |
| `docs/api.md` | Create new |

`tutorials/installation.md` is left as-is (not deleted) since it may be linked from elsewhere; the README will link to the new `docs/` versions.

---

## 1. README.md — Targeted Fixes

The README is public-facing and needs to be accurate first.

### Changes

**Dependencies section** — add `fypp`, `cmake`, note that MPI comes bundled with HPC SDK.

**Usage — CUDA Fortran only** — replace the 5-step block:

Old step 2 ("Edit mod_globals.f90 — Choose equation type `id_visc`…") is wrong; those are now `config.fypp` settings.  
Old step 4 (`make`) is wrong; build is now CMake.

Replace with:
```
1. cd 3D_solver/NSTGV
2. Edit config.fypp — choose VISC, SCHEME, ORDER, BC_X/Y/Z
3. Edit mod_globals.f90 — grid size, physical parameters, thread-block sizes
4. cmake -B build && cmake --build build -j
5. cd build && bash ../calc.sh          # mpirun -n 2 ./a.out
```

**Add curvilinear solver** — one-sentence mention with its two cases (NACA, CORN).

**Add docs links** — table or bullet list pointing to `docs/quickstart.md`, `docs/installation.md`, `docs/configuration.md`, `docs/theory.md`, `docs/api.md`.

**No structural changes** — keep badges, images, Discretization, Validations, Publication, License sections as-is.

---

## 2. docs/quickstart.md — New

Target: someone who already has HPC SDK installed and wants to run their first simulation in under 10 minutes.

Sections:
1. Prerequisites (HPC SDK, CMake 3.18+, 2 MPI ranks / 1 GPU)
2. Clone the repo
3. Pick a case — recommend `3D_solver/NSTGV` (NS Taylor-Green vortex, simplest full 3D case)
4. Build (CMake)
5. Run (`mpirun -n 2 ./a.out`)
6. Visualise output in ParaView (open `Q*.vtr` files as a series)
7. Next steps — links to configuration.md to change the scheme

---

## 3. docs/installation.md — New (replaces tutorials/installation.md)

Full installation guide for a fresh machine. CMake-based build throughout.

Sections:
1. System requirements table (Ubuntu 22.04, NVIDIA Ampere+, CUDA 12.x, HPC SDK 24–25, CMake 3.18+, Python 3.10+)
2. Install NVIDIA HPC SDK
3. Install fypp (`pip install fypp` — required by CMake)
4. Clone repository
5. (Optional) Install Python package: `pip install -e ".[dev]"`
6. Verify: `pytest ouxsbli/tests/ -v`
7. Build a test case (CMake workflow, using `3D_solver/NSTGV`)
8. Run (`cd build && mpirun -n 2 ./a.out`)
9. Available cases — three complete tables (2D, 3D, curvilinear) matching CLAUDE.md

Fix errors in `tutorials/installation.md`:
- `make clean && make` → CMake
- Missing 2D cases: OS, SBLI
- Missing 3D cases: IVST, STZ
- Scheme/physics table references `mod_globals.f90` — update to note those are generated from `config.fypp`

---

## 4. docs/configuration.md — New

Reference guide for the compile-time configuration system.

Sections:
1. **How it works** — fypp template expansion, config.fypp → .f90.fypp → build/*.f90, no runtime branching
2. **config.fypp reference table** — all variables (VISC, SCHEME, ORDER, VISC_ORDER, TVD, SLAU_VARIANT, RESCALE, COMMZ, RK, RESTART, BC_X/Y/Z, ORDER_IO) with allowed values and effect — taken from CLAUDE.md
3. **mod_globals.f90 reference** — grid size, domain lengths, physical params, GPU thread-block sizes
4. **mod_constant.f90 dispatch encoding** — id_visc / id_scheme / id_accuracy kind-dispatch table (from CLAUDE.md)
5. **Rebuild workflow** — edit config.fypp → `cmake --build build -j` (CMake detects change)
6. **Worked example** — switch NSTGV from KEEP to Hybrid scheme, step by step

---

## 5. docs/theory.md — New

Numerical methods reference for researchers.

Sections:
1. **Governing equations** — compressible Euler / NS / LES, conservative form, Q array layout
2. **Convective schemes**
   - KEEP — kinetic energy and entropy preserving, best for smooth vortical flows
   - SLAU / HRSLAU2 — low-dissipation AUSM, compressible turbulence with shocks
   - Roe — approximate Riemann solver, hypersonic/strong discontinuities
   - Hybrid — KEEP ↔ SLAU switched by Ducros sensor
3. **Viscous schemes** — Gaitonde-Visbal 2nd order; ME4-Base 4th order compact
4. **Time integration** — TVD-RK3 (3-stage), classical RK4
5. **MPI decomposition** — x-direction default (1D), z-direction COMMZ with overlap
6. **Key references** — Lusher & Sandham (2021), Hatayama et al. (2026), and scheme sources

---

## 6. docs/api.md — New

Python API reference for `ouxsbli`.

Sections:
1. **Overview** — what the Python API does (parametric case setup, build, run)
2. **Install** — `pip install -e ".[dev]"`
3. **Case class** — constructor parameters table (source, workdir, visc, scheme, accuracy, nx/ny/nz, etc.)
4. **Methods** — `case.build()`, `case.run(nranks=2)`
5. **Full example** — clone → Case(...) → build() → run() → read VTK with vtk_reader
6. **Test suite** — how to run `pytest ouxsbli/tests/`, what each test file covers

---

## Verification

After implementation, check:
1. All CMake build commands are syntactically correct (`cmake -B build && cmake --build build -j`)
2. All case names match directory names on disk (DHIT, ETGV, IVST, KHI, NSTGV, SBLI, STZ, TBL / BL, DSL, EVC, OS, SBLI, ST / NACA, CORN)
3. config.fypp variable names and values match `3D_solver/NSTGV/config.fypp`
4. Python API examples match `ouxsbli/case.py` and `ouxsbli/patcher.py`
5. README.md links to all new docs/ files (relative paths)
6. No remaining references to `make` or `make clean` in the new docs (except curvilinear Makefile where it is still correct)
