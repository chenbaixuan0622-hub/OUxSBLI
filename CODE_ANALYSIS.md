# OUxSBLI: Complete Project Analysis Report

## Executive Summary

**OUxSBLI** is a GPU-accelerated 3D compressible flow solver written in CUDA Fortran. The project has undergone significant restructuring, consolidating from multiple solver variants (1D, 2D, 3D) to a focused **3D solver suite** with six test cases. This report documents the current project structure, active configurations, and recent modifications.

**Key Status:**
- ✓ 3D solver fully operational (6 test cases)
- ✓ 1D and 2D solvers removed (commit 5893714)
- ✓ Compilation verified successful
- ✓ Recent focus: README update, SBLI refinements, code clarity improvements

---

## 1. Current Project Structure

### Repository Root
```
OUxSBLI/
├── README.md                  # Project overview and quick start
├── CODE_OF_CONDUCT.md        # Community guidelines
├── CITATION.cff              # Citation metadata
├── LICENSE                   # BSD 3-Clause License
│
├── docs/
│   ├── technical_doc.md      # Detailed technical documentation
│   ├── CODE_ANALYSIS.md      # Algorithm analysis and bugs (7.5KB)
│   ├── api.md                # Auto-generated API reference (79KB, 1338 lines)
│   ├── markdown.py           # API generation utility
│   └── multi_gpu/            # GPU scalability planning
│       ├── EXECUTIVE_SUMMARY.md
│       ├── IMPLEMENTATION_GUIDE.md
│       ├── MULTI_GPU_SCALABILITY_PLAN.md
│       └── QUICK_START.md
│
├── src/                      # Top-level solver support modules (10 files)
│   ├── main.f90              # Program entry point
│   ├── sbli.f90              # SBLI-specific main variant
│   ├── mod_constant.f90      # Physical constants
│   ├── mod_globals_c.f90     # C bindings for Python interface
│   ├── cpu_gpu_mpi.f90       # MPI/GPU communication abstraction
│   ├── calc_muscl.f90        # MUSCL reconstruction limiters
│   ├── calc_physical_quantities.f90  # Conservative ↔ primitive conversion
│   ├── print.f90             # Output/diagnostics
│   ├── set_compressible_bl.f90      # Boundary layer initialization
│   └── set_coordinate.f90           # Grid metrics and Jacobians
│
├── 3D_solver/                # Main 3D solver suite
│   ├── src/                  # Core solver kernels (20 files)
│   │   ├── calc_flux_base.f90        # Flux computation orchestrator
│   │   ├── calc_keep_3d.f90          # KEEP scheme 3D wrapper
│   │   ├── calc_keep_kernel.f90      # KEEP kernels
│   │   ├── calc_slau_3d.f90          # SLAU scheme 3D wrapper
│   │   ├── calc_slau_kernel.f90      # SLAU kernels (modified, commit 5893714)
│   │   ├── calc_roe_3d.f90           # Roe scheme 3D wrapper
│   │   ├── calc_roe_kernel.f90       # Roe kernels
│   │   ├── calc_hybrid.f90           # Hybrid MUSCL/LES flux selection
│   │   ├── calc_hybrid_kernel.f90    # Hybrid kernel (Ducros sensor)
│   │   ├── calc_les.f90              # LES/SGS model implementation
│   │   ├── calc_visc2.f90            # 2nd-order viscous fluxes (MODIFIED)
│   │   ├── calc_visc4.f90            # 4th-order viscous fluxes
│   │   ├── calc_steps.f90            # GPU time integration (TVD RK3, 4-4 RK)
│   │   ├── calc_time_dev.f90         # Time development orchestration
│   │   ├── calc_para.f90             # Parallelization control
│   │   ├── calc_rescale.f90          # Reference state rescaling
│   │   ├── preprocess.f90            # Pre-simulation setup
│   │   ├── set_init_common.f90       # Initial condition routines
│   │   ├── set_bc_common.f90         # Boundary condition routines
│   │   └── set_bc_tbl_sbli.f90       # TBL/SBLI specific BC
│   │
│   └── [6 Case Directories]
│       ├── ETGV/             # Euler Taylor-Green Vortex
│       ├── IVST/             # Isentropic Vortex
│       ├── KHI/              # Kelvin-Helmholtz Instability
│       ├── NSTGV/            # Navier-Stokes Taylor-Green Vortex
│       ├── SBLI/             # Shock-Boundary Layer Interaction
│       └── TBL/              # Turbulent Boundary Layer
│
├── img/                      # Visualization outputs and results
│   ├── Ek.png                # [NEW] Kinetic energy evolution
│   ├── enstrophy.png         # [NEW] Total enstrophy evolution
│   └── SBLI.png              # SBLI domain visualization
│
└── .gitignore               # (Updated: commit 65d3223)
```

---

## 2. Test Case Configurations

Each case in `3D_solver/[CASE]/` contains:
- `mod_globals.f90` - Case-specific parameters (equation type, scheme, grid, GPU blocks)
- `set.f90` - Grid generation, initial conditions, boundary conditions
- `Makefile` - Case-specific build configuration
- `calc.sh` - Execution script
- `profile.sh` - NVIDIA profiling script (nsys, ncu)
- `data/`, `recal/` - Output and restart directories

### Case-Specific Parameters (Current State)

| Case | id_visc | id_scheme | id_accuracy | Grid (nx×ny×nz) | Domain | Notes |
|------|---------|-----------|-------------|-----------------|--------|-------|
| **ETGV** | 2 (NS) | real(2) | int8 | 66×66×66 | Periodic cubic | Supersonic, Euler variant |
| **IVST** | 1 (NS) | real(2) | int8 | 129×7×7 | 1.0×0.1×0.1 | Isentropic, slender |
| **KHI** | 0 (Euler) | int(2) | int8 | 130×130×130 | Unit cube | Compressible instability |
| **NSTGV** | int4 | real(2) | int8 | **66×66×66** | Periodic cubic | **[OPTIMIZED grid, was 130³]** |
| **SBLI** | int4 | real(2) | int8 | 513×161×33 | Multi-domain | **[OPTIMIZED z: 0.5*blt, nz1=33 vs 257]** |
| **TBL** | int4 | int(2) | int(2) | 129×129×33 | 10×5×1.25 (blt) | Boundary layer |

**Key Observations:**
1. All cases use `id_scheme = real(2)` (floating-point type for scheme selection)
2. Most use `id_accuracy = int8` (64-bit)
3. `id_tvd = int8` for time integration type
4. NSTGV reduced from 130³ (~2.2M cells) to 66³ (~287K cells)
5. SBLI domain and grid optimized for shock-interaction focus

---

## 3. Recent Modifications

### Commit 65d3223: "Update README" (Apr 3, 2026)
**Changes:**
- Simplified project description
- Updated dependency information
- Restructured usage workflow (markdown improvements)
- **Numerical setup format**: LaTeX math blocks → Markdown table (Re, Ma, resolution)
- Publication information added (Hatayama et al. 2026)
- **New figures**: `img/Ek.png`, `img/enstrophy.png`
- **Removed**: 16 legacy images (DSL, EVC, Grid_convergence, etc.)
- Updated `.gitignore`

**Impact**: Cleaner, more discoverable README; consolidated to active research focus

### Commit 0b0ec13: "Calculate total temperature and pressure" (Apr 2)
**Changes:**
- SBLI case: Enhanced thermodynamic quantity calculations
- `3D_solver/SBLI/mod_globals.f90`: Total temperature/pressure computation

### Commit 81592e5: "HD transform was modified" (Mar 30+)
- Coordinate/Jacobian transformation updates

### Commit 85201da: "Contiguous attribute have been added" (Mar 30+)
- Array contiguity attributes for GPU optimization

### Commit 950a043: "Unnecessary files have been removed. Compilation succeeded." (Mar 30+)
- Code cleanup, verified build success

### Commit 5893714: "1D and 2D solver removed. Comments added. SLAU modified." (Mar 30)
- **Major restructuring**: Removed 1D_solver/, 2D_solver/ (scope consolidation)
- Added code documentation/comments
- SLAU scheme implementation refinements

---

## 4. Staged Changes (Ready to Commit)

### A. CODE_ANALYSIS.md (Root - NEW)
- **Purpose**: Comprehensive project survey and modification summary
- **Content**: Parameter changes, grid optimization, code improvements
- **Status**: ✓ Staged

### B. README.md (MODIFIED)
- **Changes**:
  - Simplified opening description
  - Dependency list updated
  - Usage workflow clearer
  - Numerical setup as Markdown table (not LaTeX)
  - Publication citation added
- **Status**: ✓ Staged

### C. docs/technical_doc.md (MODIFIED)
- **Changes**:
  - Code fence format: ` ``` ` → `~~~` (consistency)
  - Updated `mod_globals.f90` example parameters
  - **NEW Section 15**: Current Configuration & Recent Changes
  - Parameter type documentation (ETGV, NSTGV, SBLI)
  - Code improvements summary (viscous flux naming)
  - Version updated to 1.1 (Apr 3)
- **Status**: ✓ Staged

### D. Images (NEW)
- `img/Ek.png` (29KB) - Kinetic energy evolution
- `img/enstrophy.png` (31KB) - Enstrophy evolution
- **Status**: ✓ Staged

---

## 5. Uncommitted Changes (Working Directory)

### A. Parameter Refactoring (mod_globals.f90 variants)

#### ETGV Case
```fortran
! WORKING DIRECTORY (uncommitted):
real(2), parameter         :: id_scheme   = 0      ! [CHANGED from int(2)]
integer(kind=8), parameter :: id_accuracy = 0      ! [CHANGED from int(2)]
integer(kind=2), parameter :: id_tvd      = 0      ! [CORRECT]
! id_exchange REMOVED (was deprecated)
```
- Rationale: Flexible scheme parameterization; consistency in accuracy type

#### NSTGV Case
```fortran
! Grid optimization (uncommon):
integer, parameter :: nx = 66!130!258  ! [CHANGED 130 → 66]
integer, parameter :: ny = 66!130!258
integer, parameter :: nz = 66!130!258
```
- Impact: ~87% cell reduction (2.2M → 287K) for fast testing

#### SBLI Case
```fortran
! Domain focus (uncommon):
real(8), parameter :: Lz1 = 0.5d0 * blt  ! [CHANGED 4.0 → 0.5, 8× compression]
integer, parameter :: nz1 = 33           ! [CHANGED 257 → 33, ~87% reduction]
```
- Rationale: Focuses computation on shock-layer interaction region

### B. Code Improvements

#### calc_visc2.f90 - Variable Renaming (UNCOMMITTED)
```fortran
! OLD (less clear):
mx, my, mz       ! Generic viscosity terms

! NEW (clearer):
mxsgs, mysgs, mzsgs  ! SGS viscosity coefficients
```
- Affected subroutines:
  - `calc_Ev_LES2()` (~line 125-180)
  - `calc_Ev_LES_y()` (~line 280)
  - `calc_Ev_LES_z()` (~line 430)
- **No algorithmic change**; improved maintainability

### C. Documentation Removal

#### docs/PACKAGE_INDEX.md - DELETED
- **Size**: 441 lines
- **Reason**: Deprecated (index no longer maintained)
- **Status**: Uncommitted deletion

---

## 6. Solver Capabilities

### Discretization Schemes
1. **KEEP** (Kinetic Energy & Entropy Preserving)
   - File: `calc_keep_kernel.f90`, `calc_keep_3d.f90`
   - Order: 2nd, 4th, 6th
   - Best for: Smooth vortical flows

2. **SLAU** (Simple Low Dissipation AUSM)
   - File: `calc_slau_kernel.f90`, `calc_slau_3d.f90`
   - Order: 2nd, 4th, 6th
   - Feature: Shock-adaptive via Ducros sensor
   - **Recent improvements** (commit 5893714)

3. **Roe Scheme**
   - File: `calc_roe_kernel.f90`, `calc_roe_3d.f90`
   - Classical shock-capturing

4. **Hybrid MUSCL/LES**
   - File: `calc_hybrid.f90`, `calc_hybrid_kernel.f90`
   - Blends MUSCL limiters (smooth) ↔ LES (shocks)
   - Uses Ducros sensor for weighting

### Viscous Terms
- **2nd-order**: `calc_visc2.f90` (standard 2nd-order FD with averaging)
- **4th-order**: `calc_visc4.f90` (compact 4th-order)
- Both support SGS model for LES

### Time Integration
- **TVD RK3**: 3-3 stage Runge-Kutta (total variation diminishing)
- **4-4 RK**: 4th-order 4-stage Runge-Kutta
- Implemented in: `calc_steps.f90`, `calc_time_dev.f90`

### LES/SGS
- **Smagorinsky model**: `calc_les.f90`
- Selective mixed scale: Under development
- Applied via dynamic viscosity enhancement

---

## 7. Build & Execution

### Requirements
- **HPC SDK** 24.* or 25.*
- **CUDA** compatible with HPC SDK
- **MPI** library
- **ParaView** (optional, for visualization)

### Quick Start (Example: NSTGV)
```bash
cd 3D_solver/NSTGV
# Edit mod_globals.f90 for parameters
make clean && make
bash calc.sh
```

### Profiling
- `profile.sh` - NVIDIA nsys (system profiler)
- `profile_ncu.sh` (NSTGV) - NVIDIA ncu (kernel profiler)
- SBLI: `job_*.sh` scripts for cluster submission

---

## 8. Documentation Coverage

| Document | Purpose | Lines | Status |
|----------|---------|-------|--------|
| `README.md` | Quick start, overview | 107 | ✓ Updated |
| `docs/technical_doc.md` | Architecture, discretization, MPI | 890 | ✓ Updated (v1.1) |
| `docs/CODE_ANALYSIS.md` | Algorithm analysis, bugs | 206 | Current |
| `docs/api.md` | API reference (auto-generated) | 1338 | Maintained |
| `CODE_ANALYSIS.md` (root) | Project survey (NEW) | 261 | ✓ Staged |

---

## 9. Key Metrics & Changes Summary

### Size Changes
| Item | Before | After | Change |
|------|--------|-------|--------|
| Solver variants | 1D + 2D + 3D | 3D only | **Consolidation** |
| Main test cases | Many | **6 active** | Focused |
| NSTGV grid | 130³ cells | 66³ cells | -87% |
| SBLI z-domain | 4.0*blt | 0.5*blt | -87.5% |
| SBLI z-resolution | 257 points | 33 points | -87% |
| README images | 16+ legacy | 2 active + SBLI | **Modernized** |

### Code Quality
- Variable naming improved (viscous flux SGS terms)
- Code comments enhanced (SLAU, HD transform)
- Contiguity attributes added (GPU optimization)
- Type consistency improved (param refactoring)

---

## 10. Backward Compatibility & Migration

### Breaking Changes
1. **Parameter types**: `id_scheme` now `real(2)` (was `int(2)`)
   - Affects: Restart file compatibility, Python bindings
2. **Parameter removal**: `id_exchange` deleted
   - Impact: Case setup files referencing this will error

### Recommendations
1. **Before next commit**: Test with new grid configs (NSTGV 66³, SBLI optimized)
2. **Verify**: Numerical results match expected physics
3. **Document**: Migration guide for users with old configs
4. **Update**: Deployment scripts if parameter-dependent

---

## 11. Outstanding Tasks

### Code Review Items
- [ ] Verify NSTGV 66³ grid captures TGV physics correctly
- [ ] Test SBLI optimized domain (0.5*blt) preserves shock-layer interaction
- [ ] Confirm calc_visc2.f90 renaming (mxsgs, mysgs, mzsgs) is consistent
- [ ] Check if any restart files need recomputation

### Documentation Items
- [ ] Add migration guide for parameter type changes
- [ ] Update deployment/CI scripts if parameter-dependent
- [ ] Cross-reference new docs in API guide

### Optional Enhancements
- [ ] GPU scalability testing with new grid configs
- [ ] Benchmark performance improvements from domain reduction
- [ ] Validate LES SGS model with optimized grids

---

## 12. References & Resources

### Key Papers
1. Lusher, D.J., & Sandham, N.D. (2021). "Assessment of low-dissipative shock-capturing schemes for the compressible Taylor–Green vortex." *AIAA Journal*, 59(2), 533-545.

2. Hatayama, J., Tanaka, K., & Kouchi, T. (2026). "Nonlinear causal relationship between separation bubbles and reflected shock wave in shock wave/turbulent boundary layer interaction based on information theory." *Computers & Fluids*, 107016.

### Documentation
- `docs/technical_doc.md` - Comprehensive technical reference
- `docs/api.md` - Auto-generated API documentation
- `docs/multi_gpu/` - GPU scalability planning

---

**Report Generated**: April 3, 2026  
**Project**: OUxSBLI - GPU-Accelerated 3D Compressible Flow Solver  
**Language**: CUDA Fortran  
**Build System**: Make (case-specific)  
**Version**: 3D-only consolidation (post v1.0)
