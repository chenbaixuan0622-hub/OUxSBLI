# Specification: Fypp Macro Refactoring for `3D_solver/src/`

## 1. Context & Motivation

The original `calc_flux_base.f90` spanned 678 lines of deeply nested, dispatch-driven code controlled by five compile-time sentinel parameters (`id_visc`, `id_scheme`, `id_accuracy`, `id_bc_x/y/z`). Even though Fortran evaluates expressions like `if (kind(id_accuracy) == 8)` at compile time, the source code visually mimics runtime branching. This design forced the maintenance of four nearly identical convective scheme subroutines and three nearly identical viscous subroutines, which historically led to latent bugs (e.g., incorrect `==8` vs `>=4` conditions and missing boundary conditions in `koff` paths).

By introducing **Fypp** as a preprocessor, we transition to a **single source of truth** architecture. A localized `config.fypp` per case dictates the macro flags, which then automatically generate a flat, case-specific `calc_flux_base.f90` (~150 lines instead of 678). This approach completely eliminates dead-code branches, eradicates subroutine duplication, completely removes the old legacy `id_*` constants across the entire codebase, and implements highly scannable macro loops for directional physics ($x, y, z$).

---

## 2. Refactoring Scope & File Map

**Scope:** `3D_solver/` only. `3D_solver_curv/` and `2D_solver/` are out of scope for this refactoring.

Every file leveraging Fypp features must use the `*.fypp` extension (or `config.fypp` for configurations). The legacy `id_*` runtime-like parameters are completely eliminated from the codebase, including the GPU kernel definitions.

| File Path | Action | Description |
|:---|:---|:---|
| `3D_solver/*/config.fypp` | **Create** | New per-case macro configuration file. |
| `3D_solver/*/mod_globals.f90` | **Rename → *.f90.fypp | Includes `config.fypp`. Completely removes `id_*` constants. |
| `3D_solver/src/calc_flux_base.f90` | **Rename → *.f90.fypp | Includes `calc_EFG`, `calc_EFG_halo`, `calc_conv_EF`, `calc_conv_G_koff` wrappers, all generated via VISC/SCHEME macros. |
| `3D_solver/src/preprocess.f90` | **Rename → *.f90.fypp | Replaces old compile-time memory allocations with macro blocks. |
| `3D_solver/src/calc_*_kernel.f90` | **Rename → *.f90.fypp | Removes `id_accuracy` arguments; internalizes order with macros. Integrate subroutines for x, y, and z direction.|
| `3D_solver/*/Makefile` | **Update** | Implements `.f90.fypp` pattern matching rules and dependencies. |
| `3D_solver/src/calc_time_dev.f90` | **Edit** | Removes `id_visc` from `calc_EFG` / `calc_EFG_halo` calls. |

---

### Zdec/Halo-Exchange Architecture Note

Commit d44763c (2026-05-24) introduced a **two-phase flux dispatch** within the z-direction halo update cycle:

1. **Full-domain path** (`k_lo_G == 1 .and. k_hi_G == nz - 1`): Convective and viscous fluxes computed over all interior points; `calc_EFG` dispatches uniformly.

2. **Zdec interior path** (restricted k-range): Interior z-slices skip halo computation; the convective flux dispatch branches into:
   - `calc_conv_EF` wrapper: x + y convective flux (full domain)
   - `calc_conv_G_koff` wrapper: z convective flux restricted to `[k_lo_G, k_hi_G]`
   
   Viscous fluxes use `_koff` variants with explicit k-range arguments.

The Fypp templates in `calc_flux_base.f90.fypp` (Step 5) must generate both the unified `calc_EFG` interface and these specialized wrappers.

## 3. Fypp Variable Mapping

| Fypp Variable | Meaning | Allowed Values | Replaces Legacy Mechanism |
|:---|:---|:---|:---|
| `VISC` | Viscosity model | `'Euler'`, `'NS'`, `'LES'` | `kind(id_visc)` (2 / 4 / 8) |
| `SCHEME` | Convective scheme | `'KEEP'`, `'SLAU'`, `'Roe'`, `'Hybrid'` | Data type of `id_scheme` |
| `ORDER` | Spatial accuracy order | `2`, `4`, `6` | `kind(id_accuracy)` (2 / 4 / 8) |
| `ORDER_IO` | Stencil half-width | `ORDER // 2 - 1` ($0, 1, 2$) | Derived from `kind(id_accuracy)/3` |
| `BC_X` | Non-periodic $x$-boundary? | `True` / `False` | Logical parameter `id_bc_x` |
| `BC_Y` | Non-periodic $y$-boundary? | `True` / `False` | Logical parameter `id_bc_y` |
| `BC_Z` | Non-periodic $z$-boundary? | `True` / `False` | Logical parameter `id_bc_z` |
| `RESCALE` | Inflow rescaling status | `True` / `False` | `kind(id_rescale)` (4 $\rightarrow$ ON, 2 $\rightarrow$ OFF) |
| `RK` | Runge-Kutta order | `3`, `4` | `kind(id_RungeKutta)` (2 $\rightarrow$ RK3, 4 $\rightarrow$ RK4) |
| `RESTART` | Execution mode | `True` / `False` | `kind(id_recal)` (4 $\rightarrow$ Restart, 2 $\rightarrow$ Init) |
| `TVD` | TVD Limiter selection | `'none'`, `'minmod'`, `'muscl4'`| `kind(id_tvd)` (2 / 4 / 8) |
| `SLAU_VARIANT`| SLAU sub-type variant | `'SLAU'`, `'HRSLAU2'` | `kind(id_slau)` (2 $\rightarrow$ SLAU, 4 $\rightarrow$ HRSLAU2) |

---

## 4. Implementation Steps

### Step 1: Activate venv
Ensure Fypp is accessible in your build environment:
```bash
source ~/htymenv_new/bin/activate
fypp --version  # Verify installation (3.x expected)

```

### Step 2: Configure Per-Case `config.fypp`

Create a `config.fypp` within each of the 8 case directories (NSTGV, ETGV, IVST, KHI, SBLI, TBL, STZ, DHIT) matching its physical profile. Example for `3D_solver/SBLI/config.fypp`:

```fortran
! 3D_solver/SBLI/config.fypp
#:set VISC         = 'NS'      ! Navier-Stokes simulation
#:set SCHEME       = 'SLAU'    ! SLAU convective flux
#:set ORDER        = 6         ! 6th-order spatial accuracy
#:set BC_X         = True      ! Wall/Inflow/Outflow conditions
#:set BC_Y         = True
#:set BC_Z         = False     ! Periodic in spanwise z
#:set RESCALE      = True      ! Inflow rescaling enabled
#:set RK           = 3         ! TVD-RK3 time advancement
#:set RESTART      = False     ! Initialize fields from scratch
#:set TVD          = 'none'    ! No limiters needed
#:set SLAU_VARIANT = 'HRSLAU2' ! High-resolution SLAU2 variant
#:set ORDER_IO     = ORDER // 2 - 1

```

For the **STZ case** (smaller domain, 2nd-order max):

```fortran
! 3D_solver/STZ/config.fypp
#:set VISC         = 'Euler'    ! Euler simulation
#:set SCHEME       = 'SLAU'     ! SLAU convective flux
#:set ORDER        = 2          ! 2nd-order (6th-order unsupported for this domain size)
#:set BC_X         = False      ! Periodic in x
#:set BC_Y         = False      ! Periodic in y
#:set BC_Z         = True       ! Shock tube; nonperiodic z-boundaries
#:set RESCALE      = False      ! No rescaling needed
#:set RK           = 3          ! TVD-RK3 time advancement
#:set RESTART      = False      ! Initialize from scratch
#:set TVD          = 'muscl4'   ! MUSCL 4th-order limiter
#:set SLAU_VARIANT = 'SLAU'     ! Standard SLAU variant
#:set ORDER_IO     = ORDER // 2 - 1

```

### Step 3: Purge Constants from `mod_globals.f90.fypp`

Rename `mod_globals.f90` $\rightarrow$ `mod_globals.f90.fypp` and strip out all legacy `id_*` variables. Only regular grid dimensions, physics metrics, and GPU block configurations remain.

```fortran
#:include 'config.fypp'
module mod_globals
  use cudafor
  implicit none

  ! ── Legacy id_* parameters have been fully purged ──
  ! All compile-time logic is now handled via config.fypp macro flags.

  integer, parameter :: sp = kind(1.d0)
  real(sp), parameter :: threshold = 0.4_sp
  ! ... (Keep grid variables, thread configurations, and physical coefficients)
end module mod_globals

```

### Step 4: Streamline GPU Kernels (`calc_*_kernel.f90.fypp`)

By processing the kernel files through Fypp, we drop the `id_accuracy` interface parameter completely. Subroutines optimize purely for the targeted compilation order.

```fortran
#:include 'config.fypp'
module calc_keep_kernel
  implicit none
contains
  attributes(global) subroutine calc_keep_x(nx, ny, nz, Q, T, E)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), device :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)

    ! Only the required spatial scheme loop is generated—no generic type dispatch.
    #:if ORDER == 6
      ! Dedicated 6th-order numerical stencil computations
    #:elif ORDER == 4
      ! Dedicated 4th-order numerical stencil computations
    #:else
      ! Dedicated 2nd-order numerical stencil computations
    #:endif
  end subroutine calc_keep_x
end module calc_keep_kernel

```

### Step 5: Macro Loop Integration in `calc_flux_base.f90.fypp`

Instead of duplicating structural patterns across $x, y,$ and $z$, we consolidate code execution using a directional dictionary array (`DIRS`) inside a macro `#:for` loop.

```fortran
#:include 'config.fypp'
!> Generated from calc_flux_base.f90.fypp: VISC=${VISC}$, SCHEME=${SCHEME}$, ORDER=${ORDER}$

module calc_flux_base
  use mod_globals, only : sp, blocks, threads, &
  & blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  
  ! Dynamic module loading based on case settings
  #:if SCHEME == 'KEEP'
    use calc_keep_kernel
  #:elif SCHEME == 'SLAU'
    use calc_slau_kernel
  #:elif SCHEME == 'Hybrid'
    use calc_hybrid_kernel
    use calc_hybrid
  #:endif

  implicit none
  private
  public calc_EFG, calc_EFG_halo, init_sensor

  integer, parameter :: io_fb      = ${ORDER_IO}$
  integer, parameter :: overlap_fb = ${ORDER_IO}$ + 1
  
  ! ── Directional Metadata Mapping ──────────────────────────────────────────
  #:set DIRS = [ &
    {'d': 'x', 'f': 'E', 'bc': BC_X, 'bl': 'blocksE', 'th': 'threadsE', 'inv': 'inv_dx, inv_dy, inv_dz', 'k_lo_off': 'overlap_fb+2'}, &
    {'d': 'y', 'f': 'F', 'bc': BC_Y, 'bl': 'blocksF', 'th': 'threadsF', 'inv': 'inv_dy, inv_dx, inv_dz', 'k_lo_off': 'overlap_fb+2'}, &
    {'d': 'z', 'f': 'G', 'bc': BC_Z, 'bl': 'blocksG', 'th': 'threadsG', 'inv': 'inv_dx, inv_dy, inv_dz', 'k_lo_off': 'overlap_fb+1'}  &
  ]
  ! ──────────────────────────────────────────────────────────────────────────

contains

  ! ... (init_sensor and launch_Ducros declarations)

  ! ── CORE DISPATCH: calc_EFG (generated per VISC model) ────────────────────
  #:for VISC_VAL in ['Euler', 'NS', 'LES']
  #:if VISC == VISC_VAL
  subroutine calc_EFG_${VISC_VAL}$(nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G, k_lo_G, k_hi_G)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in) :: inv_dx, inv_dy, inv_dz
    real(8), intent(in), device :: Jacobian(nx,ny,nz), QJ(5,nx,ny,nz), Q(nx,5,ny,nz), T(nx,ny,nz)
    #:if VISC_VAL == 'Euler'
    real(8), intent(in), device :: mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    #:elif VISC_VAL == 'NS'
    real(8), intent(in), device :: mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1)
    #:else
    real(8), intent(in), device :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    #:endif
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1)
    integer, intent(in), value :: k_lo_G, k_hi_G

    if (k_lo_G == 1 .and. k_hi_G == nz - 1) then
      ! ─── Full Domain Execution Phase ──────────────────────────────────────
#:if SCHEME in ('SLAU', 'Hybrid')
      call launch_Ducros(nx, ny, nz, inv_dx, inv_dy, inv_dz, Q)
#:endif

      ! Convective Flux Generation Loop
#:for dir in DIRS
  #:if dir['bc']
      call calc_${SCHEME.lower()}$_${dir['d']}$<<<${dir['bl']}$,${dir['th']}$>>>(nx, ny, nz, Q, ${('T, ' if SCHEME=='KEEP' else 'sensor, ')}$${dir['f']}$)
  #:else
      call calc_${SCHEME.lower()}$_${dir['d']}$_in<<<${dir['bl']}$,${dir['th']}$>>>(nx, ny, nz, Q, ${('T, ' if SCHEME=='KEEP' else 'sensor, ')}$${dir['f']}$)
  #:endif
#:endfor

      ! Viscous Flux Generation Loop
#:if VISC != 'Euler'
  #:for dir in DIRS
    #:set V = '' if VISC=='NS' else '_LES'
    #:set MQC = '' if VISC=='NS' else 'mut, qc2, '
    #:if ORDER >= 4
      #:if dir['bc']
      call calc_${dir['f']}$v${V}$4<<<${dir['bl']}$v,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$)
      #:else
      call calc_${dir['f']}$v${V}$4_in<<<${dir['bl']}$v,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$)
      #:endif
    #:else
      call calc_${dir['f']}$v${V}$2<<<${dir['bl']}$v,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$)
    #:endif
  #:endfor
#:endif

    else
      ! ─── Zdec Interior Execution Phase ────────────────────────────────────
      ! x and y loops compute full domain; z limits calculations to koff bounds
#:for dir in DIRS
  #:if dir['d'] != 'z'
    #:if dir['bc']
      call calc_${SCHEME.lower()}$_${dir['d']}$<<<${dir['bl']}$,${dir['th']}$>>>(nx, ny, nz, Q, ${('T, ' if SCHEME=='KEEP' else 'sensor, ')}$${dir['f']}$)
    #:else
      call calc_${SCHEME.lower()}$_${dir['d']}$_in<<<${dir['bl']}$,${dir['th']}$>>>(nx, ny, nz, Q, ${('T, ' if SCHEME=='KEEP' else 'sensor, ')}$${dir['f']}$)
    #:endif
  #:else
      block
        integer :: n_G
        type(dim3) :: bG
        n_G = k_hi_G - k_lo_G + 1
        bG = dim3(blocksG%x, blocksG%y, (n_G + threadsG%z - 1) / threadsG%z)
        call calc_${SCHEME.lower()}$_z_in_koff<<<bG,threadsG>>>(nx, ny, nz, Q, ${('T, ' if SCHEME=='KEEP' else 'sensor, ')}$G, k_lo_G, k_hi_G)
      end block
  #:endif
#:endfor

      ! Viscous Interior koff Tracking
#:if VISC != 'Euler'
  #:for dir in DIRS
    #:set V = '' if VISC=='NS' else '_LES'
    #:set MQC = '' if VISC=='NS' else 'mut, qc2, '
    #:set BL_INT = 'blocksEv_int' if dir['d'] in ('x', 'y') else 'blocksGv_int'
    #:if ORDER >= 4
      #:set SUFFIX = '_koff' if dir['bc'] else '_in_koff'
      call calc_${dir['f']}$v${V}$4${SUFFIX}$<<<${BL_INT}$,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, ${dir['k_lo_off']}$, nz-overlap_fb-1)
    #:else
      call calc_${dir['f']}$v${V}$2_koff<<<${BL_INT}$,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, ${dir['k_lo_off']}$, nz-overlap_fb-1)
    #:endif
  #:endfor
#:endif
    end if
  end subroutine calc_EFG_${VISC_VAL}$
  #:endif
  #:endfor

  ! Single-procedure interface: only calc_EFG_${VISC}$ was generated above.
  interface calc_EFG
    module procedure calc_EFG_${VISC}$
  end interface calc_EFG

  ! ── HALO DISPATCH: calc_EFG_halo (generated per VISC model) ───────────────
  ! Same #:for VISC_VAL pattern as calc_EFG — mu/mut/qc2 dimensions are VISC-specific.
  !
  ! Halo flux computation ranges (in terms of overlap_fb = ORDER_IO + 1):
  !
  !  | Flux Type    | Low halo range                  | High halo range                      |
  !  |---|---|---|
  !  | Convective G | [overlap_fb : 2*overlap_fb-1]   | [nz-2*overlap_fb+1 : nz-overlap_fb]  |
  !  | Viscous E, F | [2 : overlap_fb+1]              | [nz-overlap_fb : nz-1]               |
  !  | Viscous G    | [overlap_fb : 2*overlap_fb]     | [nz-2*overlap_fb : nz-overlap_fb]    |

  #:for VISC_VAL in ['Euler', 'NS', 'LES']
  #:if VISC == VISC_VAL
  subroutine calc_EFG_halo_${VISC_VAL}$(nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in) :: inv_dx, inv_dy, inv_dz
    real(8), intent(in), device :: Jacobian(nx,ny,nz), QJ(5,nx,ny,nz), Q(nx,5,ny,nz), T(nx,ny,nz)
    #:if VISC_VAL == 'Euler'
    real(8), intent(in), device :: mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    #:elif VISC_VAL == 'NS'
    real(8), intent(in), device :: mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1)
    #:else
    real(8), intent(in), device :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    #:endif
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1)
    ! ... (refresh ghost Q at z-halos, then compute convective G and viscous E/F/G
    !      in the koff ranges tabulated above)
#:if VISC != 'Euler'
  #:for dir in DIRS
    #:set V = '' if VISC=='NS' else '_LES'
    #:set MQC = '' if VISC=='NS' else 'mut, qc2, '
    #:if ORDER >= 4
      #:set SUFFIX = '_koff' if dir['bc'] else '_in_koff'
      call calc_${dir['f']}$v${V}$4${SUFFIX}$<<<blocks${dir['f']}$v_halo,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, ${'2' if dir['d'] in ('x','y') else '1'}$, ${dir['k_lo_off']}$)
      call calc_${dir['f']}$v${V}$4${SUFFIX}$<<<blocks${dir['f']}$v_halo,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, nz-overlap_fb, nz-1)
    #:else
      call calc_${dir['f']}$v${V}$2_koff<<<blocks${dir['f']}$v_halo,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, ${'2' if dir['d'] in ('x','y') else '1'}$, ${dir['k_lo_off']}$)
      call calc_${dir['f']}$v${V}$2_koff<<<blocks${dir['f']}$v_halo,${dir['th']}$v>>>(nx, ny, nz, ${dir['inv']}$, Q, T, mu, ${MQC}$${dir['f']}$, nz-overlap_fb, nz-1)
    #:endif
  #:endfor
#:endif
  end subroutine calc_EFG_halo_${VISC_VAL}$
  #:endif
  #:endfor

  interface calc_EFG_halo
    module procedure calc_EFG_halo_${VISC}$
  end interface calc_EFG_halo
end module calc_flux_base

```

### Step 6: Refactor `preprocess.f90.fypp`

Simplify the dynamic allocation arrays and setup routines by dropping runtime checks.

```fortran
! allocate_device_mem macro segments
#:if VISC == 'Euler'
  allocate(T(nx,ny,nz), mu(1,1,1), mut(1,1,1), qc2(1,1,1), stat=ierr)
#:elif VISC == 'NS'
  allocate(T(nx,ny,nz), mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1), stat=ierr)
#:else
  allocate(T(nx,ny,nz), mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz), stat=ierr)
#:endif

! pre_calc segment
  overlap = ${ORDER_IO + 1}$

```

### Step 7: Update Build Infrastructure (Makefile Pattern Rules)

Add `.f90.fypp → .f90` pattern matching rules to each case's `Makefile` to automatically trigger Fypp preprocessing. This preserves the existing `mpif90 -cuda -acc -fast -gpu=ptxinfo,rdc,lto` flags and dependency chain.

Add these rules to each `3D_solver/<CASE>/Makefile`:

```makefile
# Fypp preprocessing rules
FYPP = fypp
FYPP_FLAGS = -I$(CURDIR)

# Rule for generating .f90 from .f90.fypp in the current case directory
%.f90: %.f90.fypp config.fypp
	$(FYPP) $(FYPP_FLAGS) $< $@

# Rule for generating .f90 from .f90.fypp in the shared src directory
../src/%.f90: ../src/%.f90.fypp config.fypp
	$(FYPP) $(FYPP_FLAGS) -I$(CURDIR) $< $@
```

Example: To build `3D_solver/SBLI`, simply run:

```bash
cd 3D_solver/SBLI
make clean && make -j
```

The Makefile will automatically preprocess `.f90.fypp` files (both in the case directory and in `../src/`) to `.f90` before compilation, with dependencies on `config.fypp` ensuring regeneration when configuration changes.


### Step 8: Call-Site Adjustments (`calc_time_dev.f90`)

Strip out the legacy `id_visc` tracker from all calling locations:

```fortran
! Old Call Syntax:
! call calc_EFG(id_visc, nx, ny, nz, ...)

! New Standard Call Syntax:
call calc_EFG(nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G, k_lo_G, k_hi_G)
call calc_EFG_halo(nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)

```

---

## 5. Verification & Testing

Verify correctness step-by-step for a specific targeted configuration (e.g., `SBLI`):

```bash
cd 3D_solver/SBLI
make clean

# Test isolated manual parsing before a full compilation
fypp -I. mod_globals.f90.fypp mod_globals.f90
fypp -I. ../src/calc_flux_base.f90.fypp calc_flux_base.f90

# Confirm structural checks are clean
grep -n 'id_visc' calc_flux_base.f90  # Should return empty
grep -n 'kind(id' calc_flux_base.f90  # Should return empty

# Verify ORDER_IO formula is correct
grep -n 'ORDER_IO' config.fypp       # Should read: ORDER // 2 - 1

# Build and run tests
make -j
bash calc.sh

# Test the zdec/halo-exchange path with multiple MPI ranks
mpirun -n 2 a.out                    # Exercises interior/halo phase transition

# Verify halo boundary fluxes via ParaView inspection or Fortran verification routine
# (Check that the computed z-flux results match single-rank baseline where applicable)

```

