# Plan: Fix Bugs + Move id_* Backward-Compat Layer to `mod_constant.f90.fypp`

## Context

Two tasks:
1. **Fix bugs** in staged files (`calc_slau_kernel.f90.fypp`, `calc_visc4.f90.fypp`, `calc_visc4_internal.f90.fypp`, `src/mod_constant.f90`).
2. **Move the backward-compat `id_*` block** out of `mod_globals.f90.fypp` into a new `src/mod_constant.f90.fypp`. All files that currently `use mod_globals, only : id_*` must switch to `use mod_constant, only : id_*`.

---

## Part 1 — Bug Fixes

### Bug 1 — `calc_slau_kernel.f90.fypp:272`: typo `fnz` → `fdz`

In `calc_slau_z`, ORDER==4 branch, the `w` component uses undefined `fnz`:
```fortran
! Wrong:
call delta4(fnz,   w(idx-1:idx+2),   wl,   wr(idx_r))
! Fix:
call delta4(fdz,   w(idx-1:idx+2),   wl,   wr(idx_r))
```

### Bug 2 — `calc_visc4.f90.fypp:55` and `calc_visc4_internal.f90.fypp:55`: truncated `intege` line

Both files share the same stray `intege` on line 55 inside the `#:if is_koff` block of the x-direction subroutine. Delete it in both files.

### Bug 3 — `calc_visc4.f90.fypp`: missing `#:endfor` + wrong `end subroutine` names

The file has `#:for is_koff in [False, True]` at lines 25, 205, 384 with **no `#:endfor` anywhere** — Fypp treats them as nested, generating 2+4+8 copies instead of 2+2+2. Additionally the `end subroutine` names are missing `${suf}$`.

| Line | Current | Fix |
|------|---------|-----|
| 200 | `end subroutine calc_Ev4` | `end subroutine calc_Ev4${suf}$` then `#:endfor` |
| 379 | `end subroutine calc_Fv4` | `end subroutine calc_Fv4${suf}$` then `#:endfor` |
| 554 | `end subroutine calc_Gv4` | `end subroutine calc_Gv4${suf}$` then `#:endfor` |

### Bug 4 — `calc_visc4_internal.f90.fypp`: same missing `#:endfor` + wrong names (including missing `_in`)

| Line | Current | Fix |
|------|---------|-----|
| 129 | `end subroutine calc_Ev4_in` | `end subroutine calc_Ev4${suf}$_in` then `#:endfor` |
| 237 | `end subroutine calc_Fv4` | `end subroutine calc_Fv4${suf}$_in` then `#:endfor` |
| 343 | `end subroutine calc_Gv4` | `end subroutine calc_Gv4${suf}$_in` then `#:endfor` |

---

## Part 2 — Move id_* to `mod_constant.f90.fypp`

### Step 1 — Create `src/mod_constant.f90.fypp`

Convert the existing `src/mod_constant.f90` to a Fypp file. The new file:
- Adds `#:include 'config.fypp'` at the top
- Removes the `use mod_globals, only : ..., id_visc, id_accuracy, ...` lines for all id_*
- Keeps `use mod_globals, only : R, gamma, Pr, dimension` (needed for Cp, gamma_1, and `hoge`)
- Defines all id_* using Fypp macros (moved from `mod_globals.f90.fypp`):

```fortran
#:include 'config.fypp'
module mod_constant
  use cudafor
  use mod_globals, only : R, gamma, Pr, dimension
  implicit none
  ! physical constants (unchanged)
  real(8), parameter :: gamma_1        = gamma - 1.d0
  ...
  ! Fypp-generated kind-dispatch parameters
  #:if VISC == 'Euler'
  integer(2), parameter :: id_visc = 1
  #:elif VISC == 'NS' and ORDER >= 4
  integer(4), parameter :: id_visc = 2
  ...
  #:endif
  #:if ORDER == 2
  integer(2), parameter :: id_accuracy = 0
  ...
  #:endif
  #:if TVD == 'none'
  integer(2), parameter :: id_tvd = 0
  ...
  #:endif
  #:if SLAU_VARIANT == 'SLAU'
  integer(2), parameter :: id_slau = 0
  #:elif SLAU_VARIANT == 'HRSLAU2'
  integer(4), parameter :: id_slau = 0
  #:endif
  integer(2), parameter :: id_gpumpi = 0
  #:if RESCALE
  integer(4), parameter :: id_rescale = 0
  #:else
  integer(2), parameter :: id_rescale = 0
  #:endif
  #:if RESTART
  integer(4), parameter :: id_recal = 0
  #:else
  integer(2), parameter :: id_recal = 0
  #:endif
  #:if RK == 3
  integer(2), parameter :: id_RungeKutta = 0
  #:elif RK == 4
  integer(4), parameter :: id_RungeKutta = 0
  #:elif RK == 8
  integer(8), parameter :: id_RungeKutta = 0
  #:endif
  #:if SCHEME == 'KEEP'
  integer(2), parameter :: id_scheme = 0
  #:elif SCHEME == 'SLAU'
  real(2),    parameter :: id_scheme = 0
  #:elif SCHEME == 'Roe'
  real(4),    parameter :: id_scheme = 0
  #:elif SCHEME == 'Hybrid'
  real(8),    parameter :: id_scheme = 0
  #:endif
contains
  subroutine check()
    integer :: hoge = dimension + id_visc &
                      + id_accuracy + id_tvd + id_slau &
                      + id_rescale + id_gpumpi + id_RungeKutta + id_recal
  end subroutine check
end module mod_constant
```

Note: `id_scheme` is re-introduced here (for `3D_solver_curv` which still uses runtime dispatch). It is NOT re-added to `mod_globals.f90.fypp`.

### Step 2 — Remove the backward-compat block from all 8 `mod_globals.f90.fypp` files

Remove the entire section between the comments `! Fypp-generated kind-dispatch parameters (backward compatibility layer)` and the end of the `id_RungeKutta` block (inclusive). Keep `id_bc_x/y/z` (logical boundary flags — these are NOT id_* dispatch parameters and stay in mod_globals).

Affected files (8 total):
- `3D_solver/NSTGV/mod_globals.f90.fypp`
- `3D_solver/ETGV/mod_globals.f90.fypp`
- `3D_solver/IVST/mod_globals.f90.fypp`
- `3D_solver/KHI/mod_globals.f90.fypp`
- `3D_solver/DHIT/mod_globals.f90.fypp`
- `3D_solver/TBL/mod_globals.f90.fypp`
- `3D_solver/SBLI/mod_globals.f90.fypp`
- `3D_solver/STZ/mod_globals.f90.fypp`

### Step 3 — Update all callers: `use mod_globals → use mod_constant` for id_*

Each file below imports one or more id_* from `mod_globals`; move those to `use mod_constant`:

**Fypp source files (`.f90.fypp`):**
| File | Variables to move |
|------|-------------------|
| `3D_solver/src/calc_visc4.f90.fypp:6` | `id_visc` |
| `3D_solver/src/calc_visc4_internal.f90.fypp:6` | `id_visc` |
| `3D_solver/src/calc_slau_kernel.f90.fypp:3` | `id_slau` |
| `3D_solver/src/calc_slau_kernel_internal.f90.fypp:3` | `id_slau` |
| `3D_solver/src/calc_hybrid_kernel.f90.fypp:3` | `id_slau` |
| `3D_solver/src/calc_hybrid_kernel_internal.f90.fypp:3` | `id_slau` |

**Permanent `.f90` files:**
| File | Variables to move |
|------|-------------------|
| `src/print.f90:4` | `id_accuracy` |
| `src/calc_muscl.f90:143,156` | `id_tvd` |
| `src/sbli.f90:5` | `id_RungeKutta, id_rescale, id_recal` |
| `3D_solver/src/calc_rescale.f90:4` | `id_gpumpi, id_recal` |
| `3D_solver/src/calc_visc4_les_internal.f90:5` | `id_visc` |
| `3D_solver/src/calc_roe_kernel_internal.f90:7` | `id_accuracy` |
| `3D_solver/NSTGV/set.f90:8,17,51` | `id_accuracy` |
| `3D_solver/DHIT/set.f90:9,18,28` | `id_accuracy` |
| `3D_solver/DHIT/set_init_dhit.f90:3` | `id_accuracy` |
| `3D_solver/STZ/set.f90:10,69` | `id_accuracy` |
| `3D_solver/TBL/set.f90:4` | `id_rescale` |
| `3D_solver/SBLI/set.f90:4` | `id_rescale` |
| `3D_solver_curv/src/calc_flux_base_curv.f90:5` | `id_scheme, id_accuracy, id_slau` |
| `3D_solver_curv/src/calc_slau_kernel_curv.f90:4` | `id_slau` |
| `3D_solver_curv/src/calc_hybrid_kernel_curv.f90:4` | `id_accuracy, id_slau` |
| `3D_solver_curv/src/main_curv.f90:4` | `id_RungeKutta, id_rescale, id_recal` |
| `3D_solver_curv/src/calc_time_dev_curv.f90:5` | `id_visc` |
| `3D_solver_curv/src/preprocess_curv.f90:31,81` | `id_visc, id_accuracy` |
| `2D_solver/EVC/set.f90:2` | `id_accuracy` |
| `2D_solver/src/calc_hybrid_kernel_internal.f90:9` | `id_accuracy, id_slau` |

Also check `main.f90` files (STZ and other case-specific) for `id_RungeKutta`, `id_rescale`, `id_recal` imports from `mod_globals`.

### Step 4 — Update Makefiles to generate `mod_constant.f90` from fypp

Add a rule to cover `../../src/` files (for all 3D_solver, 3D_solver_curv, 2D_solver case Makefiles):

```makefile
../../src/%.f90: ../../src/%.f90.fypp config.fypp
	$(FYPP) $(FYPP_FLAGS) $< $@
```

The existing `../src/%.f90` rule (already in place) covers `3D_solver/src/` files. This new rule covers the shared `src/` directory. Add to all case Makefiles (8 for 3D_solver, plus curv and 2D cases).

---

## Verification

```bash
# Check id_* are gone from mod_globals, present in mod_constant
cd 3D_solver/NSTGV
fypp -I. mod_globals.f90.fypp /dev/stdout | grep "id_visc\|id_accuracy\|id_tvd"  # must be empty
fypp -I. ../../src/mod_constant.f90.fypp /dev/stdout | grep "integer.*id_visc"    # must show definition

# Check visc4 generates correct subroutine names (x/y/z + koff variants)
fypp -I. ../src/calc_visc4.f90.fypp /dev/stdout | grep "subroutine calc_"
# expected: calc_Ev4, calc_Ev4_koff, calc_Fv4, calc_Fv4_koff, calc_Gv4, calc_Gv4_koff

# Compile
make clean && make -j

# For a case with rescaling (TBL uses id_rescale)
cd ../TBL && make clean && make -j
```
