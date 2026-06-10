# Plan: Move id_* declarations from mod_globals to mod_constant

## Context

`id_*` kind-dispatch parameters (`id_visc`, `id_scheme`, `id_accuracy`, `id_tvd`, `id_slau`, `id_rescale`, `id_gpumpi`, `id_RungeKutta`, `id_recal`) are currently declared in every per-case `mod_globals.f90` file. All callers already import them from `mod_constant` (not `mod_globals`). For 3D Cartesian solver cases, `src/mod_constant.f90.fypp` already defines all `id_*` via fypp conditionals driven by per-case `config.fypp` — the `mod_globals.f90` declarations are now dead. For 2D and 3D curvilinear solver cases, which currently lack fypp infrastructure, fypp support must be added so they can also generate per-case `mod_constant.f90` from the shared template.

---

## 1. Update `src/mod_constant.f90.fypp` — add VISC_ORDER

The 2D BL case uses 6th-order convective stencils but 2nd-order viscous stencils. The current template ties the `id_visc` value (1=2nd-order, 2=4th-order viscous) to `ORDER` (convective order). Add an optional `VISC_ORDER` variable that defaults to `ORDER` so these can differ.

**File:** `src/mod_constant.f90.fypp`

After line 1 (`#:include 'config.fypp'`), add:
```fypp
#:if not defined('VISC_ORDER')
  #:set VISC_ORDER = ORDER
#:endif
```

Then in the `id_visc` block (lines 26–34), replace every `ORDER` reference with `VISC_ORDER`:
```fypp
#:if VISC == 'NS' and VISC_ORDER >= 4
integer(4), parameter :: id_visc = 2
#:elif VISC == 'NS'
integer(4), parameter :: id_visc = 1
#:elif VISC == 'LES' and VISC_ORDER >= 4
integer(8), parameter :: id_visc = 2
#:elif VISC == 'LES'
integer(8), parameter :: id_visc = 1
```

---

## 2. Add `config.fypp` to 2D solver cases (6 new files)

Each file goes in `2D_solver/<CASE>/config.fypp`. Variables derived from the current `mod_globals.f90` kind values.

**`2D_solver/BL/config.fypp`:**
```fypp
#:set VISC          = 'NS'
#:set SCHEME        = 'KEEP'
#:set ORDER         = 6
#:set VISC_ORDER    = 2
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(BL: `id_visc=int(4)=1` → NS + 2nd-order viscous; `id_accuracy=int(8)` → 6th-order convective; `id_RungeKutta=int(2)` → TVD3)

**`2D_solver/DSL/config.fypp`:**
```fypp
#:set VISC          = 'Euler'
#:set SCHEME        = 'KEEP'
#:set ORDER         = 6
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 4
```
(DSL: `id_visc=int(2)` → Euler; `id_RungeKutta=int(4)` → RK4)

**`2D_solver/EVC/config.fypp`:**
```fypp
#:set VISC          = 'Euler'
#:set SCHEME        = 'KEEP'
#:set ORDER         = 6
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 4
```
(EVC: same as DSL; `id_RungeKutta=int(4)` → RK4)

**`2D_solver/OS/config.fypp`:**
```fypp
#:set VISC          = 'Euler'
#:set SCHEME        = 'SLAU'
#:set ORDER         = 6
#:set TVD           = 'muscl4'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(OS: `id_visc=int(2)` → Euler; `id_scheme=real(2)` → SLAU; `id_tvd=int(8)` → muscl4; `id_slau=int(4)` → HRSLAU2)

**`2D_solver/ST/config.fypp`:**
```fypp
#:set VISC          = 'NS'
#:set SCHEME        = 'KEEP'
#:set ORDER         = 6
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(ST: `id_visc=int(4)=2` → NS + 4th-order viscous; ORDER=6 ≥ 4 → value=2 ✓)

**`2D_solver/SBLI/config.fypp`:**
```fypp
#:set VISC          = 'NS'
#:set SCHEME        = 'SLAU'
#:set ORDER         = 6
#:set TVD           = 'muscl4'
#:set SLAU_VARIANT  = 'HRSLAU2'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(SBLI: `id_visc=int(4)=2` → NS + 4th-order viscous; `id_scheme=real(2)` → SLAU; `id_tvd=int(8)` → muscl4)

---

## 3. Add `config.fypp` to 3D curvilinear cases (2 new files)

**`3D_solver_curv/NACA/config.fypp`:**
```fypp
#:set VISC          = 'NS'
#:set SCHEME        = 'Hybrid'
#:set ORDER         = 2
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'SLAU'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(NACA: `id_visc=int(4)=0` → NS 2nd-order viscous; ORDER=2 < 4 → value=1; behavior identical to value=0)

**`3D_solver_curv/CORN/config.fypp`:**
```fypp
#:set VISC          = 'Euler'
#:set SCHEME        = 'Hybrid'
#:set ORDER         = 2
#:set TVD           = 'none'
#:set SLAU_VARIANT  = 'SLAU'
#:set RESCALE       = False
#:set RESTART       = False
#:set RK            = 3
```
(CORN: `id_visc=int(2)` → Euler)

---

## 4. Update 2D Makefiles — add fypp rule for mod_constant.f90 (6 files)

For each of `2D_solver/BL/Makefile`, `2D_solver/DSL/Makefile`, `2D_solver/EVC/Makefile`, `2D_solver/OS/Makefile`, `2D_solver/ST/Makefile`, `2D_solver/SBLI/Makefile`:

**Add before `.SUFFIXES`:**
```makefile
FYPP       = fypp
FYPP_FLAGS = -I$(CURDIR)

mod_constant.f90: ../../src/mod_constant.f90.fypp config.fypp
	$(FYPP) $(FYPP_FLAGS) $< $@
```

No changes needed to the `mod_constant.o: mod_globals.mod` dependency line (the generated file still `use mod_globals, only : R, gamma, Pr, dimension`).

---

## 5. Update 3D curvilinear Makefiles — add fypp rule (2 files)

For `3D_solver_curv/NACA/Makefile` and `3D_solver_curv/CORN/Makefile`:

**Add before `.SUFFIXES`:**
```makefile
FYPP       = fypp
FYPP_FLAGS = -I$(CURDIR)

mod_constant.f90: ../../src/mod_constant.f90.fypp config.fypp
	$(FYPP) $(FYPP_FLAGS) $< $@
```

---

## 6. Remove `id_*` declarations from all `mod_globals.f90` files

Remove all lines declaring: `id_visc`, `id_scheme`, `id_accuracy`, `id_tvd`, `id_slau`, `id_rescale`, `id_gpumpi`, `id_recal`, `id_RungeKutta`. Keep `dimension`, `sp`, `threshold`, `gamma`, `R`, `Pr`, `Prt`, `blt`, grid/mesh/GPU params, physical properties, etc.

Also remove the associated comment blocks (the `!!!!!!!!!!!!` tables describing these variables) only if they pertain solely to `id_*` variables. If comment blocks describe other remaining parameters, keep them.

**3D Cartesian (8 files):**
- `3D_solver/DHIT/mod_globals.f90`
- `3D_solver/ETGV/mod_globals.f90`
- `3D_solver/IVST/mod_globals.f90`
- `3D_solver/KHI/mod_globals.f90`
- `3D_solver/NSTGV/mod_globals.f90`
- `3D_solver/SBLI/mod_globals.f90`
- `3D_solver/STZ/mod_globals.f90`
- `3D_solver/TBL/mod_globals.f90`

**3D curvilinear (2 files):**
- `3D_solver_curv/NACA/mod_globals.f90`
- `3D_solver_curv/CORN/mod_globals.f90`

**2D solver (6 files):**
- `2D_solver/BL/mod_globals.f90`
- `2D_solver/DSL/mod_globals.f90`
- `2D_solver/EVC/mod_globals.f90`
- `2D_solver/OS/mod_globals.f90`
- `2D_solver/ST/mod_globals.f90`
- `2D_solver/SBLI/mod_globals.f90`

**ncu profiling (1 file — no Makefile, just remove for consistency):**
- `ncu/mod_globals.f90`

---

## 7. Update `.gitignore`

Append to `.gitignore`:
```
2D_solver/*/mod_constant.f90
3D_solver_curv/*/mod_constant.f90
```

---

## Notes on `src/mod_constant.f90`

This file (gitignored, on disk) currently imports `id_*` from `mod_globals`. Once all 2D and curvilinear cases generate their own local `mod_constant.f90` via fypp, make finds the local file first (before `../../src/mod_constant.f90` in vpath), so `src/mod_constant.f90` becomes unreachable and irrelevant. No code change needed; leave it in place.

---

## Verification

1. `cd 3D_solver/DHIT && make clean && make` — should build without errors; `mod_globals.f90` no longer defines any `id_*`.
2. `cd 2D_solver/BL && make clean && make` — fypp generates `mod_constant.f90`; verify `id_visc = integer(4), value=1` in the generated file.
3. `cd 2D_solver/EVC && make clean && make` — fypp generates `mod_constant.f90`; verify `id_RungeKutta = integer(4)` (RK4).
4. `cd 3D_solver_curv/NACA && make clean && make` — fypp generates `mod_constant.f90`.
5. `pytest ouxsbli/tests/` — all integration tests pass.
6. `grep -r "id_visc\|id_scheme\|id_accuracy\|id_tvd\|id_slau\|id_rescale\|id_gpumpi\|id_recal\|id_RungeKutta" */*/mod_globals.f90 --include="*.f90"` — no matches (all removed).
