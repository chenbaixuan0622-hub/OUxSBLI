# Bug-Fix Plan for OUxSBLI

## Context

The user asked for a thorough bug search across the Cartesian solver (`3D_solver/`) and the new curvilinear solver (`3D_solver_curv/`), plus updates to `CLAUDE.md` and new docs if needed. Three parallel Explore agents produced initial candidate bugs; each was independently verified by reading the actual source files. Two bugs were confirmed; several agent claims were disproved.

---

## Bug 1 — CRITICAL: Missing `/dz` in SGS heat flux, z-direction (LES only)

**File:** `3D_solver_curv/src/calc_visc2_curv.f90`, line 711  
**Subroutine:** `calc_Gv_LES2_curv`

### Root cause
In the curvilinear solver the ξ and η computational directions use unit spacing (Δξ = Δη = 1), so finite differences in those directions need no dimensional scaling. The z-direction, however, uses the physical spacing `dz` that is passed in as an argument. The molecular heat-flux term on the same z-face (line 736) correctly divides by `dz`:

```fortran
viscous_work = Cp_over_Pr * mu_f * dTdz / dz + ...
```

But the SGS (subgrid-scale) heat-flux at line 711 omits the `/dz`:

```fortran
! CURRENT (buggy)
Hsgs = -0.5d0*(mut(i,j,k)+mut(i,j,k+1))*(H2-H1)/Prt

! CORRECT
Hsgs = -0.5d0*(mut(i,j,k)+mut(i,j,k+1))*(H2-H1)/dz/Prt
```

`(H2-H1)` is a dimensionless enthalpy difference (computational units); dividing by `dz` converts it to a physical gradient. Without the division the SGS heat flux is `dz` times too large.

Countercheck: the same subroutine's stress-tensor SGS terms all correctly divide by `dz` (e.g. `muzsgs = mut_f_z*(u(idx+1)-u(idx))/dz` at line 723).  
The analogous ξ- and η-direction Hsgs in `calc_Ev_LES2_curv` (line ≈401) and `calc_Fv_LES2_curv` (line ≈565) are correct because Δξ = Δη = 1.

### Fix
```fortran
! line 711 – add /dz before /Prt
Hsgs = -0.5d0*(mut(i,j,k)+mut(i,j,k+1))*(H2-H1)/dz/Prt
```

---

## Bug 2 — MINOR: Block-local variables written uninitialized to shared memory (hybrid kernel, KEEP branch)

**File:** `3D_solver/src/calc_hybrid_kernel.f90`, lines 132–137 (and the analogous y/z kernels)  
**Subroutines:** `calc_hybrid_x6`, `calc_hybrid_x4`, `calc_hybrid_y6`, `calc_hybrid_y4`, `calc_hybrid_z6`, `calc_hybrid_z4`

### Root cause
Inside the `block` construct (line 71–138), `rhol, ul, vl, wl, pl` are declared block-local. They are only assigned in the SLAU branch (`fdx > threshold`). In the KEEP branch these variables remain uninitialized. After `syncthreads()` at line 132, ALL threads—including KEEP-branch threads—write their (possibly uninitialized) `rhol, ul, vl, wl, pl` to the shared-memory arrays at `(it,jt,kt)`:

```fortran
call syncthreads()
rho(it,jt,kt) = rhol   ! uninitialized for KEEP-branch threads
  u(it,jt,kt) = ul
  ...
```

### Why no computation error occurs
After the write, SLAU threads at position `(it,jt,kt)` read back only their OWN position from shared memory (line 141), and KEEP threads do not enter the SLAU block (line 139 guard). Neighboring SLAU threads already captured their right-state scalars (`rho_r`, etc.) BEFORE `syncthreads()`. Therefore no SLAU thread ever reads garbage from a KEEP thread's position.

Despite being harmless at runtime, this is undefined behavior per the Fortran standard.

### Fix (optional — low priority)
Initialize the block-local variables before the `if` chain, or assign them explicitly inside the KEEP branch:

```fortran
! Add after line 72 in each affected kernel:
rhol = 0.d0; ul = 0.d0; vl = 0.d0; wl = 0.d0; pl = 0.d0
```

---

## Disproved agent claims (for reference)

| Claim | Verdict |
|---|---|
| `calc_keep_x2`: passing `u` as 6th KEEP arg is wrong | Disproved. KEEP's 6th arg is normal velocity `uu`; for x-direction, `uu = u`. Higher-order variants use `associate(uu => u)`. |
| `calc_Ev2`: `(my1-my2)` formula is physically wrong | Disproved. Algebraic expansion shows equivalence to the explicit LES form. |
| `calc_Ev_LES2`: `Q(i,2)**2` is ρ²u² | Disproved. `calc_quantities_T_3D` is called first, converting QJ → Q (primitive), so Q(i,2) = u. |

---

## CLAUDE.md Updates

Add a new section documenting:

1. **`3D_solver_curv/`** — curvilinear 2D O-grid solver; per-case config in `3D_solver_curv/<CASE>/`; uses `main_curv.f90` as entry point; Jacobian convention: `Jacobian(i,j) = 1/(J_2D·dz)`, stored quantity is `QJ = Q_physical × J_2D × dz`.

2. **`load_smem_visc2.f90`** — new untracked file in `3D_solver/src/` providing async pipeline shared-memory load helpers (`pipelineMemcpyAsync`/`pipelineCommit`/`pipelineWaitPrior`). Requires `wmma` module.

3. **LES limitation**: `id_visc = kind8` (LES) is NOT yet implemented in the curvilinear solver (`3D_solver_curv/`). Only Euler (`integer(2)`) and NS (`integer(4)`) are supported via `calc_flux_base_curv.f90`.

4. **`id_scheme` kind table correction**: fix the comment in `mod_globals.f90` — `integer(2)` maps to KEEP; `real(2)` maps to SLAU; `real(8)` maps to Hybrid.

---

## docs/ File to Create

**`docs/bugs.md`** — document the two bugs above (critical and minor), with file paths, line numbers, root cause, and the fix.

---

## Execution Order

1. Fix Bug 1: edit `3D_solver_curv/src/calc_visc2_curv.f90` line 711.
2. Fix Bug 2 (optional): add initializations in `3D_solver/src/calc_hybrid_kernel.f90`.
3. Update `CLAUDE.md`.
4. Create `docs/bugs.md`.

---

## Verification

- Rebuild the curvilinear NACA case: `cd 3D_solver_curv/NACA && make clean && make`
- Run with LES mode (`id_visc = integer(8)` in mod_globals) to exercise the fixed Hsgs path.
- Run a short simulation (`nt = 10`) and confirm `G(5,...)` energy flux is finite and matches order-of-magnitude of `E(5,...)` and `F(5,...)`.
- For the hybrid-kernel UB fix: rebuild `3D_solver/SBLI` or `3D_solver/TBL` (Hybrid scheme cases) and confirm no change in outputs (it's a no-op fix).
