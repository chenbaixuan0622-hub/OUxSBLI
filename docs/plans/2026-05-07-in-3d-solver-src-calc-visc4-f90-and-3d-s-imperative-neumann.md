# Plan: Pipeline txx/txy/txz in calc_visc4 using pipelineMemcpyAsync

## Context

Currently `calc_visc4.f90` and `calc_visc4_internal.f90` load **7 shared memory arrays** upfront (u, v, w, uy/vy/uz/wz or their y/z analogues) before computing txx, txy, txz in sequence. The three stress components are independent:

| Component | Velocity smem | Gradient smem |
|-----------|--------------|---------------|
| txx (x-dir) | u | vy, wz |
| txy (x-dir) | v | uy |
| txz (x-dir) | w | uz |

By staging the loads — one velocity component + its gradients per stage — shared memory can be reduced from **7 arrays → 3 arrays**, boosting occupancy. `mu3` and boundary scalars live in registers and persist across stages.

---

## Critical Files

- [3D_solver/src/load_smem_visc4.f90](3D_solver/src/load_smem_visc4.f90) — shared-memory loading helpers (modify)
- [3D_solver/src/calc_visc4.f90](3D_solver/src/calc_visc4.f90) — 6 global kernels (modify)
- [3D_solver/src/calc_visc4_internal.f90](3D_solver/src/calc_visc4_internal.f90) — 3 interior-only kernels (modify)
- [3D_solver/src/calc_visc_me4_base.f90](3D_solver/src/calc_visc_me4_base.f90) — `calc_tau_straight`, `calc_tau_cross` (read-only)

---

## Implementation

### 1. `load_smem_visc4.f90` — replace 3 subroutines with 9 stage subroutines

Replace `load_smem_visc4_x/y/z` with three stage subroutines per direction. Update `public` list accordingly.

**x-direction signatures** (y/z analogous):

```fortran
! Stage 1: async load u → vel; sync compute vy → g1, wz → g2
attributes(device) subroutine load_smem_visc4_x_s1(it, jt, kt, j, k, nx, ny, nz, inv_dy, inv_dz, Q, vel, g1, g2)

! Stage 2: async load v → vel; sync compute uy → g1
attributes(device) subroutine load_smem_visc4_x_s2(it, jt, kt, j, k, nx, ny, nz, inv_dy, Q, vel, g1)

! Stage 3: async load w → vel; sync compute uz → g1
attributes(device) subroutine load_smem_visc4_x_s3(it, jt, kt, j, k, nx, ny, nz, inv_dz, Q, vel, g1)
```

**Stage mapping per direction:**

| Direction | Stage 1 (vel, g1, g2) | Stage 2 (vel, g1) | Stage 3 (vel, g1) |
|-----------|----------------------|--------------------|-------------------|
| x (`_Ev`) | u, vy=dv/dy, wz=dw/dz | v, uy=du/dy | w, uz=du/dz |
| y (`_Fv`) | v, wz=dw/dz, ux=du/dx | u, vx=dv/dx | w, vz=dv/dz |
| z (`_Gv`) | w, ux=du/dx, vy=dv/dy | u, wx=dw/dx | v, wy=dw/dy |

**Non-TMA body pattern (identical for all 9 subroutines, only component indices differ):**
```fortran
! Stage n: load Q(:, comp_vel, ...) → vel; compute g1 [and g2 for s1]
do ii = ...
  if (in_bounds) call pipelineMemcpyAsync(vel(idx), Q(i, comp_vel, j, k))
enddo
call pipelineCommit()
do ii = ...
  ! compute g1 = gradient of comp_g1 in perp direction (from Q, no smem)
  ! compute g2 = (s1 only) gradient of comp_g2 in second perp direction
enddo
call pipelineWaitPrior(0)
call syncthreads()
```

**TMA body pattern** (`#if _CUDA_ARCH_ >= 900`): use one `integer(8), shared :: barrier` per subroutine; re-init barrier at entry with `barrier_init`, issue single `tma_bulk_load` for the one velocity component, compute gradients synchronously, then `syncthreads + barrier_arrive + barrier_try_wait_sleep`.

---

### 2. `calc_visc4.f90` — reduce smem + 3-stage structure in each kernel

**Shared memory change (all 6 kernels):**
```fortran
! BEFORE (7 arrays):
real(8), dimension(...), shared :: u, v, w, uy, vy, uz, wz  ! or ux,vx,vz,wz etc.

! AFTER (3 arrays):
real(8), dimension(...), shared :: smem_a  ! vel: u→v→w across stages
real(8), dimension(...), shared :: smem_b  ! g1:  vy→uy→uz across stages
real(8), dimension(...), shared :: smem_c  ! g2:  wz (stage 1 only)
```

**Variable scope change:** Move `mu3(3)` (and `mut3(3)` for LES) and the boundary scalars `mx, muy, mvy, muz, mwz` from inside nested `block` constructs to the **subroutine declaration** area so they persist across the 3 pipeline stages.

**Kernel body structure (shown for `calc_Ev4`; same pattern for Fv4, Gv4 and LES variants):**

```fortran
! common setup
it=threadIdx%x; jt=threadIdx%y; kt=threadIdx%z
j = (blockIdx%y-1)*blockDim%y + jt + 1
k = (blockIdx%z-1)*blockDim%z + kt + 1
offset_yz = (jt-1)*sx + (kt-1)*sx*sy

! ── STAGE 1 ──────────────────────────────────────────────
call load_smem_visc4_x_s1(it,jt,kt, j,k, nx,ny,nz, dy,dz, Q, smem_a,smem_b,smem_c)
i = (blockIdx%x-1)*blockDim%x + it
idx = it + offset_yz
if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
if (interior) then
    mu3(:) = ...          ! global mu, no smem
    kTx3(:) = ...; kTx = flux4(kTx3)    ! global T + register mu3
    call calc_tau_straight(mu3, smem_a(idx-2:idx+3), smem_b(...), smem_c(...), dx(i), txx, utxx)
else
    mx = 0.5*(mu(i,j,k)+mu(i+1,j,k))
    kTx = Cp_over_Pr*mx*(-T(i,j,k)+T(i+1,j,k))*dx(i)
    ! compute muy, mvy, muz, mwz from Q (unchanged from current code)
    mux = mx*(-Q(i,2,j,k)+Q(i+1,2,j,k))*dx(i)   ! Q replaces smem_a(idx)
    txx = two_third*(2.*mux - mvy - mwz)
    utxx = 0.5*(Q(i,2,j,k)+Q(i+1,2,j,k))*txx
endif
call syncthreads()   ! guard smem_a/b/c before stage 2 overwrites

! ── STAGE 2 ──────────────────────────────────────────────
call load_smem_visc4_x_s2(it,jt,kt, j,k, nx,ny,nz, dy, Q, smem_a,smem_b)
if (interior) then
    call calc_tau_cross(mu3, smem_a(idx-2:idx+3), smem_b(...), dx(i), txy, vtxy)
else
    mvx = mx*(-Q(i,3,j,k)+Q(i+1,3,j,k))*dx(i)
    txy = muy + mvx
    vtxy = 0.5*(Q(i,3,j,k)+Q(i+1,3,j,k))*txy
endif
call syncthreads()   ! guard before stage 3

! ── STAGE 3 ──────────────────────────────────────────────
call load_smem_visc4_x_s3(it,jt,kt, j,k, nx,ny,nz, dz, Q, smem_a,smem_b)
if (interior) then
    call calc_tau_cross(mu3, smem_a(idx-2:idx+3), smem_b(...), dx(i), txz, wtxz)
else
    mwx = mx*(-Q(i,4,j,k)+Q(i+1,4,j,k))*dx(i)
    txz = mwx + muz
    wtxz = 0.5*(Q(i,4,j,k)+Q(i+1,4,j,k))*txz
endif

! ── STORE ────────────────────────────────────────────────
E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
```

For **LES variants**: additionally lift `mut3(3)` to subroutine scope (stage 1 computes it alongside `mu3`); `Hsgs` stays register, computed in stage 3 (only uses global T, u/v/w registers from smem_a).

For **`calc_Fv4`/`calc_Fv_LES4`**: same 3-stage pattern; stages load v/u/w and their gradients as in the table above.

For **`calc_Gv4`/`calc_Gv_LES4`**: same; stages load w/u/v.

---

### 3. `calc_visc4_internal.f90` — same smem reduction, simpler (no boundary path)

Each of `calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`:
- Replace 7 smem arrays with 3 (`smem_a`, `smem_b`, `smem_c`)
- Lift `mu3(3)` to subroutine scope
- Same 3-stage structure; no else-branch needed
- `if (...interior check...)` guard wraps the tau calls inside each stage
- Two `syncthreads()` calls between stages (after stage 1 tau, after stage 2 tau)
- Final store moved outside the inner `block`

---

## Key Invariants

- `pipelineWaitPrior(0)` inside each `load_smem_visc4_*_sN` subroutine ensures the async load committed in that subroutine is complete before returning.
- The `syncthreads()` at the **end of each load subroutine** ensures smem is visible to all threads before they read it.
- The `syncthreads()` called by the **global kernel between stages** ensures all threads finish reading stage N's smem before any thread begins overwriting it in stage N+1.
- `mu3`, `kTx`, and boundary scalars (`mx`, `muy`, …) are registers — they survive across `syncthreads()` boundaries.
- Boundary threads read velocity face values directly from `Q(i, comp, j, k)` instead of smem; this is equivalent since smem was a direct copy from Q.

---

## Verification

```bash
cd 3D_solver/NSTGV   # exercises NS viscous path (id_visc = integer(2))
make clean && make
bash calc.sh
# compare VTK output against reference run
bash profile.sh      # verify occupancy increase from smem reduction
```

Also verify in a case using `calc_visc4_internal` (e.g. KHI with periodic BCs).
