# Plan: Fix Ghost-Cell Jacobian in Periodic BC (`set_bc_cyclic_x`)

## Context

The WING O-grid solver uses Jacobian-scaled conservative variables: `QJ = q / J_2D`.
The periodic xi-direction BC (`set_bc_cyclic_x`) currently does a simple copy:

```fortran
Q(1,:,:,:)  = Q(nx-1,:,:,:)   ! left ghost ← last interior
Q(nx,:,:,:) = Q(2,:,:,:)      ! right ghost ← first interior
```

This is only correct if `J(1,j) = J(nx-1,j)` and `J(nx,j) = J(2,j)`.

**However, `set_metrics_curv` (`src/set_coordinate.f90` lines 524–551) uses a
one-sided stencil for the xi-boundary ghost cells:**

```fortran
! i=1  (forward difference):
xxi = x_phys(2,n) - x_phys(1,n)
! i=nx (backward difference):
xxi = x_phys(nx,n) - x_phys(nx-1,n)
```

With the periodic convention `x(1)=x(nx-1)` and `x(nx)=x(2)`, near the
trailing edge `x(2) ≈ x(nx-1) ≈ chord`, so `xxi ≈ 0` at the ghost cells.
This makes `J(1,j)` and `J(nx,j)` nearly degenerate — and definitely **not
equal** to `J(nx-1,j)` and `J(2,j)` (which use central differences).

As a result the simple copy `Q(1)=Q(nx-1)` is wrong: it produces
`QJ(1) = q(1)/J(1) ≠ q(nx-1)/J(nx-1) = QJ(nx-1)`, corrupting the periodic
seam every time-step.

---

## Fix

Add a block at the **end of `set_metrics_curv`** in `src/set_coordinate.f90`
to overwrite ghost-cell metrics with their interior periodic counterparts:

```fortran
! Enforce xi-periodic consistency: ghost cell metrics = interior periodic pair
do n = 1, ny
  Jac(1,n)   = Jac(nx-1,n);  Jac(nx,n)  = Jac(2,n)
  xi_x(1,n)  = xi_x(nx-1,n); xi_x(nx,n) = xi_x(2,n)
  xi_y(1,n)  = xi_y(nx-1,n); xi_y(nx,n) = xi_y(2,n)
  eta_x(1,n) = eta_x(nx-1,n); eta_x(nx,n)= eta_x(2,n)
  eta_y(1,n) = eta_y(nx-1,n); eta_y(nx,n)= eta_y(2,n)
enddo
```

**Why this is the right location:**
- `set_metrics_curv` is called once in `set_metrics` (`3D_solver/WING/set.f90`)
  before the time loop — the corrected arrays are then copied to the GPU once.
- Ghost cells represent the same physical location as their interior counterparts,
  so their metrics should be identical. The one-sided stencil was only a
  numerical artefact of not knowing the periodic neighbour.
- This also fixes `eta_x(1,j=1)` and `eta_x(nx,j=1)` used in the wall BC
  loop (`set_bc`) for the corner cells i=1 and i=nx.
- With `J(1)=J(nx-1)` and `J(nx)=J(2)` enforced, the simple QJ copy in
  `set_bc_cyclic_x` becomes exactly correct with no signature change required.

**`set_metrics_curv` is only called for the WING/curvilinear case** —
Cartesian cases use `set_Jacobian_xy2/3` instead, so this change is safe.

---

## Files to Modify

| File | Change |
|------|--------|
| `src/set_coordinate.f90` | Add periodic ghost-cell metric copy at end of `set_metrics_curv` (after line 551) |

No other files need modification.

---

## Verification

```bash
cd 3D_solver/WING
make clean && make
bash calc.sh
```

Check in ParaView that the wake seam shows no discontinuity in density or
pressure, and that the periodic BC no longer generates spurious artifacts at
i=1 / i=nx.
