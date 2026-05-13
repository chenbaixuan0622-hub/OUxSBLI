# Plan: Move boundary/ghost-cell metric computation from `set_metrics_curv` into `set.f90`

## Context

`set_metrics_curv` in `src/set_coordinate.f90` currently ends with three extra sections (lines 425–480) beyond the interior cell computation:

1. **eta-boundary cells** (lines 425–441): 1-sided stencils at j=1 and j=ny
2. **xi-boundary cells** (lines 442–470): 1-sided stencils at i=1 and i=nx
3. **Ghost cell enforcement** (lines 471–480): copies xi-ghost metrics to periodic counterparts
   (`Jac(1,n)=Jac(nx-1,n)`, `Jac(nx,n)=Jac(2,n)`, etc.)

Section 3 contains an embedded NACA comment ("Q(1)=Q(nx-1), Q(nx)=Q(2) … near the TE") and
explicitly assumes a periodic xi direction (O-grid). This is **correct for NACA but wrong for CORN**,
which has non-periodic inlet/outlet xi boundaries. Currently CORN's `set_metrics` wrapper calls
`set_metrics_curv` and silently receives incorrect periodic ghost metric values.

The user's assessment is correct: all three sections are case-specific and should live in each
case's `set.f90`.

## Files to Modify

| File | Action |
|------|--------|
| `src/set_coordinate.f90` | Remove lines 425–480 from `set_metrics_curv` (keep only face normals + interior cells) |
| `3D_solver_curv/CORN/set.f90` | Add eta-boundary + xi-boundary metric blocks inside `set_metrics`; omit ghost-cell periodic enforcement |

`3D_solver_curv/NACA/set.f90` is incomplete (untracked, no `set_metrics` yet); it is left for when
NACA is finished. A note will be added in the truncated `set_metrics_curv` docstring.

## Implementation Steps

### Step 1 — Truncate `set_metrics_curv` in `src/set_coordinate.f90`

Remove lines 425–480 so the subroutine ends after the interior-cell loop (line 424).
Update the leading comment to reflect that boundary/ghost cells are the caller's responsibility.

Result: `set_metrics_curv` computes only
- xi-face normals (`n_xi_x`, `n_xi_y`)
- eta-face normals (`n_eta_x`, `n_eta_y`)
- Interior cell-centre metrics (`xi_x`, `xi_y`, `eta_x`, `eta_y`, `Jac`) for `m=2..nx-1, n=2..ny-1`

All boundary/corner cells are left uninitialised for the caller to fill.

### Step 2 — Extend `set_metrics` in `3D_solver_curv/CORN/set.f90`

After `call set_metrics_curv(...)`, add (inline, no new subroutine):

```fortran
! --- eta-boundary cells: j=1 (lower wall, 1-sided forward) and j=ny (upper wall, backward) ---
do m = 2, nx-1
  xxi  = 0.5d0*(x_phys(m+1,1)-x_phys(m-1,1))
  yxi  = 0.5d0*(y_phys(m+1,1)-y_phys(m-1,1))
  xeta = x_phys(m,2)-x_phys(m,1);  yeta = y_phys(m,2)-y_phys(m,1)
  J2          = xxi*yeta - xeta*yxi
  Jac_cpu(m,1) = J2
  xi_x_cpu(m,1) =  yeta/J2;  xi_y_cpu(m,1) = -xeta/J2
  eta_x_cpu(m,1)= -yxi /J2; eta_y_cpu(m,1) =  xxi /J2
  xxi  = 0.5d0*(x_phys(m+1,ny)-x_phys(m-1,ny))
  yxi  = 0.5d0*(y_phys(m+1,ny)-y_phys(m-1,ny))
  xeta = x_phys(m,ny)-x_phys(m,ny-1);  yeta = y_phys(m,ny)-y_phys(m,ny-1)
  J2           = xxi*yeta - xeta*yxi
  Jac_cpu(m,ny) = J2
  xi_x_cpu(m,ny) =  yeta/J2;  xi_y_cpu(m,ny) = -xeta/J2
  eta_x_cpu(m,ny)= -yxi /J2; eta_y_cpu(m,ny) =  xxi /J2
enddo
! --- xi-boundary cells: i=1 (inlet, 1-sided fwd) and i=nx (outlet, backward), all j ---
do n = 1, ny
  xxi = x_phys(2,n)-x_phys(1,n);  yxi = y_phys(2,n)-y_phys(1,n)
  if (n == 1) then
    xeta = x_phys(1,2)-x_phys(1,1);  yeta = y_phys(1,2)-y_phys(1,1)
  elseif (n == ny) then
    xeta = x_phys(1,ny)-x_phys(1,ny-1);  yeta = y_phys(1,ny)-y_phys(1,ny-1)
  else
    xeta = 0.5d0*(x_phys(1,n+1)-x_phys(1,n-1))
    yeta = 0.5d0*(y_phys(1,n+1)-y_phys(1,n-1))
  endif
  J2          = xxi*yeta - xeta*yxi
  Jac_cpu(1,n) = J2
  xi_x_cpu(1,n) =  yeta/J2;  xi_y_cpu(1,n) = -xeta/J2
  eta_x_cpu(1,n)= -yxi /J2; eta_y_cpu(1,n) =  xxi /J2
  xxi = x_phys(nx,n)-x_phys(nx-1,n);  yxi = y_phys(nx,n)-y_phys(nx-1,n)
  if (n == 1) then
    xeta = x_phys(nx,2)-x_phys(nx,1);  yeta = y_phys(nx,2)-y_phys(nx,1)
  elseif (n == ny) then
    xeta = x_phys(nx,ny)-x_phys(nx,ny-1);  yeta = y_phys(nx,ny)-y_phys(nx,ny-1)
  else
    xeta = 0.5d0*(x_phys(nx,n+1)-x_phys(nx,n-1))
    yeta = 0.5d0*(y_phys(nx,n+1)-y_phys(nx,n-1))
  endif
  J2           = xxi*yeta - xeta*yxi
  Jac_cpu(nx,n) = J2
  xi_x_cpu(nx,n) =  yeta/J2;  xi_y_cpu(nx,n) = -xeta/J2
  eta_x_cpu(nx,n)= -yxi /J2; eta_y_cpu(nx,n) =  xxi /J2
enddo
! Note: no periodic ghost-cell enforcement — CORN has Dirichlet inlet (i=1) and
! zero-gradient outlet (i=nx) in set_bc; those flow-state copies dominate anyway.
```

`set_metrics` already receives `x_phys_g`, `y_phys_g` as `intent(in)` and the output arrays
as allocatables, so no signature change is needed.

Declare `real(8) :: xxi, yxi, xeta, yeta, J2` as locals inside `set_metrics`.

## Verification

1. Confirm `set_metrics_curv` ends at the closing `enddo` of the interior loop (after `Jac(m,n) = J2`).
2. Confirm CORN `set_metrics` fills every cell of `Jac_cpu`, `xi_x_cpu`, `xi_y_cpu`, `eta_x_cpu`,
   `eta_y_cpu` — previously corner cells (e.g., `Jac_cpu(1,1)`) were left as whatever 1-sided
   stencil computed then overwritten by periodic copies; after this change they hold correct 1-sided
   values.
3. Build CORN: `cd 3D_solver_curv/CORN && make clean && make` — must compile without errors.
4. Short smoke-run: `bash calc.sh` for a few steps; check that printed residuals or VTK output are
   physically consistent (no NaN, pressure positive).
