# Fix: High-order `set_bc` in `3D_solver/STZ/set.f90`

## Context

`overlap_fb = kind(id_accuracy)/3 + 1` controls how many ghost cells each boundary needs:
- 2nd order (kind=2): `overlap_fb = 1` — current code is correct
- 4th order (kind=4): `overlap_fb = 2` — needs 2 ghost cells; current code only fills 1
- 6th order (kind=8): `overlap_fb = 3` — needs 3 ghost cells

In STZ, the z-boundary already loops over `ovlp` (correct). The x and y slip walls hardcode a single ghost cell and need the same treatment.

---

## Ghost-cell formulas

Throughout, `ovlp = kind(id_accuracy)/3 + 1`.

**Slip wall**, normal-velocity component `nc` negated (wall at low end of index `d`):
```fortran
do m = 0, ovlp-1
  QJ(..., l, d=ovlp-m,    ...) = sign(l,nc) * QJ(..., l, d=ovlp+1+m, ...)
  ! sign = -1 for l==nc; +1 for l!=nc
enddo
! High end: ghost d=nd-ovlp+1+m  ←  interior d=nd-ovlp-m (same sign pattern)
```

**No-slip wall** (all velocity antisymmetric, scalars symmetric):
```fortran
do m = 0, ovlp-1
  QJ(i,1,ovlp-m,k)  =  QJ(i,1,ovlp+1+m,k)   ! rho
  QJ(i,2,ovlp-m,k)  = -QJ(i,2,ovlp+1+m,k)   ! rho*u
  QJ(i,3,ovlp-m,k)  = -QJ(i,3,ovlp+1+m,k)   ! rho*v
  QJ(i,4,ovlp-m,k)  = -QJ(i,4,ovlp+1+m,k)   ! rho*w
  QJ(i,5,ovlp-m,k)  =  QJ(i,5,ovlp+1+m,k)   ! rho*E
enddo
```
This replaces the current `rho*u=0` approach with the antisymmetric image approach, which correctly gives zero velocity at the wall face for all orders (face value = average = 0).

**Zero-gradient / Neumann** (outflow, symmetric image):
```fortran
do m = 0, ovlp-1
  QJ(..., ovlp-m, ...) = QJ(..., ovlp+1+m, ...)
enddo
```

**Riemann far-field** (multiple ghost cells): read from first interior cell `ny-ovlp`, compute Riemann state once, then apply to all `ovlp` ghost cells `j=ny-ovlp+1..ny` using their respective Jacobians.

**z-periodic** (generalised from 6th-order hardcode):
```fortran
do m = 1, ovlp
  QJ(i,l,j,m)            = QJ(i,l,j,nz-2*ovlp+m)
  QJ(i,l,j,nz-ovlp+m)    = QJ(i,l,j,ovlp+m)
enddo
```
Verify for ovlp=3: m=1→k=1←k=nz-5, m=2→k=2←k=nz-4, m=3→k=3←k=nz-3 (matches current hardcode).

---

## Change — `3D_solver/STZ/set.f90`

### x slip wall (lines ~78–92): replace with ovlp-loop
```fortran
!$cuf kernel do(3)<<<*,*>>>
do k = 1, nz
  do j = 1, ny
    do m = 0, ovlp-1
      QJ(ovlp-m,      1,j,k) =  QJ(ovlp+1+m,   1,j,k)
      QJ(ovlp-m,      2,j,k) = -QJ(ovlp+1+m,   2,j,k)
      QJ(ovlp-m,      3,j,k) =  QJ(ovlp+1+m,   3,j,k)
      QJ(ovlp-m,      4,j,k) =  QJ(ovlp+1+m,   4,j,k)
      QJ(ovlp-m,      5,j,k) =  QJ(ovlp+1+m,   5,j,k)
      QJ(nx-ovlp+1+m, 1,j,k) =  QJ(nx-ovlp-m,  1,j,k)
      QJ(nx-ovlp+1+m, 2,j,k) = -QJ(nx-ovlp-m,  2,j,k)
      QJ(nx-ovlp+1+m, 3,j,k) =  QJ(nx-ovlp-m,  3,j,k)
      QJ(nx-ovlp+1+m, 4,j,k) =  QJ(nx-ovlp-m,  4,j,k)
      QJ(nx-ovlp+1+m, 5,j,k) =  QJ(nx-ovlp-m,  5,j,k)
    enddo
  enddo
enddo
```

### y slip wall (lines ~94–108): replace with ovlp-loop
```fortran
!$cuf kernel do(3)<<<*,*>>>
do k = 1, nz
  do i = 1, nx
    do m = 0, ovlp-1
      QJ(i,1,ovlp-m,     k) =  QJ(i,1,ovlp+1+m,   k)
      QJ(i,2,ovlp-m,     k) =  QJ(i,2,ovlp+1+m,   k)
      QJ(i,3,ovlp-m,     k) = -QJ(i,3,ovlp+1+m,   k)
      QJ(i,4,ovlp-m,     k) =  QJ(i,4,ovlp+1+m,   k)
      QJ(i,5,ovlp-m,     k) =  QJ(i,5,ovlp+1+m,   k)
      QJ(i,1,ny-ovlp+1+m,k) =  QJ(i,1,ny-ovlp-m,  k)
      QJ(i,2,ny-ovlp+1+m,k) =  QJ(i,2,ny-ovlp-m,  k)
      QJ(i,3,ny-ovlp+1+m,k) = -QJ(i,3,ny-ovlp-m,  k)
      QJ(i,4,ny-ovlp+1+m,k) =  QJ(i,4,ny-ovlp-m,  k)
      QJ(i,5,ny-ovlp+1+m,k) =  QJ(i,5,ny-ovlp-m,  k)
    enddo
  enddo
enddo
```

No change to the z-boundary section (already uses `ovlp` loop).

---

## File Modified

`3D_solver/STZ/set.f90` — x and y slip wall kernels: replace single-cell assignments with ovlp-loop.

---

## Verification

```bash
cd 3D_solver/STZ && make clean && make && bash calc.sh
pytest ouxsbli/tests/test_st.py
```
