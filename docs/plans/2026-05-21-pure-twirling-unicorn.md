# Fix MPI deadlock in STZ case (pre_calc_zdec)

## Context

The STZ case (z-direction Sod shock tube) uses 4 MPI ranks: even ranks (0, 2) compute on GPU, odd ranks (1, 3) do file I/O. After fixing the even/odd pattern in the previous session, a deadlock remains during the t=0 snapshot step.

## Root Cause

**File:** `3D_solver/src/preprocess.f90`, `pre_calc_zdec` subroutine, lines 207–210.

Execution sequence:

| Compute rank 0 | IO rank 1 |
|---|---|
| `print_vtk(0,...,myrank+1=1,...)` — **writes file directly, no MPI** | — |
| `MPI_SEND(ke0, tag=1, dest=1)` | `MPI_RECV(ke0, tag=1, src=0)` ← receives |
| `MPI_SEND(entropy0, tag=1, dest=1)` | `MPI_RECV(entropy0, tag=1, src=0)` ← receives |
| `MPI_BARRIER` ← **stuck here** | `send_recv_for_print_odd(..., step=0, ...)` posts `MPI_IRECV(rho1d/p1d/v1d)` → `MPI_WAITALL` ← **stuck here forever** |

The IO rank calls `send_recv_for_print_odd` (from `print.f90:311`) expecting the compute rank to send `rho1d`, `p1d`, `v1d` for the step=0 snapshot. But the compute rank already wrote the file directly via `print_vtk` and never sends those arrays. Neither rank reaches the barrier → deadlock.

## Fix

**One file changes:** `3D_solver/src/preprocess.f90`

In `pre_calc_zdec`, replace the direct `print_vtk` call with MPI sends of `rho1d`/`p1d`/`v1d` to the IO rank, so `send_recv_for_print_odd` on the IO side can receive them and write the file. The `ke0`/`entropy0` sends must stay **before** the array sends to match the IO rank's receive order.

### Change at lines 167 and 207–210

**Add declarations** (line 167, after existing `integer i, j, k, l, ierr`):
```fortran
integer ireq3(3), istat3(MPI_STATUS_SIZE, 3)
```

**Replace** (lines 207–210):
```fortran
! OLD — writes file directly, never sends rho1d/p1d/v1d:
call print_vtk(0, nx, ny, nz, myrank+1, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
```

**With** (ke0/entropy0 sent first so IO rank can unblock, then rho1d/p1d/v1d):
```fortran
! NEW — IO rank receives ke0/entropy0 first, then rho1d/p1d/v1d via send_recv_for_print_odd:
call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
call MPI_ISEND(rho1d, nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(1), ierr)
call MPI_ISEND(p1d,   nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(2), ierr)
call MPI_ISEND(v1d,   nx*ny*nz*3, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(3), ierr)
call MPI_WAITALL(3, ireq3, istat3, ierr)
```

No changes needed to IO rank code (`else` branches in `calc_time_dev.f90`) — that logic is already correct.

## Why tags are safe

- `ke0`/`entropy0` use tag `myrank+1` (e.g. tag=1 from rank 0 → rank 1). IO rank receives with tag `myrank=1`. ✓
- `rho1d`/`p1d`/`v1d` also use tag `myrank+1=1`. IO rank's `send_recv_for_print_odd3` receives with tag `myrank=1`. ✓
- No collision: `MPI_SEND(ke0)` is blocking and completes before `MPI_ISEND(rho1d)` starts; MPI preserves per-tag ordering between the same rank pair.

## Verification

```bash
cd 3D_solver/STZ
make clean && make
bash calc.sh
```

Expected: all 4 ranks complete without hang; `data/1/` and `data/3/` contain `Q00000.vtr` through `Q00020.vtr`. Open in ParaView (combine z-slabs) — final snapshot should show the classic Sod profile at T=0.2.
