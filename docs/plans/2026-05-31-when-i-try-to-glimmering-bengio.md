# Plan: Fix TGV Case Deadlock

## Context

When running the TGV cases (NSTGV, ETGV) with `mpirun -n 2`, the solver deadlocks immediately during initialization. The root cause is a MPI message tag mismatch in `pre_calc` in `3D_solver/src/preprocess.f90.fypp`.

## Root Cause

### The architecture

The solver uses a compute/IO rank pairing:
- **Even rank (0)**: does GPU computation; sends field data to the odd rank for I/O
- **Odd rank (1)**: does no GPU computation; receives field data and writes VTK output

During initialization, the even rank calls `pre_calc`, which sends `ke0`, `entropy0`, `rho1d`, `p1d`, `v1d` to the odd rank. The odd rank then calls `send_recv_for_print_odd` (→ `send_recv_for_print_odd3`) to receive those fields.

### The tag mismatch

`pre_calc` (lines 97–101 of `3D_solver/src/preprocess.f90.fypp`) sends all five messages using `tag = myrank+1`:

```fortran
call MPI_SEND (ke0,      ..., myrank+1, myrank+1, ...)   ! tag = 1
call MPI_SEND (entropy0, ..., myrank+1, myrank+1, ...)   ! tag = 1
call MPI_ISEND(rho1d,    ..., myrank+1, myrank+1, ...)   ! tag = 1  ← wrong
call MPI_ISEND(p1d,      ..., myrank+1, myrank+1, ...)   ! tag = 1  ← wrong
call MPI_ISEND(v1d,      ..., myrank+1, myrank+1, ...)   ! tag = 1  ← wrong
```

`send_recv_for_print_odd3` (lines 319–327 of `src/print.f90`) expects **distinct** tags for the three field arrays:

```fortran
MPI_IRECV(rho1d, ..., myrank-1, 3*myrank-2, ...)   ! tag = 1
MPI_IRECV(p1d,   ..., myrank-1, 3*myrank-1, ...)   ! tag = 2  ← no match ever sent
MPI_IRECV(v1d,   ..., myrank-1, 3*myrank,   ...)   ! tag = 3  ← no match ever sent
```

**What happens at runtime (rank 0 → rank 1, nranks=2):**

| Message | Sent tag | Expected tag | Result |
|---------|----------|--------------|--------|
| ke0     | 1        | 1            | ✓ received |
| entropy0| 1        | 1            | ✓ received |
| rho1d   | 1        | 1            | ✓ received |
| p1d     | 1        | **2**        | ✗ rank 1 blocks forever |
| v1d     | 1        | **3**        | ✗ rank 1 blocks forever |

Rank 1 hangs on `MPI_WAITALL` waiting for tags 2 and 3. Rank 0 hangs on `MPI_WAITALL` for its ISENDs of p1d and v1d (tag=1) to complete — which requires rank 1 to receive them, but rank 1 is stuck on tag=2. Classic deadlock.

The time loop's `send_recv_for_print_even3` already uses the correct distinct tags (lines 279–285 of `src/print.f90`):
```fortran
MPI_ISEND(rho1d, ..., myrank+1, 3*(myrank+1)-2, ...)   ! tag = 1
MPI_ISEND(p1d,   ..., myrank+1, 3*(myrank+1)-1, ...)   ! tag = 2
MPI_ISEND(v1d,   ..., myrank+1, 3*(myrank+1),   ...)   ! tag = 3
```

`pre_calc` must use the same convention.

## Fix

**File:** `3D_solver/src/preprocess.f90.fypp`, lines 99–101

Change:
```fortran
call MPI_ISEND(rho1d, nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(1), ierr)
call MPI_ISEND(p1d,   nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(2), ierr)
call MPI_ISEND(v1d,   nx*ny*nz*3, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(3), ierr)
```

To:
```fortran
call MPI_ISEND(rho1d, nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
call MPI_ISEND(p1d,   nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
call MPI_ISEND(v1d,   nx*ny*nz*3, MPI_REAL4, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
```

This makes the initialization send tags match what `send_recv_for_print_odd3` expects — the same convention already used in `send_recv_for_print_even3`.

## Verification

```bash
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j
cd build && mpirun -n 2 ./a.out
```

Expected: no hang, VTK output appears in `build/data/`, `ke0` and `entropy0` log files written by rank 1.

Also run the ETGV case to verify:
```bash
cd 3D_solver/ETGV
cmake -B build && cmake --build build -j
cd build && mpirun -n 2 ./a.out
```

Run the regression test suite:
```bash
pytest ouxsbli/tests/test_etgv.py -v
```
