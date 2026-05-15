# Staged Changes Review

## Context

The staged changes refactor `calc_steps.f90`, `preprocess.f90`, `calc_time_dev.f90`, and `calc_para.f90` to precompute the per-cell metric products `dtdxdy`, `dtdydz`, `dtdzdx` on the CPU and pass them as 2D device arrays to the RK kernels — eliminating per-thread arithmetic in the hot CUDA loops. `calc_para.f90` also replaces implicit array assignment with explicit `cudaMemcpy`/`cudaMemcpyAsync` to pipeline the host↔device transfers around MPI calls. `mod_globals.f90` changes are a temporary test configuration (small grid, 1 step).

---

## Bugs Found

### Bug 1 — **Critical** — Double `deallocate` in `RungeKutta_3rd` (`calc_time_dev.f90:311-313`)

```fortran
! NEW (buggy)
deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, xix, etay, zetaz, Jacobian, dtdxdy, dtdydz, dtdzdx)
deallocate(dtdxdy, dtdydz, dtdzdx)   ! ← leftover line, double-free at runtime
```

`dtdxdy/dtdydz/dtdzdx` appear in both `deallocate` statements. The second `deallocate` is a cut-paste artifact and must be removed. The other three `RungeKutta_*` routines have a single correct `deallocate` statement.

---

### Bug 2 — **Critical** — Invalid CUDA stream handles in `calc_para.f90` (lines ~205-210)

```fortran
stat = cudaMemcpyAsync(Qr1d_right, Qr_right, ..., cudaMemcpyHostToDevice, 1)  ! stream=1
stat = cudaMemcpyAsync(Qr1d_left,  Qr_left,  ..., cudaMemcpyHostToDevice, 2)  ! stream=2
stat = cudaStreamSynchronize(1)
stat = cudaStreamSynchronize(2)
```

Integer literals `1` and `2` are used as `cudaStream_t` handles. Valid stream handles can only be obtained from `cudaStreamCreate`. In CUDA, the legacy default stream handle is `1` (`cudaStreamLegacy`), so both calls may silently redirect to the same stream, defeating the async intent. The fix is:

```fortran
integer(cuda_stream_kind) :: stream1, stream2
stat = cudaStreamCreate(stream1)
stat = cudaStreamCreate(stream2)
! ... use stream1, stream2 ...
stat = cudaStreamSynchronize(stream1)
stat = cudaStreamSynchronize(stream2)
stat = cudaStreamDestroy(stream1)
stat = cudaStreamDestroy(stream2)
```

These should be declared and used within `exchange_cyclic`.

---

### Bug 3 — **Wrong scheme selection** — `mod_globals.f90` (`id_scheme` kind change)

```fortran
- real(2), parameter :: id_scheme = 0   ! → KEEP scheme
+ integer(2), parameter :: id_scheme = 0
```

Per CLAUDE.md, only `real(2)` and `real(4)`/`real(8)` are defined scheme selectors. `integer(2)` is not a valid scheme kind and will dispatch incorrectly in `calc_flux_base.f90`. This is likely a typo — should remain `real(2)`.

---

### Bug 4 — **Wrong BCs for NSTGV** — `mod_globals.f90`

```fortran
- logical, parameter :: id_bc_x = .false.   ! periodic
+ logical, parameter :: id_bc_x = .true.    ! non-periodic wall BC
```

TGV uses periodic BCs in all directions. `.true.` activates wall/non-periodic BCs which will corrupt the solution. This seems an accidental change bundled with the small-grid test config (`nx=66`, `nt=1`). Should revert to `.false.`.

---

### Minor — Unused `use cooperative_groups` in `calc_steps.f90` (line 3)

```fortran
+ use cooperative_groups
```

This module is imported but never used in any visible kernel. May cause a compile error or warning depending on HPC SDK version. Remove unless cooperative-group functionality is planned.

---

## Optimization Observations (already in the staged changes — good)

1. **Precomputed 2D metric arrays** (`dtdxdy`, `dtdydz`, `dtdzdx`): eliminates 6 multiplications and 3 additions per thread per kernel call from the hot inner loop. Reduction in arithmetic plus better memory access pattern (2D lookup vs 1D gather-and-multiply).

2. **Local scalar caching before the `do l` loop** in all kernels:
   ```fortran
   dtdxdy_tmp = dtdxdy(i,j)   ! load once, reuse 5×
   dtdydz_tmp = dtdydz(j,k)
   dtdzdx_tmp = dtdzdx(i,k)
   ```
   Avoids repeated global memory reads across the 5-variable loop.

3. **`calc_step4` register temp instead of redundant write**: `R = Rs(i,l,j,k) + R` then `Rs(i,l,j,k) = 0.d0` avoids the intermediate global-memory write-then-read that the old code did.

4. **Async H2D upload in `calc_para.f90`**: the intent to overlap two H2D copies on separate streams is sound — the fix is just to create real stream handles (Bug 2).

---

## Additional Optimization Ideas

### A. Overlap D2H and H2D in `calc_para.f90`

Currently `cudaMemcpy` (synchronous, stream 0) is used for D2H before each MPI call. If the `flatten` kernel has completed, both D2H copies could also be made async:

```
Stream A: flatten → cudaMemcpyAsync D2H Qs_left → (MPI blocks host, not device)
```

MPI blocks the CPU thread, not the GPU, so a cudaMemcpyAsync D2H that finishes before `MPI_SENDRECV` returns is fine, but there's no benefit unless GPU work can be overlapped simultaneously. Given that `flatten` must finish before any D2H anyway, synchronous D2H is correct here.

### B. Persistent streams for `calc_para.f90`

Creating/destroying streams per call to `exchange_cyclic` adds overhead. If this subroutine is called every time step, the two streams should be created once (at `pre_calc` or init time) and reused across time steps.

### C. `dtdzdx` naming vs. semantics

Note that `dtdzdx(i,k) = 0.25 * dt * dz(k)*dx(i)` — this is the area used with the **y-flux** `F`, which is the y-face area. The name `dtdzdx` literally reads "dt × dz × dx", which is the z-x plane area used for the y-flux divergence. Consistent with the original variable names but worth confirming against the flux sign convention.

---

## Files to Fix

| File | Change needed |
|------|--------------|
| [3D_solver/src/calc_time_dev.f90](3D_solver/src/calc_time_dev.f90) | Remove duplicate `deallocate(dtdxdy, dtdydz, dtdzdx)` (~line 313) |
| [3D_solver/src/calc_para.f90](3D_solver/src/calc_para.f90) | Declare `integer(cuda_stream_kind)` stream handles; call `cudaStreamCreate`/`Destroy` |
| [3D_solver/NSTGV/mod_globals.f90](3D_solver/NSTGV/mod_globals.f90) | Revert `id_scheme` to `real(2)`; revert `id_bc_x/y/z` to `.false.` |
| [3D_solver/src/calc_steps.f90](3D_solver/src/calc_steps.f90) | Remove unused `use cooperative_groups` |
