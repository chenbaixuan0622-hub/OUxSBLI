# Performance Analysis: `calc_step2_3` in `calc_steps.f90`

## Context

Profiling with NCU on the 513³ NSTGV grid (stable branch) identifies `calc_steps_calc_step2_3_` as the top bottleneck at **12.41 ms per call** (2 calls per TVD RK3 timestep). The pipeline branch's memcpy optimization has not touched this kernel — it shows identical metrics in both branches.

---

## Root Cause Diagnosis

| Metric | Value | Meaning |
|--------|-------|---------|
| DRAM Throughput | **80%** | Nearly saturating memory bandwidth |
| SM (Compute) Throughput | **13%** | Severely memory-bound; GPU is mostly idle |
| L2 Hit Rate | **33%** | 67% of L2 misses go straight to DRAM |
| Scheduler "No Eligible" | **87%** | Warps stall waiting for DRAM loads |
| Achieved Occupancy | **57%** | Register pressure limits blocks/SM to 10 |
| Registers/Thread | **44** | With 128 threads/block → block limit = 10 blocks/SM |

**Primary cause: AoS flux array layout with stride-5 warp access**

Flux arrays are declared as `E(5, nx-1, ny-2, nz-2)` (component-first). In Fortran column-major storage, consecutive flux cells `E(:, i, j, k)` and `E(:, i+1, j, k)` are 5 doubles (40 bytes) apart. Since all 32 threads in a warp vary in `i`, each component load is a **stride-5 access** across the warp.

Per warp per component-face read:
- **AoS (current):** threads span 32×40 = 1,280 bytes → ~10 L2 sectors fetched, 2 useful → **5× over-fetch**
- **SoA:** threads span 32×8 = 256 bytes → 2 L2 sectors fetched, 2 useful → **no over-fetch**

`calc_step2_3` reads 30 component-face values from E/F/G per grid point. With 5× over-fetch, total DRAM traffic is ~3.5× more than necessary. The 512³-scale arrays (5.4 GB each) completely exceed L2 capacity, so every miss goes to DRAM — there is no caching benefit to hide the over-fetch.

**Secondary cause: Register pressure limits warp parallelism**

44 registers/thread × 128 threads/block = 5,632 registers/block. With 65,536 registers/SM, only 11 blocks fit → Block Limit Registers = 10. The resulting 57% occupancy (40 warps/SM) is not enough to hide 400–600 cycle DRAM latency, causing 87% scheduler stalls.

---

## Suggestions (ranked by estimated impact)

### Suggestion 1 — SoA flux array layout: `E(5,N...)` → `E(N...,5)` [Highest impact, ~3-4× speedup on this kernel]

Change all three flux array declarations from component-first to component-last:

```fortran
! Before (AoS — component fastest, stride-5 across warp)
E(5, nx-1, ny-2, nz-2)
F(5, nx-2, ny-1, nz-2)
G(5, nx-2, ny-2, nz-1)

! After (SoA — spatial index fastest, stride-1 across warp)
E(nx-1, ny-2, nz-2, 5)
F(nx-2, ny-1, nz-2, 5)
G(nx-2, ny-2, nz-1, 5)
```

Access in `calc_R` changes from `E(comp, i, j, k)` → `E(i, j, k, comp)`. All producers of E/F/G also write in the new order: `E(i_flux, j, k, comp) = ...` instead of `E(comp, i_flux, j, k) = ...`.

**Files requiring changes:**

| File | Role |
|------|------|
| `3D_solver/src/calc_steps.f90` | `calc_R` — flux consumer |
| `3D_solver/src/calc_keep_kernel.f90` | KEEP flux producer |
| `3D_solver/src/calc_keep_3d.f90` | KEEP 3D dispatcher |
| `3D_solver/src/calc_slau_kernel.f90` | SLAU flux producer |
| `3D_solver/src/calc_slau_3d.f90` | SLAU 3D dispatcher |
| `3D_solver/src/calc_hybrid_kernel.f90` | Hybrid flux producer |
| `3D_solver/src/calc_hybrid.f90` | Hybrid dispatcher |
| `3D_solver/src/calc_visc2.f90` (if it uses E/F/G) | Viscous flux |
| `3D_solver/src/calc_visc4.f90` | 4th-order viscous |
| `3D_solver/src/calc_flux_base.f90` | Dispatch + allocation calls |
| `src/cpu_gpu_mpi.f90` | Device array allocation |

**Effect on producers (e.g., KEEP kernel):** Currently the KEEP kernel writes all 5 components at once for a single cell: `E(1,i,j,k) = ...; E(2,i,j,k) = ...; ...`. With SoA, it writes `E(i,j,k,1) = ...`. The write access by a single thread is now non-sequential within the E array (components 1–5 are separated by `(nx-1)×(ny-2)×(nz-2)` elements). However, the write occupancy concern is less critical than read coalescing because writes are posted to L2 and combined. Alternatively, the producers could use a register-resident intermediate and issue a single 40-byte store using `type` or vectorized writes.

**Estimated speedup on `calc_step2_3`:** ~3.5× (from ~12.4 ms to ~3.5 ms per call), based on reducing the E/F/G over-fetch ratio from 5× to 1×.

---

### Suggestion 2 — Reduce register count via `attributes(maxregcount=N)` [Quick win, ~15-25% speedup]

The compiler allocates 44 registers/thread but only 32 are needed to increase block occupancy from 10 to 16 blocks/SM (100% theoretical). Add the CUDA Fortran register cap directive to the step kernels:

```fortran
attributes(global) subroutine calc_step2_3(...)
  !$cuf kernel do <<< grid, threads >>>   ! (existing)
  ! OR: add compiler directive at module/routine level
```

In NVHPC CUDA Fortran, per-routine register limiting is set via compile flag `-gpu=maxregcount:32` in the Makefile for the step file only, or using a compiler hint. The NSTGV `Makefile` compiles all sources with the same flags; consider splitting `calc_steps.f90` into its own compile rule with an additional `-gpu=maxregcount:32` flag.

Alternatively, manual register reduction: the `R(5)` array and temporaries like `coef3_dtdxdy/dtdydz/dtdzdx` are candidates — they can be computed inline without storing in named variables, potentially reducing register pressure.

**Tradeoff:** Risk of register spilling to local memory (currently 0 spill). If spilling occurs, local memory is cached in L1 (fast), so net effect is usually still positive. Profile after applying to verify.

**Estimated speedup on `calc_step2_3`:** ~15-25% (reduces the 87% No-Eligible stall rate by providing more in-flight warps to hide DRAM latency).

---

### Suggestion 3 — Mixed precision (real(4)) for E, F, G flux arrays [~2× reduction in bandwidth, medium implementation effort]

The E/F/G flux arrays are computed fresh every RK stage and never accumulated long-term. Single-precision storage halves their size:

- Each flux array: `real(8)` → `real(4)`: 5.4 GB → 2.7 GB
- All three: 16.2 GB → 8.1 GB
- L2 working set halved → L2 hit rate improves from 33% to potentially 50-60%

Implementation: declare E/F/G as `real(4)` device arrays in `cpu_gpu_mpi.f90`; cast to `real(8)` when reading in `calc_R`:

```fortran
! In calc_R, reading as real(4):
real(4), intent(in), device, contiguous :: E(nx-1,ny-2,nz-2,5)
R(1) = dtdydz * (-dble(E(i,j,k,1)) + dble(E(i+1,j,k,1)))
```

**Tradeoff:** Introduces single-precision rounding into the RK update. Acceptable for many CFD applications (turbulence, SBLI) but needs validation. The Q array (solution vector) remains `real(8)`.

**Estimated speedup on `calc_step2_3`:** ~1.5-2× combined with Suggestion 1.

---

## Combined Expected Impact

Applying Suggestions 1 + 2:
- `calc_step2_3`: ~12.4 ms → ~3.0-4.0 ms per call
- Total per TVD RK3 timestep (2 calls): ~24.8 ms → ~6-8 ms

This would likely move the bottleneck to the `calc_visc4_internal_calc_gv4_in_` kernel (currently 13.47 ms/call in pipeline), which has separate optimization opportunities.

---

## Verification

```bash
cd 3D_solver/NSTGV
make clean && make
bash calc.sh          # confirm simulation runs correctly
bash test_cicd.sh     # all 12 build targets pass
pytest ouxsbli/tests/ # correctness regression (especially test_etgv.py)
bash profile.sh       # re-profile to measure speedup
```
