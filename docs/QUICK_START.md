# Quick Start: Multi-GPU Optimization Checklist

## 🎯 Executive Summary

This checklist enables overlapped computation and communication in OUxSBLI for efficient multi-GPU execution.

**Expected Benefit**: 
- 2 GPU: +7-10% efficiency
- 4 GPU: +15-20% efficiency  
- 8 GPU: +25-30% efficiency

**Time Investment**: 
- PoC (proof-of-concept): 1-2 weeks
- Full implementation: 6-8 weeks

---

## 📋 Phase 1: Minimal PoC (Week 1-2)

Minimal changes to demonstrate overlap concept:

- [ ] **Diagnostic**: Run nsys profile on current code
  ```bash
  cd 3D_solver/ETGV
  nsys profile -t cuda,mpi mpirun -np 4 ./a.out
  ```
  
- [ ] **Create `calc_halo_exchange.f90`** (Module 2.2 in IMPLEMENTATION_GUIDE.md)
  - ~80 lines of code
  - Non-blocking MPI Isend/Irecv wrapper
  
- [ ] **Add minimal modification to `calc_time_dev.f90`** (Template 2 in IMPLEMENTATION_GUIDE.md)
  ```fortran
  ! Insert before calc_step:
  call post_halo_z(myrank, nranks, nx, ny, nz, 3, QJ, halo_state, ierr)
  
  ! Replace set_bc with wait:
  call calc_step1(...)  ! ← GPU computes while MPI active
  call wait_halo_z(myrank, nx, ny, nz, 3, QJ2, halo_state, ierr)
  ```

- [ ] **Test**: Compare results with original (should be identical)
  ```bash
  mpirun -np 4 ./a.out.old > results.old
  mpirun -np 4 ./a.out.new > results.new
  diff results.old results.new
  ```

- [ ] **Profile**: Measure speedup
  ```bash
  time mpirun -np 4 ./a.out.old  # Note: seconds
  time mpirun -np 4 ./a.out.new
  ```

**Expected Outcome**: ~5-10% speedup with minimal code

---

## 📋 Phase 2: Domain Splitting (Week 3-4)

Add interior/boundary region separation:

- [ ] **Create `calc_domain.f90`** (Module 2.1 in IMPLEMENTATION_GUIDE.md)
  - ~120 lines
  - Define interior vs. boundary regions
  - Test on all ranks

- [ ] **Create `calc_overlap_manager.f90`** (Module 2.3 in IMPLEMENTATION_GUIDE.md)
  - ~80 lines
  - Track timing per phase
  - Report overlap efficiency

- [ ] **Split KEEP flux kernels** (Step 4.1 in IMPLEMENTATION_GUIDE.md)
  - Add `calc_keep_x6_interior`, `calc_keep_y6_interior`, `calc_keep_z6_interior`
  - ~100 lines per direction (copy + modify bounds)
  - Test interior-only computation produces same fluxes

- [ ] **Do similar for SLAU kernels** (if using SLAU scheme)
  - `calc_slau_x6_interior`, etc.

- [ ] **Test interior kernels**
  ```bash
  # Compare flux values on interior with original
  # Should be bitwise identical
  ```

**Expected Outcome**: ~15-20% speedup with proper interior/boundary separation

---

## 📋 Phase 3: Integration & Optimization (Week 5-8)

Full time-stepping restructure with all tricks:

- [ ] **Restructure main RK loop** (Step 3.1 in IMPLEMENTATION_GUIDE.md)
  - 3-phase pattern: compute interior → post async → compute boundary → wait
  - ~50 lines modification

- [ ] **Optimize GPU memory layout**
  - Ensure boundary planes are contiguous
  - Enable 1-shot `cudaMemcpy`

- [ ] **Add CUDA streams** (Module Step 10 in IMPLEMENTATION_GUIDE.md)
  - Overlap kernel launches with data transfers

- [ ] **Persistent MPI** (Optional, for repeated patterns)
  - Use `MPI_Send_init` / `MPI_Start` instead of `MPI_Isend`

- [ ] **Test on varied GPU counts**
  - 2 GPU ✓
  - 4 GPU ✓
  - 8 GPU ✓
  - 16 GPU (if available) ✓

- [ ] **Profile & document**
  - nsys: Verify GPU kernel overlaps with MPI
  - ncu: Check GPU utilization > 90%
  - Create performance report

**Expected Outcome**: 25-30% speedup at scale

---

## 🔧 Quick Reference: File Changes

```
NEW FILES TO CREATE:
├── src/calc_domain.f90                    (Step 2.1)
├── src/calc_halo_exchange.f90             (Step 2.2)
├── src/calc_overlap_manager.f90           (Step 2.3)
├── src/calc_flux_interior.f90             (Step 4.1)
├── src/calc_flux_boundary.f90             (Step 4.1)
└── docs/IMPLEMENTATION_GUIDE.md           (Already created)

EXISTING FILES TO MODIFY:
├── 3D_solver/src/calc_time_dev.f90        (Step 3.1, ~80 lines changed)
├── src/calc_flux_base.f90                 (Step 4.2, ~50 lines added)
├── src/calc_keep_kernel.f90               (Step 4.1, ~100 lines added)
┐── src/calc_slau_kernel.f90              (Step 4.1, ~100 lines added)
└── Makefile                               (Add new source files, <10 lines)

DOCUMENTATION:
├── docs/MULTI_GPU_SCALABILITY_PLAN.md     (Strategic overview)
├── docs/IMPLEMENTATION_GUIDE.md           (Detailed implementation)
└── docs/QUICK_START.md                    (This file)
```

---

## 🚀 Starting Commands

### Option A: Copy-Paste PoC
```bash
# 1. Back up current working version
cd /home/jhatayama/OUxSBLI
git checkout -b feature/overlapped-compute

# 2. Create minimal halo exchange module
cp docs/IMPLEMENTATION_GUIDE_src/calc_halo_exchange.f90 src/

# 3. Modify calc_time_dev.f90 (use sed or manual edit)
# ... see Step 3.1 template ...

# 4. Test
cd 3D_solver/ETGV
make clean && make
mpirun -np 4 ./a.out
```

### Option B: Guided Development
1. Start with IMPLEMENTATION_GUIDE.md Part 1 (Diagnostic)
2. Follow each step sequentially
3. Test after each commit
4. Profile regularly

---

## 📊 Performance Checklist

After each phase, validate:

- [ ] **Correctness**
  ```bash
  diff <(standard_output) <(optimized_output)  # Should be very close or identical
  ```

- [ ] **GPU Utilized**
  ```bash
  nsys profile -t cuda ...
  # Check: Kernel execution continuous (no long idle gaps)
  ```

- [ ] **MPI Not Blocking**
  ```bash
  # In timeline, see MPI calls (blue bars) overlapped with kernels (green bars)
  ```

- [ ] **Speedup Measured**
  ```bash
  time mpirun -np 4 ./a.out.original
  time mpirun -np 4 ./a.out.optimized
  # Compute: (t_original - t_optimized) / t_original * 100 = X% speedup
  ```

---

## 🐛 Troubleshooting Map

| Problem | Solution | See |
|---------|----------|-----|
| Code crashes after changes | Check domain bounds in interior/boundary split | Part 7, Issue 1 |
| No speedup observed | Verify overlap fraction > 30% in logs | Part 7, Issue 3 |
| Out of memory errors | Reduce pinned buffer sizes or use CUDA-aware MPI | Part 7, Issue 2 |
| Different results | Ensure interior/boundary kernels return same fluxes | Unit test 5.1 |
| MPI hangs | Verify MPI_Isend/Irecv match (send size = recv size) | Step 3.1 |

---

## 📈 Expected Output

### After PoC (Week 2):
```
Rank 0: Interior=10.2ms, Boundary=4.5ms, Comm=8.1ms
  Estimated overlap fraction: 55.6%
  → ~5-10% speedup
```

### After Full Optimization (Week 8):
```
Rank 0: Interior=9.8ms, Boundary=4.2ms, Comm=8.0ms
  Estimated overlap fraction: 85.2%
  → 25-30% speedup at scale
```

---

## 📚 Documentation Structure

```
docs/
├── api.md                              (API reference - auto-generated)
├── technical_doc.md                    (Existing technical documentation)
├── MULTI_GPU_SCALABILITY_PLAN.md       (Strategic architecture + 10 steps)
├── IMPLEMENTATION_GUIDE.md             (Detailed code walkthroughs + templates)
└── QUICK_START.md                      (This file - checklist + commands)
```

**Read in order**:
1. This file (QUICK_START.md) - 5 min overview
2. MULTI_GPU_SCALABILITY_PLAN.md - Understand strategy (30 min)
3. IMPLEMENTATION_GUIDE.md - Implement step-by-step (ongoing)

---

## 🎓 Key Concepts

**Interior Compute** (safe, no neighbor data):
- Fluxes at non-boundary locations
- Start immediately after previous step
- Can run in parallel with MPI

**Boundary Regions** (depends on neighbor data):
- Fluxes at domain edges (z-direction boundaries)
- Must wait for MPI recv before computing
- But can overlap with interior compute

**Async Halo Exchange** (non-blocking MPI):
- `MPI_Isend` starts transfer, returns immediately
- GPU continues compute while network transfers data
- `MPI_Wait` blocks only if transfer not complete

**Overlap Fraction**:
- `(boundary_compute_time)` / `(halo_wait_time)`
- Best case: 100% (all boundary compute happens during MPI wait)
- Realistic: 50-90% depending on problem size

---

## ✅ Success Criteria

Your optimization is successful when:

1. ✅ **Correctness**: Solution bitwise identical or within machine epsilon
2. ✅ **Speedup**: Measured 10%+ improvement on 4+ GPU case
3. ✅ **Scaling**: Linear weak scaling up to 16 GPUs
4. ✅ **GPU Util**: >90% GPU utilization in profiler timeline
5. ✅ **No Hangs**: All tests run to completion
6. ✅ **Documented**: Code comments explain overlap logic

---

## 📞 Getting Help

**If stuck on**:
- Module creation → See Part 2 (calc_domain.f90 example)
- Kernel splitting → See Part 4 (interior/boundary signature)
- Debugging → See Part 7 (troubleshooting issues)
- Performance → See Part 5 & 6 (validation/profiling)

**Generate diagnostic info**:
```bash
# Send input to forum/issue:
nsys profile -t cuda,mpi -o diag_ mpirun -np 4 ./a.out
gzip diag_*.nsys-rep
# Attach to issue with error message
```

---

**Version**: 1.0 (March 31, 2026)  
**Status**: Ready for implementation  
**Next**: Choose PoC or Guided Development path above
