# Bug Comparison Analysis Plan: Non-Curv vs Curv Variants

## Task Overview
Compare 9 non-curv reference files with their curv variants to find bugs introduced or present.

## Reference Files (Non-Curv - Original)
1. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_flux_base.f90
2. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_hybrid.f90
3. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_hybrid_kernel.f90
4. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_keep_kernel.f90
5. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_slau_kernel.f90
6. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_steps.f90
7. /home/jhatayama/ouxsbli/wing/OUxSBLI/3D_solver/src/calc_time_dev.f90
8. /home/jhatayama/ouxsbli/wing/OUxSBLI/src/print.f90
9. /home/jhatayama/ouxsbli/wing/OUxSBLI/src/set_coordinate.f90

## Curv Variants (Derived - May contain bugs)
- calc_flux_base_curv.f90
- calc_hybrid_curv.f90
- calc_hybrid_kernel_curv.f90
- calc_keep_kernel_curv.f90
- calc_slau_kernel_curv.f90
- calc_steps_curv.f90
- calc_time_dev_curv.f90
- print_curv.f90
- (no set_coordinate_curv.f90 found - may not be applicable)

## Analysis Strategy

### Phase 1: Identify Key Differences
For each file pair, use diff to identify:
- Added/removed subroutines
- Modified function signatures
- Changed array dimensions
- Altered calculation logic
- Variable type changes

### Phase 2: High-Risk Areas to Investigate
Based on code structure analysis, focus on:
- **calc_flux_base**: Flux calculation dispatch, array dimensions
- **calc_hybrid**: Ducros sensor computation (lines 12-73)
- **calc_hybrid_kernel**: Hybrid scheme blending logic, sensor thresholds
- **calc_keep_kernel**: KEEP scheme kernels with shared memory
- **calc_slau_kernel**: SLAU scheme with shared memory arrays
- **calc_steps**: Runge-Kutta time integration, accumulation logic
- **calc_time_dev**: Time stepping loop, flux integration
- **print**: Output data dimensions, array indexing

### Phase 3: Specific Bug Patterns to Check
1. **Array indexing mismatches** - Off-by-one errors in curvilinear code
2. **Shared memory array sizing** - Different sizes in x, y, z for curv vs Cartesian
3. **Flux array dimensions** - E(5,nx-1,ny-2,nz-2) vs modified for curv
4. **Jacobian usage** - Different metric handling in curvilinear
5. **Sensor/Threshold application** - Different values/thresholds in curv
6. **Runge-Kutta coefficients** - Any sign/magnitude changes
7. **Boundary conditions** - Different BC handling in curvilinear
8. **Variable type promotions** - Integer vs Real confusion

## Expected Findings
Since curv variants handle curvilinear grids instead of Cartesian:
- Jacobian metrics likely added/modified
- Flux arrays may have different shapes
- Sensor calculations may be different
- Boundary condition handling likely differs
- Coordinate metric terms (xi, eta, zeta) may affect calculations

## Deliverable Format
For each bug found:
1. File name
2. Line numbers
3. Nature of bug (logic error, type mismatch, array bounds, etc.)
4. Code snippet comparison
5. Impact assessment
6. Suggested fix
