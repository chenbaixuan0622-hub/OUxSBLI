# OUxSBLI Project Context

## 1. Introduction & Overview

**OUxSBLI** stands for **O**kayama **U**niversity × **SBLI** (Shock-Boundary Layer Interaction). It is a GPU-accelerated Computational Fluid Dynamics (CFD) solver designed to study compressible flows, particularly shock-boundary layer interactions, using modern GPU computing and high-order finite-difference methods.

The project represents an effort of students at the University of Okayama to develop state-of-the-art tools for studying complex aerodynamic and hypersonic phenomena.

---

## 2. Scientific Background

### 2.1 Governing Equations

OUxSBLI solves the **3D compressible Navier-Stokes equations** (or Euler equations in inviscid mode):

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial \rho \mathbf{u}}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) + \nabla p = \nabla \cdot \boldsymbol{\tau}$$

$$\frac{\partial E}{\partial t} + \nabla \cdot (E \mathbf{u} + p \mathbf{u}) = \nabla \cdot (\boldsymbol{\tau} \cdot \mathbf{u}) + \nabla \cdot \mathbf{q}$$

Where:
- $\rho$: fluid density
- $\mathbf{u}$: velocity vector
- $p$: pressure
- $E$: total energy per unit volume
- $\boldsymbol{\tau}$: viscous stress tensor
- $\mathbf{q}$: heat flux vector

### 2.2 Research Focus: Shock-Boundary Layer Interaction

**SBLI** occurs when a shock wave impinges on a boundary layer, creating complex three-dimensional separated flow regions. This phenomenon is critical in:

- **Hypersonic Aerodynamics**: Shock-induced separation on aircraft fuselages & wings
- **Rocket Engines**: Supersonic nozzle shock structures
- **Supersonic Inlets**: Flow control and starting problems
- **Ramjet/Scramjet Propulsion**: Inlet supersonics flow

**Challenges**:
- Unsteady shock motion
- Turbulence generation and amplification
- Multi-scale interactions (shock curvature ↔ boundary layer thickness)
- Accurate prediction requires DNS or LES with shock-capturing schemes

---

## 3. Computational Approach

### 3.1 Why GPU Acceleration?

**Problem Scale**:
- Production SBLI simulations: 20M-200M grid points
- Hundreds to thousands of time steps
- Multiple schemes/accuracies for convergence studies

**Traditional CPU Performance**:
- Single Intel Xeon: ~200-500 GFLOPS
- Full SBLI simulation: weeks to months

**GPU Acceleration**:
- Single RTX 4090: ~1000+ TFLOPS (40× speedup)
- SBLI simulation: hours to days
- Enables rapid iteration for algorithm development

### 3.2 High-Order Methods Rationale

**Why not 2nd-order?**
- 2nd-order requires very fine grids for accuracy
- Smooth flows: KEEP/SLAU(4th) can use 5× coarser grids than 2nd-order
- DNS/LES: High accuracy needed for turbulence spectrum

**Efficiency Trade-off**:
- Wider stencil (more memory bandwidth)
- Larger timestep (higher accuracy → larger CFL)
- Net result: 2-4× speedup vs. 2nd-order at equivalent accuracy

### 3.3 Scheme Selection Philosophy

| Scheme | Purpose | Accuracy | Cost | Best For |
|--------|---------|----------|------|----------|
| **KEEP** | Kinetic energy preserving | 2nd/4th/6th | Low | Smooth turbulence, DNS |
| **SLAU** | Low dissipation, shock-adaptive | 2nd/4th/6th | Medium | SBLI, mixed smooth+shock |
| **Roe** | Classical, robust | 2nd order | Low | Hypersonic, extreme shocks |
| **Hybrid** | Adaptive switching | Variable | Variable | General-purpose |

---

## 4. Use Cases

### 4.1 Fundamental Studies

#### Taylor-Green Vortex (TGV)
```
Purpose: Validate code, verify convergence rates
Setup: Periodic domain, smooth initial vortex
Cases:
  - ETGV: Euler TGV (compressible, inviscid)
  - NSTGV: Navier-Stokes TGV (with viscosity, LES)
Output: Kinetic energy decay curves, spectral analysis
```

**Why TGV?**
- Analytical solution exists (before shock formation)
- Clean test of advection, viscosity, LES models
- Independent of scheme (all should show same decay rate)
- Sensitivity to grid resolution and accuracy order

#### Isentropic Vortex (IVST)
```
Purpose: Validate convection schemes
Setup: Convecting vortex in uniform compressible flow
Target: Low numerical error over many wavelengths
```

#### Kelvin-Helmholtz Instability (KHI)
```
Purpose: Test turbulence transition, mixing layer physics
Setup: Two counter-flowing streams, shear-induced rolls
Features: Vortex rollup, pairing, small-scale generation
```

### 4.2 Applied Research

#### Turbulent Boundary Layer (TBL)
```
Purpose: Study wall-bounded turbulence
Setup: Channel or flat-plate configuration
Conditions: Mach 2-5, Reτ ~ 200-1000
Output: Mean profiles, Reynolds stress, spectra, structures
```

#### Shock-Boundary Layer Interaction (SBLI)
```
Purpose: Primary focus - understand SBLI physics
Setup: Oblique shock impinging on turbulent BL
Conditions: Shock-induced separation region
Physics: Unsteady shock motion, vorticity generation, BL destabilization
Output:
  - Mean flow topology (separation bubble)
  - Pressure-induced fluctuations
  - Supersonic coherent structures
  - Shock motion statistics
```

---

## 5. Technical Capabilities

### 5.1 Discretization Methods

**Spatial (Convection)**:
- KEEP scheme: Kinetic energy & entropy preserving (no artificial dissipation)
- SLAU scheme: Low-dissipation AUSM with shock-adaptive interpolation
- Roe scheme: Classical flux-splitting for shocks
- Accuracy: 2nd, 4th, 6th-order finite-difference

### 5.2 Parallel Computing

**GPU Support**:
- NVIDIA HPC SDK (nvfortran compiler)
- CUDA Fortran GPU kernels
- Shared memory optimization for stencil operations
- Coalesced global memory access patterns

**MPI Parallelization**:
- 1D domain decomposition (z-direction)
- GPU-aware MPI (direct GPU-to-GPU communication)
- Asynchronous ghost region exchange
- Fallback to pinned-memory MPI for compatibility

**Scalability**:
- Single GPU: Full problem on RTX 4090 (8GB → 256M+ points)
- Multi-GPU: Up to 8 GPUs per node (typical HPC configs)
- Production: 100+ ranks on GPU clusters

### 5.3 Test Cases

| Case | Equations | Grid Size | Domain | Input Files |
|------|-----------|-----------|--------|-------------|
| ETGV | Euler | 513³ | $[0,2\pi]^3$ | 3D_solver/ETGV/ |
| NSTGV | NS | 513³ | $[0,2\pi]^3$ | 3D_solver/NSTGV/ |
| IVST | Euler | 256³ | $[0,25]×[0,25]×[0,50]$ | 3D_solver/IVST/ |
| KHI | NS | 256³ | Domain-specific | 3D_solver/KHI/ |
| TBL | NS | 256×512×128 | Channel geometry | 3D_solver/TBL/ |
| SBLI | NS (dual-region) | 512×350×256 | BL + shock | 3D_solver/SBLI/ |

---

## 6. Development History & Status

### 6.1 Evolution

**Phase 1: Foundation (2020-2021)**
- Initial CUDA Fortran port of serial CFD solver
- 2nd-order schemes only
- Single-GPU support
- TGV test cases

**Phase 2: Expansion (2021-2023)**
- Added 4th/6th-order KEEP scheme
- SLAU scheme with shock sensing
- Multi-GPU with MPI
- LES turbulence model (under dev)
- SBLI dual-region solver

**Phase 3: Optimization (2023-2026)**
- Refactored kernels: calc_keep_kernel_internal, calc_slau_kernel_internal
- Performance profiling and tuning
- Memory efficiency improvements
- Extended accuracy validation

### 6.2 Recent Improvements (April 2026)

**Kernel Optimizations**:
- Unified multi-accuracy KEEP kernels (code consolidation)
- Smart shared-memory reuse in SLAU (stack frame avoidance)
- Wiggle detector for shock-adaptive dissipation
- 10-15% speedup on typical workloads

**Code Quality**:
- Removed unused variables
- Cleaned up build scripts
- Improved documentation

---

## 7. Scientific Goals

### 7.1 Short-term (0-2 years)

- ✅ Validate KEEP/SLAU schemes on smooth test cases (TGV, IVST)
- ✅ Demonstrate SBLI capability and accuracy improvements over industrial codes
- 🔄 Extend LES model validation on TBL cases
- 🔄 Multi-GPU scaling studies
- Algorithm advancement
- Large-scale production simulations (100+ GPU hours)
- Integration with aerodynamic design optimization

---

## 8. Why OUxSBLI?

### 8.1 Advantages vs. Alternatives

| Feature | OUxSBLI | Commercial (ANSYS FUN3D/CFX) | OpenFOAM | FLASH |
|---------|---------|--------------------------------|----------|-------|
| **GPU Native** | ✅ Full CUDA | ⚠️ Partial | ❌ Limited | ⚠️ Some |
| **High-Order** | ✅ 2-6th | ⚠️ up to 3rd | ⚠️ up to 2nd | ✅ AMR-native |
| **Compressible RMM** | ✅ KEEP/SLAU | ⚠️ Roe-only | ❌ No | ✅ PPM-MUSCL |
| **API Exposed** | ✅ Open-source | ❌ Closed | ✅ Open | ✅ Open |
| **Time-to-Solution** | ✅ Hours | ❌ Days-weeks | ⚠️ Days | ⚠️ Variable |
| **Cost** | ✅ Free | ❌ $$$$$ | ✅ Free | ✅ Free |

---

## 9. Community & Collaboration

### 9.1 Current Users

- University of Okayama : Primary developers & users

### 9.2 Contributing & Getting Involved

**For Researchers**:
1. Clone/fork the repository
2. Set up CUDA Fortran environment (HPC SDK 24.* or 25.*)
3. Run test cases on your GPU hardware
4. Contribute new schemes/cases via pull requests

**For Students**:
- Excellent vehicle for learning:
  - GPU computing (CUDA Fortran)
  - CFD methods (schemes, time-stepping)
  - HPC optimization (memory, communication)
  - Compressible flow physics

---

## 10. Computational Requirements & Hardware

### 10.1 Minimum Setup

```
GPU: RTX 4060 (8 GB) or equivalent
CPU: AMD or Intel (for compilation & I/O)
RAM: 32 GB (host)
MPI: NVIDIA HPC-X (GPU-aware MPI)
Storage: 128 GB (for medium simulations)
```

---

## 11. Getting Started

### 11.1 Quick Start

```bash
# Clone repository
git clone https://github.com/htymjun/OUxSBLI.git
cd OUxSBLI

# Load HPC SDK
module load hpc-sdk/24.x

# Compile test case
cd 3D_solver/ETGV
make

# Run simulation
bash calc.sh

# View output
paraview data/Q*.vtr
```

### 11.2 Documentation

- **technical_doc.md**: Detailed CFD/algorithm documentation
- **architecture.md**: System design and data flow
- **api.md**: Function/subroutine signatures
- **this file (context.md)**: Project background and research context

---

## 12. Support & Contact

For questions, bug reports, or collaboration inquiries:

**Email**: (To be provided by institution)  
**Issues/Bugs**: GitHub issues (https://github.com/htymjun/OUxSBLI/issues)  
**Discussions**: GitHub discussions or mailing list

---

**Document Version**: 1.0 (April 5, 2026)  
**Last Updated**: Project context and scientific background comprehensive review  
**Status**: Active development with regular performance optimizations
