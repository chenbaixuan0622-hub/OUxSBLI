# Markdown Documentation Rules for OUxSBLI

This document defines strict rules for maintaining three distinct markdown documentation files: `architecture.md`, `context.md`, and `technical_doc.md`. Each serves a specific purpose and audience.

---

## Overview Table

| Aspect | architecture.md | context.md | technical_doc.md |
|--------|-----------------|------------|------------------|
| **Primary Audience** | AI agents | AI agents | Human developers |
| **Purpose** | System structure & design | Implementation details | Technical reference |
| **Format** | Compact, structured | Detailed organized | Prose + diagrams |
| **Update Frequency** | When architecture changes | When code updates | When methods change |
| **Include Functions?** | Module names only | Yes, all subroutines | Yes, key ones only |
| **Include Physics?** | No | Rarely | Yes, detailed |
| **Code Examples?** | Pseudocode only | Real code snippets | Equations & concepts |

---

## 1. architecture.md - Compact AI-Optimized Design Document

### 1.1 Purpose
- **Audience**: AI agents (Copilot, Claude, etc.)
- **Use Case**: AI understand system structure, navigation, and design patterns in minimal tokens
- **Goal**: Enable AI to reason about module relationships without verbose prose

### 1.2 Content Structure (Required Sections)

#### **1.2.1 Top-Level Architecture Diagram** (Required)
- Visual representation of 3-5 key layers
- Use ASCII art or text-based block diagrams
- Show data flow: Input → Processing → Output
- **Constraint**: Max 20 lines of diagram

```
EXAMPLE:
┌─────────┐
│ Input Q │
├─────────┤
│ Flux    │
│ Compute │
├─────────┤
│ BC      │
│ Update  │
├─────────┤
│Output Q │
└─────────┘
```

#### **1.2.2 Module Dependency Graph** (Required)
- List each module and its **direct dependencies only**
- Format: `module_name: [dep1, dep2, dep3]`
- Do NOT list transitive dependencies
- Keep list to ≤ 30 modules

```
EXAMPLE:
calc_time_dev.f90: [calc_physical_quantities, calc_keep_kernel_internal, calc_visc2, set_bc_common]
calc_keep_kernel_internal.f90: [mod_constant, mod_globals]
calc_slau_kernel_internal.f90: [mod_constant, mod_globals, calc_muscl]
```

#### **1.2.3 Module Categories (Tier Classification)** (Required)
- Group by tier (1-7)
- Each module: name, file path, brief 1-line purpose
- **Format**: `Tier N: [Module Name] (file_path) - Purpose`
- **Constraint**: 1-2 lines per module max

```
EXAMPLE:
Tier 1: [mod_globals] (3D_solver/ETGV/mod_globals.f90) - Grid size and scheme parameters
Tier 1: [set] (3D_solver/ETGV/set.f90) - Grid initialization and BC setup
Tier 3: [calc_keep_kernel_internal] (3D_solver/src/calc_keep_kernel_internal.f90) - Multi-accuracy KEEP flux kernels
```

#### **1.2.4 Data Structures** (Required)
- **Only** list primary data structures passed between modules
- Format: `name(dimensions) : type : description`
- **Constraint**: Max 15 data structures, only those AI needs to know

```
EXAMPLE:
Q(5,nx,ny,nz) : real(8), device : conservative variables [ρ, ρu, ρv, ρw, E]
E(5,nx-1,ny-2,nz-2) : real(8), device : x-direction convection flux
T(nx,ny,nz) : real(8), device : temperature field
detector(nx,ny,nz) : real(8), device : shock detection sensor
```

#### **1.2.5 Design Patterns** (Required)
- Max 4 key patterns
- Each: pattern name, 2-3 line description, 1 tiny code example

```
EXAMPLE:
**Pattern: Parameter-based Accuracy Selection**
  Use kind(id_accuracy) to select 2nd/4th/6th order at compile time.
  io = kind(id_accuracy) / 3 → 0, 1, or 2 for stencil size
  
**Pattern: Shared Memory Tiling**
  Load stencil data with overlap, sync threads, avoid global memory in inner loop
  Reduces bandwidth by 10-100x
```

#### **1.2.6 Execution Flow** (Required)
- Single timestep sequence using bullet points
- Show 6-10 key steps only (no sub-steps)
- Include stage count and sync points

```
EXAMPLE:
1. [Stage 1] calc_physical_quantities(Q → T, μ)
2. [Stage 1] calc_keep_x_in, _y_in, _z_in (Q, T → E, F, G)
3. [Stage 1] calc_visc2 (Q, T, μ → ΔE, ΔF, ΔG)
4. [Stage 1] RHS computation & Q update
5. [Stage 1] set_bc_cyclic + MPI_exchange (async)
6. (Repeat Stages 2, 3 for RK scheme)
```

### 1.3 Style Rules for architecture.md

- **Language**: English, imperative mood ("Add X", "Compute Y")
- **Verbosity**: Absolute minimum - 1 sentence per concept
- **Headings**: Use #, ##, ### only (max depth 3)
- **Lists**: Bullet form, no elaboration
- **Code**: Pseudocode ONLY (not real Fortran) - show intent, not implementation
- **Jargon**: OK to use (GPU, shared memory, stencil, SLAU, KEEP)
- **Line Length**: Keep under 100 chars where possible
- **NO**: 
  - Long prose paragraphs
  - References to papers
  - Detailed physics equations
  - Full function signatures
  - Historical context

### 1.4 Maintenance Rules for architecture.md

- **When to Update**:
  - New module added/removed
  - Data structure changed
  - Kernel renamed or refactored
  - Execution flow fundamentally altered

- **NOT Updated For**:
  - Bug fixes
  - Performance tweaks
  - Minor variable name changes
  - Documentation improvements alone

- **Validation**: Can condensed version fit on 2-3 printed pages? If no, trim.

---

## 2. context.md - Detailed Implementation Context for AI Agents

### 2.1 Purpose
- **Audience**: AI agents (Copilot, Claude) during deep code understanding tasks
- **Use Case**: AI needs detailed knowledge of what each function/subroutine does, why it exists, and how it connects to physics
- **Goal**: Provide AI sufficient context to modify, extend, or explain implementations

### 2.2 Content Structure (Required Sections)

#### **2.2.1 Introduction & Project Context** (Required)
- 3-4 paragraphs explaining:
  - What is OUxSBLI and its purpose
  - Primary use cases (SBLI, TBL, etc.)
  - Why this approach (GPU, high-order, etc.)
- **Target**: Help AI understand "game"

#### **2.2.2 Scientific Background (Brief)** (Required)
- Governing equations: Euler/NS in conservative form
- Key physics: shock, boundary layer, turbulence
- Why this matters: 1-2 paragraphs max
- **NO**: Derivations, deep theory, research papers

#### **2.2.3 Module-by-Module Implementation Details** (Required)
**For EACH major module (20-30 total), provide:**

```
### mod_globals (3D_solver/ETGV/mod_globals.f90)

**Purpose**: Hold simulation parameters that control problem configuration

**Key Parameters**:
- id_visc (integer): 2=Navier-Stokes, 4=Euler
- id_scheme (real): 0=KEEP, 2=SLAU, 4=Roe
- id_accuracy (parameter by kind): 2/4/8 for 2nd/4th/6th order
- nx, ny, nz: grid dimensions (INTEGER, parameter)
- Lx, Ly, Lz: domain sizes (REAL(8), parameter)
- threadsE, threadsF, threadsG (type(dim3)): GPU thread block config
- dt: timestep (REAL(8), parameter)

**Key Functions/Subroutines**: None (this is parameter-only module)

**Dependencies**: cudafor only

**Used By**: Every simulation module (calc_time_dev, calc_keep_*, etc.)

**Example Usage**:
  if (kind(id_accuracy) == 4) then
    ! 4th-order accuracy selected - use 4-point stencil
  endif

**Important Notes**:
- Changed at compile-time only (rebuilding required)
- Thread block dimensions must be tuned per GPU/problem size
- Mismatch between id_visc and scheme can cause artifacts
```

**REQUIRED fields for each module:**
1. File path
2. Purpose (1-2 sentences)
3. Key parameters/variables (if parameters module)
4. Key functions/subroutines (names + 1-line description per subroutine)
5. Dependencies (direct only)
6. Used by (which other modules call this one)
7. Example usage or pseudocode snippet
8. Implementation notes or gotchas

#### **2.2.4 Subroutine-Level Details** (Required for >10 subroutines per module)

For modules with many subroutines, add subsection:

```
### calc_time_dev.f90 Subroutines

**RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)**
- Purpose: 3-stage TVD RK3 time stepping (low-storage variant)
- Inputs: Q(5,nx,ny,nz) conservative state, grid geometry, metadata
- Outputs: Q updated by 3 RK stages
- Key Steps:
  1. Repeat 3 times: calc_physical_quantities → calc_EFG → RHS → UpdateQ
  2. After stage 1&2: reduced coefficient interpolation
  3. After stage 3: full update with TVD property guaranteed
- GPU Kernels Called: (names, no args)
- MPI Calls: Async send/recv of ghost regions
- Notes: Always stable for CFL ≤ 1.0; less accurate than RK4

**RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)**
- Purpose: 4-stage classical RK4 time stepping
- ... (similar structure)
```

#### **2.2.5 Data Flow Diagrams** (Optional but Recommended)
- Text-based diagrams showing data movement between modules
- Example: `Q → [calc_physical_quantities] → T, μ → [calc_keep_x_in] → E`
- Max 1-2 per major section

#### **2.2.6 Integration Examples** (Required)
For AI to understand how to call/modify, provide "integration snippets":

```
**Integration Example: Adding a New Boundary Condition**

1. Create set_bc_mybc.f90 with subroutines:
   - set_bc_mybc_init(nx, ny, nz, Q) : CPU version for IC
   - set_bc_mybc(myrank, nx, ny, nz, Jacobian, Q) : GPU kernel

2. Add to calc_time_dev.f90 after each stage:
   call set_bc_mybc(myrank, nx, ny, nz, Jacobian, Q)

3. Update Makefile to include set_bc_mybc.f90

4. Test: Run existing test case (TGV) and verify Q conservation
```

### 2.3 Style Rules for context.md

- **Language**: English, semi-technical, explain-to-AI tone
- **Verbosity**: Moderate - enough detail for AI to understand, but not novel-length
- **Structure**: 
  - Consistent module template (Purpose, Inputs, Outputs, Steps, Notes)
  - Subsections for subroutines when >5 per module
  - Integration examples as code-like text blocks
- **Code**: 
  - Real Fortran for subroutine signatures
  - Pseudocode/comments for logic flow
  - Full function prototype with intent/type
- **Lists**: Mix of bullets and structured tables
- **Headings**: Use #, ##, ###, #### (max 4 levels)
- **References**: OK to cite "see architecture.md for data structures"
- **Physics**: Keep minimal - focus on implementation, not derivation
- **NO**:
  - Lengthy prose paragraphs
  - Full derivations or papers
  - Historical context or evolution
  - Performance metrics (save for technical_doc.md)
  - Human-only discussion (e.g., "we found...")

### 2.4 Coverage Rules for context.md

**MUST Document:**
- All modules in 3D_solver/src/ (core solvers)
- All case-specific set.f90 files (per case: ETGV, SBLI, etc.)
- All mod_globals.f90 parameter sets
- All calc_*_kernel*.f90 files

**SHOULD Document:**
- Major time-integration functions
- All convection/viscous scheme bases
- Boundary condition implementations

**OPTIONAL:**
- Helper functions (utility_sort, print_debug, etc.)
- Auto-generated files
- Temporary/debugging code

### 2.5 Maintenance Rules for context.md

- **When to Update**:
  - New subroutine added
  - Function signature changed (intent, type, dimensions)
  - Major algorithm change in a subroutine
  - New module created

- **NOT Updated For**:
  - Internal variable renames
  - Comment-only changes
  - Whitespace/formatting
  - Minor logic refactors (same input/output behavior)

- **Validation**: Can an AI understand how to call and modify this code? Yes → OK.

---

## 3. technical_doc.md - Human-Readable Technical Reference

### 3.1 Purpose
- **Audience**: Human developers (researchers, postdocs, interns)
- **Use Case**: Understand CFD methods, numerical schemes, performance considerations
- **Goal**: Comprehensive technical reference for implementation decisions and debugging

### 3.2 Content Structure (Required Sections)

#### **3.2.1 Project Overview** (Required)
- What is OUxSBLI (3-4 paragraphs)
- Key capabilities and use cases
- Target audience and skill level
- **Tone**: Professional research

#### **3.2.2 Computational Foundations** (Required)
- Governing equations (Euler, Navier-Stokes) in conservation form
- Physical parameters (γ, R, Pr, Ma, Re, etc.)
- Initial/boundary condition types
- **Include**: LaTeX equations (allowed in technical_doc only)

#### **3.2.3 Spatial Discretization** (Required)
For EACH major scheme (KEEP, SLAU, Roe, etc.):

```
### 3.3.1 KEEP (Kinetic Energy and Entropy Preserving) Scheme

**Mathematical Formulation**:
  Split convective flux into energy-preserving form:
  E = E_adv + E_press where E_adv preserves kinetic energy invariant
  
**Order**: 2nd, 4th, 6th-order finite difference options

**Implementation Files**:
  - calc_keep_3d.f90: Flux functions KEEP2, KEEP4, KEEP6
  - calc_keep_kernel_internal.f90: GPU kernels calc_keep_x_in, etc.

**Accuracy Properties**:
  - Converges at order p for smooth flows (p=2,4,6)
  - Maintains KE invariant for inviscid flows (test: TGV)
  - Stable without explicit dissipation on smooth grids

**Performance Characteristics**:
  - Flops/point: ~200 (2nd) to ~400 (6th)
  - Memory B/W: ~100 GB/s typical on RTX 4090
  - Speedup vs 2nd-order: 0.8x (higher cost, larger timesteps compensate)

**When to Use**:
  - Clean inviscid flows (TGV, IVST)
  - DNS/LES of smooth turbulence
  - Research on scheme properties

**Gotchas**:
  - Requires (io+1) guard points on each boundary
  - Can be dissipative near shocks without sensor
  [reference to context.md for implementation details]
```

#### **3.2.4 Temporal Discretization** (Required)
- RK3-TVD: stability properties, CFL limits, TVD proof
- RK4: classical formulation, stability region, CFL ~0.8
- Rescaling variants: purpose in SBLI

#### **3.2.5 GPU Implementation Details** (Required)
- Thread organization (1D, 2D, 3D blocks)
- Shared memory strategy (padding, bank conflicts)
- Global memory access patterns (coalescing, strides)
- Example optimization techniques with code snippets

#### **3.2.6 MPI Parallelization** (Required)
- Domain decomposition strategy
- Ghost region sizes per scheme/accuracy
- Synchronization points
- Scaling limits and load imbalance concerns

#### **3.2.7 Convergence Studies** (Required)
- Grid convergence rates (graphs, tables)
- Accuracy vs. grid refinement
- Example: TGV kinetic energy decay rates

#### **3.2.8 Performance Profiling & Optimization** (Required)
- Typical timings (flop count, memory bandwidth, latency)
- Profiling with nsys, ncu commands
- Common bottlenecks
- Case studies: speedup from optimizations

#### **3.2.9 Debugging & Troubleshooting** (Required)
- Common error messages and how to fix
- Validation tests (TGV energy, shock sensor output)
- Memory issues (insufficient GPU RAM, allocations)
- MPI communication deadlocks

#### **3.2.10 Appendices** (Optional)
- Glossary of terms
- Mathematical symbol definitions
- Physical constant values
- Fortran syntax cheat sheet

### 3.3 Style Rules for technical_doc.md

- **Language**: English, technical, written for domain experts
- **Verbosity**: As detailed as needed for complete understanding
- **Tone**: Objective, educational, research-level
- **Structure**:
  - Prose paragraphs with explanations
  - Equations with justification
  - Code snippets for key algorithms
  - Figures/diagrams (ASCII in doc, reference externals if needed)
- **Math**: 
  - LaTeX allowed: `$E = mc^2$` or `$$\frac{\partial Q}{\partial t} = ...$$`
  - Define all symbols before use
  - Cite references when using classical results
- **Code**: 
  - Real Fortran when illustrating algorithms
  - Pseudocode for high-level logic
  - Include accompanying explanation
- **References**: 
  - Cite papers for mathematical schemes
  - Cross-reference context.md for implementation details
  - Cross-reference architecture.md for system structure
- **Headings**: Use #, ##, ###, #### as needed (depth 4)
- **Tables**: Use for results, comparisons, scheme properties
- **Figures**: ASCII diagrams inline, reference external image files in subdirs
- **ALLOWED**:
  - Detailed physics and derivations
  - Historical context and evolution
  - Author thoughts ("this approach avoids XYZ because...")
  - Performance comparisons
  - Research insights

### 3.4 Coverage Rules for technical_doc.md

**MUST Include:**
- Complete mathematical formulations
- All major schemes (KEEP, SLAU, Roe, etc.)
- Temporal and spatial discretization theory
- GPU/MPI algorithmic details
- Convergence and error analysis

**SHOULD Include:**
- Debugging procedures
- Performance benchmarks
- Validation test cases
- Optimization techniques

**OPTIONAL:**
- Detailed proofs (cite papers instead)
- Extensive code listings (show key snippets only)
- Extensive tables of results

### 3.5 Maintenance Rules for technical_doc.md

- **When to Update**:
  - New spatial/temporal scheme added
  - Performance characteristics change significantly
  - Mathematical formulation corrected
  - Error analysis or convergence rates change
  - Debugging procedure discovered/refined

- **NOT Updated For**:
  - Internal implementation refactors (unless changes external interface)
  - Clarifications of code logic alone (update context.md instead)
  - Minor typo fixes

- **Validation**: Can a researcher understand the math and validate implementation? Yes → OK.

---

## 4. Cross-Reference Rules

### 4.1 Linking Between Documents

**From architecture.md:**
- "See context.md for implementation details"
- "See technical_doc.md section X.X for theory"

**From context.md:**
- "Architecture overview in architecture.md"
- "Mathematical details in technical_doc.md section X"
- "For API signatures, see api.md"

**From technical_doc.md:**
- "Implementation in context.md module section"
- "System structure in architecture.md"
- "Full API signatures in api.md"

### 4.2 Avoiding Duplication

| Content | Goes To | Why |
|---------|---------|-----|
| Mathematical derivation | technical_doc.md | Humans need it, AI doesn't |
| Module dependencies | architecture.md | Compact, AI needs it |
| Function signatures | context.md + api.md | Not both; api.md is auto-generated |
| Implementation details | context.md | AI needs to understand code |
| Performance metrics | technical_doc.md | Humans interpret, discuss |
| Execution flow | architecture.md | AI needs compact form |

---

## 5. Update Workflow

### 5.1 When Making Code Changes

1. **Code Change**: Modify Fortran files
2. **Update architecture.md?** 
   - If module added/removed: YES
   - If execution flow changed: YES
   - Otherwise: NO
3. **Update context.md?** 
   - If subroutine signature changed: YES
   - If major logic changed: YES
   - If new function added: YES
   - Otherwise: NO
4. **Update technical_doc.md?** 
   - If algorithm/scheme changed: YES
   - If performance impact significant: YES
   - Otherwise: NO
5. **Update api.md?** Auto-generated - run script if system supports

### 5.2 Update Frequency Targets

- **architecture.md**: 1-2 times per year (major refactors)
- **context.md**: 4-8 times per year (code changes, new features)
- **technical_doc.md**: 2-4 times per year (algorithm updates)
- **api.md**: Per release or monthly (auto-generated)

---

## 6. Quality Checklist

### 6.1 Before Committing architecture.md

- [ ] All major modules listed (≥95% of codebase)
- [ ] Dependency graph has no cycles
- [ ] Data structures match actual arrays in code
- [ ] Design patterns are real and used in code
- [ ] Execution flow computable from pseudocode
- [ ] No prose paragraphs longer than 3 lines
- [ ] Fits on <5 printed pages
- [ ] AI can reason about structure from content

### 6.2 Before Committing context.md

- [ ] Every module in architecture.md has a section here
- [ ] Every subroutine has: Purpose, Inputs, Outputs, Dependencies
- [ ] Real Fortran signatures match actual code
- [ ] Examples compile and run (check syntax)
- [ ] No physics derivations (save for technical_doc.md)
- [ ] Consistent formatting across all module sections
- [ ] AI can understand how to call each function
- [ ] Each section fits on ≤1 page

### 6.3 Before Committing technical_doc.md

- [ ] All mathematical equations have accompanying text
- [ ] References cited for non-obvious results
- [ ] Performance claims backed by section with numbers
- [ ] Examples include validation/test description
- [ ] Convergence rates shown with error plots or tables
- [ ] GPU/MPI implementation justified with performance data
- [ ] Debugging section includes ≥3 real failure modes
- [ ] Humans can implement code from math description

---

## 7. Automation & Tools

### 7.1 Suggested Scripts

Create (optional, not required):
- `validate_architecture.py`: Check module dependencies exist in code
- `validate_context.py`: Verify function signatures match actual files
- `markdown_lint.py`: Check heading depth, line length, code block syntax

### 7.2 GitHub Integration

- Add to `.gitignore`:
  ```
  # Auto-generated only
  /docs/api.md (if auto-gen enabled)
  ```
- Pre-commit hook: Check markdown files follow these rules
- Branch protection: Require doc updates with code changes (optional)

---

## 8. Example: Consistent Documentation Set

**Scenario**: Adding a new scheme `calc_wenofv_kernel_internal.f90`

### architecture.md Addition:
```
Tier 3: [calc_wenofv_kernel_internal] (3D_solver/src/calc_wenofv_kernel_internal.f90) - 5th-order WENO finite-volume
  Dependencies: [mod_constant, calc_weno_3d]
  Kernels: calc_wenofv_x_in, calc_wenofv_y_in, calc_wenofv_z_in
```

### context.md Addition:
```
### calc_wenofv_kernel_internal.f90

**Purpose**: Implement WENO finite-volume scheme for shock-capturing with high accuracy

**Key Functions/Subroutines**:
- calc_wenofv_x_in(nx, ny, nz, Q, sensor, E): x-direction WENO-FV kernel
- calc_wenofv_y_in(nx, ny, nz, Q, sensor, F): y-direction WENO-FV kernel  
- calc_wenofv_z_in(nx, ny, nz, Q, sensor, G): z-direction WENO-FV kernel

**Dependencies**: mod_constant (γ, R, etc.), calc_weno_3d (flux functions)

**Used By**: calc_flux_base.f90 (dispatcher)

**Implementation Notes**:
- WENO-5 (5-point stencil, 5th order on smooth regions)
- Auto-detects shocks and reduces order locally
- Shared memory: stencil data + WENO weights buffer
```

### technical_doc.md Addition:
```
### 3.3.X WENO Finite-Volume Scheme

**Motivation**: KEEP/SLAU offer accuracy or stability; WENO combines both

**Mathematical Basis**: 
Weighted essentially non-oscillatory reconstruction using multiple 3-point substencils...
[equations here]

**Implementation**:
GPU kernel computes adaptive weights per cell, then reconstructs left/right edge values...

**Performance**:
- RTX 4090: ~150 GFLOP/s (more work than KEEP/SLAU)
- Requires 5-point stencil (guard zones must be ≥2 per side)
- CFL ≤ 0.8 recommended (like Roe, less stable than KEEP)

**When to Use**: Hypersonic with strong shocks; slower but more robust
```

---

## 9. Summary of Rules

| Rule | Who Enforces | When Checked |
|------|-------------|--------------|
| architecture.md ≤ 5 pages | Author | Before commit |
| context.md consistent format | Author/Reviewer | Code review |
| technical_doc.md equations cited | Author | Before commit |
| No duplication between files | Author | Manual check |
| Cross-references valid | Tools (lint) | Pre-commit hook |
| AI can parse architecture.md | Copilot test | Before merge |
| Human understands technical_doc.md | Expert review | Code review |

---

**Document Version**: 1.0 (April 5, 2026)  
**Last Updated**: Comprehensive markdown rules established
