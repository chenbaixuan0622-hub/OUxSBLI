# Plan: Shared CMakeLists for 3D solver cases

## Context

`3D_solver/ETGV/CMakeLists.txt` is the only working CMake build. It contains ~126 lines of
mostly-common logic. The user wants the common parts extracted into `3D_solver/CMakeLists.txt`
and slim per-case files prepared for TBL and SBLI.

**Key findings:**
- TBL and SBLI need two extra hand-written sources ETGV doesn't use:
  `3D_solver/src/set_init_common.f90` and `src/set_compressible_bl.f90`
- SBLI additionally needs `3D_solver/src/set_bc_tbl_sbli.f90`
- SBLI's Makefile references `sbli.o` but `sbli.f90` does not exist anywhere;
  `main.f90` (`3D_solver/src/main.f90`) is used instead
- `set_compressible_bl.f90` lives in top-level `src/` (i.e. `../../src/` from a case dir)
- `set_init_common.f90` lives in `3D_solver/src/` (i.e. `../src/` from a case dir)

---

## Structure

```
3D_solver/
├── CMakeLists.txt          ← NEW: project + macro + add_subdirectory
├── ETGV/
│   └── CMakeLists.txt      ← UPDATED: slim (set extras + call macro)
├── TBL/
│   └── CMakeLists.txt      ← NEW
└── SBLI/
    └── CMakeLists.txt      ← NEW
```

CMake `macro()` is used (not `function()`): in a macro, `CMAKE_CURRENT_SOURCE_DIR` and
`CMAKE_CURRENT_BINARY_DIR` resolve at call time to the calling subdirectory — exactly
what we need so relative paths like `../src/` work correctly per case.

---

## `3D_solver/CMakeLists.txt` (NEW)

```cmake
cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER mpif90)
endif()

project(OUxSBLI LANGUAGES Fortran)

set(SUBROUTINES "name:calc_tau_straight,calc_tau_cross,flux4")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mpreprocess -Mchkptr \
-Minline=reshape,${SUBROUTINES} -Minfo=accel,inline -gpu=ptxinfo,rdc,lto -lnvhpcwrapnvtx")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Macro: call from each case's CMakeLists.txt after setting CASE_EXTRA_SOURCES.
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory.
macro(add_ousbli_case)
  set(FYPP_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}")

  set(_PARENT_FYPP "mod_constant")
  set(_SHARED_FYPP
      "calc_visc2" "calc_visc4" "calc_visc4_internal"
      "calc_keep_kernel" "calc_keep_kernel_internal"
      "calc_slau_kernel" "calc_slau_kernel_internal"
      "calc_roe_kernel" "calc_roe_kernel_internal"
      "calc_hybrid_kernel" "calc_hybrid_kernel_internal"
      "calc_flux_base" "calc_time_dev" "preprocess")

  set(_GEN "")
  foreach(_f ${_PARENT_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${CMAKE_CURRENT_SOURCE_DIR}/config.fypp"
      COMMENT "fypp ../../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()
  foreach(_f ${_SHARED_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${CMAKE_CURRENT_SOURCE_DIR}/config.fypp"
      COMMENT "fypp ../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  set(_BASE
      "${CMAKE_CURRENT_SOURCE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/cpu_gpu_mpi.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_coordinate.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc2.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc4.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_les.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_hybrid.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_muscl.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_steps.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_rescale.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_para.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_bc_common.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/print.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/main.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_executable(${_CASE} ${_BASE} ${CASE_EXTRA_SOURCES} ${_GEN})
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()

add_subdirectory(ETGV)
add_subdirectory(TBL)
add_subdirectory(SBLI)
```

---

## `3D_solver/ETGV/CMakeLists.txt` (UPDATED — replace current content)

```cmake
set(CASE_EXTRA_SOURCES "")
add_ousbli_case()
```

---

## `3D_solver/TBL/CMakeLists.txt` (NEW)

```cmake
set(CASE_EXTRA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_compressible_bl.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_init_common.f90")
add_ousbli_case()
```

---

## `3D_solver/SBLI/CMakeLists.txt` (NEW)

```cmake
# sbli.f90 does not exist; main.f90 is used as the entry point
set(CASE_EXTRA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_bc_tbl_sbli.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_compressible_bl.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_init_common.f90")
add_ousbli_case()
```

---

## Verification

```bash
cd /home/jhatayama/ouxsbli/fypp/OUxSBLI
source ~/htymenv_new/bin/activate
module load nvhpc/25.5

cmake -S 3D_solver -B 3D_solver/build
cmake --build 3D_solver/build --target ETGV 2>&1 | tail -20
cmake --build 3D_solver/build --target TBL   2>&1 | tail -20
cmake --build 3D_solver/build --target SBLI  2>&1 | tail -20
```

Executables land at:
- `3D_solver/build/ETGV/a.out`
- `3D_solver/build/TBL/a.out`
- `3D_solver/build/SBLI/a.out`
