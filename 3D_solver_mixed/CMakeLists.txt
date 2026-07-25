cmake_minimum_required(VERSION 3.18)  # list(SORT ... COMPARE NATURAL), used below for the CUDA stub search

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER mpif90)
endif()

if(NOT DEFINED CASE_GPU_CC)
  set(CASE_GPU_CC "native")
endif()

if(NOT DEFINED CASE_MAIN_SOURCE)
  set(CASE_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/main.f90")
endif()

set(SUBROUTINES "name:calc_tau_straight_mp,calc_tau_cross_mp,interp6_mp,diff6_mp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mpreprocess -Mchkptr -Minline=reshape,${SUBROUTINES} -Minfo=accel,inline -gpu=ptxinfo,rdc,lto,cc${CASE_GPU_CC},lineinfo")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Locate libcuda.so for green_ctx_bindings_mixed.f90's driver-API calls
# (mod_streams_mixed.f90's Green Context stream setup). These are compiled
# into every case unconditionally (see the _BASE comment in
# add_ouxsbli_case() below), so this is needed for every build, not just
# OVERLAP_GC=True ones. Mirrors the same stub-library search
# 3D_solver_fltflt/CMakeLists.txt uses -- the real driver library
# isn't reliably on the default linker search path on every system (e.g.
# WSL keeps it under /usr/lib/wsl/lib, HPC login/build nodes often have no
# GPU driver installed at all), so link against the NVHPC SDK's stub
# libcuda.so at build time and let the dynamic loader resolve the real
# driver .so at runtime on whichever node actually has a GPU.
if(DEFINED ENV{NVHPC_ROOT})
  file(GLOB _CUDA_STUB_DIRS "$ENV{NVHPC_ROOT}/cuda/*/targets/*/lib/stubs")
  # Green Contexts need CUDA 12.4+; prefer the newest bundled stub rather
  # than whichever the filesystem glob happens to list first.
  list(SORT _CUDA_STUB_DIRS COMPARE NATURAL ORDER DESCENDING)
endif()
find_library(CUDA_DRIVER_LIBRARY
  NAMES cuda
  HINTS ${_CUDA_STUB_DIRS}
  PATHS /usr/lib/wsl/lib /usr/lib64 /usr/lib/x86_64-linux-gnu)

# Macro: call from each case's CMakeLists.txt after setting CASE_EXTRA_SOURCES.
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory.
macro(add_ouxsbli_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src" "-I${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src" "-DMIXED_FORK=True")

  set(_PARENT_FYPP "mod_constant" "calc_muscl")
  set(_SHARED_FYPP
      "load_smem_visc_cent_mp"
      "calc_visc_high_internal_mp"
      "calc_slau_kernel_internal_acc"
      "calc_flux_base_mixed"
      "calc_time_dev_mixed")
  set(_BASELINE_FYPP
      "preprocess"
      "calc_slau_kernel"
      "calc_slau_kernel_internal")

  set(_GEN "")
  foreach(_f ${_PARENT_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  foreach(_f ${_SHARED_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_visc_cent_mp.f90.fypp" "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  foreach(_f ${_BASELINE_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../../3D_solver/src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  # green_ctx_bindings_mixed.f90/mod_streams_mixed.f90 (and the CUDA driver
  # lib below) are ALWAYS compiled/linked, regardless of the case's
  # config.fypp OVERLAP_GC value. OVERLAP_GC only gates fypp `#:if` content
  # inside calc_flux_base_mixed.f90.fypp/calc_time_dev_mixed.f90.fypp (whether
  # calc_conv_visc_paired_mixed/init_streams are compiled in) -- it must NOT
  # also gate which source files CMake puts in the target, because that
  # decision is only re-evaluated at `cmake -B` configure time. Editing
  # config.fypp and re-running `cmake --build` alone (the normal workflow for
  # every other flag here -- SCHEME/RESCALE/COMMZ) does NOT rerun that
  # configure-time parse, so a CMake-level OVERLAP_GC gate would silently
  # leave mod_streams_mixed.f90 out of the build after flipping the flag,
  # producing "Unable to open MODULE file mod_streams_mixed.mod" the moment
  # the freshly fypp-regenerated calc_flux_base_mixed.f90 tries to `use` it.
  # Compiling these two small files unconditionally costs nothing and keeps
  # OVERLAP_GC consistent with how every other config.fypp flag behaves.
  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/cpu_gpu_mpi.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_coordinate.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_physical_quantities_mixed.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/print.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/calc_hybrid.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/calc_steps.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/set_bc_common.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/green_ctx_bindings_mixed.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/mod_streams_mixed.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_executable(${_CASE} ${_BASE} ${CASE_MAIN_SOURCE} ${CASE_EXTRA_SOURCES} ${_GEN})
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  # CUDA driver lib (libcuda.so / stub), needed by green_ctx_bindings_mixed.f90's
  # Green Context calls -- distinct from the runtime API libs nvfortran's
  # -cuda flag already links. Linked unconditionally for the same reason the
  # source files above are always compiled -- see the comment above _BASE.
  if(CUDA_DRIVER_LIBRARY)
    target_link_libraries(${_CASE} PRIVATE ${CUDA_DRIVER_LIBRARY})
  else()
    target_link_options(${_CASE} PRIVATE -lcuda)
  endif()
  if (DEFINED CASE_LINK_OPTIONS)
    target_link_options(${_CASE} PRIVATE ${CASE_LINK_OPTIONS})
  endif()
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
