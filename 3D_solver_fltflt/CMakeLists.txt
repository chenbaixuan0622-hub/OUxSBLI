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

set(SUBROUTINES "name:calc_tau_straight,calc_tau_cross,flux4,fltflt_two_sum,fltflt_fast_two_sum,fltflt_two_prod_fma,fltflt_add_same_sign,fltflt_fma,fltflt_fma_approx,fltflt_abs,fltflt_sqrt,fltflt_sqrt_fast,fltflt_norm3d,fltflt_dot2,fltflt_dot3,fltflt_add3,fltflt_add4,fltflt_add5,fltflt_add_mul,fltflt_sub_mul,fltflt_add_add_mul,fltflt_add_abs,fltflt_sub_abs,fltflt_square,fltflt_recip,fltflt_shfl_down,fltflt_shfl_xor,add_ff_ff,add_ff_r4,add_r4_ff,add_ff_r8,add_r8_ff,neg_ff,sub_ff_ff,sub_ff_r4,sub_r4_ff,sub_ff_r8,sub_r8_ff,mul_ff_ff,mul_ff_r4,mul_r4_ff,mul_ff_r8,mul_r8_ff,div_ff_ff,div_ff_r4,div_r4_ff,div_ff_r8,div_r8_ff,eq_ff_ff,eq_ff_r4,eq_r4_ff,ne_ff_ff,ne_ff_r4,ne_r4_ff,lt_ff_ff,lt_ff_r4,lt_r4_ff,gt_ff_ff,gt_ff_r4,gt_r4_ff,le_ff_ff,le_ff_r4,le_r4_ff,ge_ff_ff,ge_ff_r4,ge_r4_ff,fltflt_min_ff_ff,fltflt_min_ff_r4,fltflt_min_r4_ff,fltflt_max_ff_ff,fltflt_max_ff_r4,fltflt_max_r4_ff")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mpreprocess -Mchkptr -Minline=reshape,${SUBROUTINES} -Minfo=accel,inline -gpu=ptxinfo,rdc,lto,cc${CASE_GPU_CC},lineinfo -lnvhpcwrapnvtx")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Locate libcuda.so for green_ctx_bindings_ff.f90's driver-API calls
# (mod_streams_ff.f90's Green Context stream setup). The real driver
# library isn't reliably on the default linker search path on every
# system -- e.g. WSL keeps it under /usr/lib/wsl/lib, and HPC cluster
# login/build nodes often have no real GPU driver installed at all. The
# NVHPC SDK ships its own stub libcuda.so precisely for this: link
# against the stub at build time (it only needs to satisfy the linker),
# and the real driver .so is resolved at runtime by the dynamic loader
# on whichever node actually has a GPU. Falls back to a bare "-lcuda"
# link flag if nothing is found, preserving the previous behavior on
# systems where the real driver is already on the default search path.
if(DEFINED ENV{NVHPC_ROOT})
  file(GLOB _CUDA_STUB_DIRS "$ENV{NVHPC_ROOT}/cuda/*/targets/*/lib/stubs")
  # An HPC SDK install bundles multiple CUDA toolkit versions side by side
  # (e.g. 11.8 and 12.5); Green Contexts need 12.4+, so the search must
  # prefer the newest stub rather than whichever the filesystem glob
  # happens to list first -- otherwise find_library can silently pick an
  # older stub whose libcuda.so lacks cuGreenCtxCreate & friends, which
  # only fails much later at link time with an obscure undefined-reference
  # error instead of a clear "not found".
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
  set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src" "-I${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src" "-DFLTFLT_FORK=True")

  set(_PARENT_FYPP "mod_constant" "calc_muscl")
  set(_SHARED_FYPP
      "calc_keep_kernel_internal_ff"
      "calc_slau_kernel_internal_ff"
      "calc_hybrid_kernel_ff"
      "calc_hybrid_kernel_internal_ff"
      "calc_flux_base_ff"
      "calc_time_dev_ff")
  set(_BASELINE_FYPP
      "calc_visc2"
      "calc_visc_high"
      "calc_visc_high_internal"
      "load_smem_visc_cent"
      "preprocess")

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
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
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

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/cpu_gpu_mpi.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_coordinate.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/fltflt.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/green_ctx_bindings_ff.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/mod_streams_ff.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_hybrid_ff.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/load_smem_visc2.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/load_smem_visc_me4_base.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/calc_steps.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/set_bc_common.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/print.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_executable(${_CASE} ${_BASE} ${CASE_MAIN_SOURCE} ${CASE_EXTRA_SOURCES} ${_GEN})
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  # CUDA driver lib (libcuda.so / stub), needed by green_ctx_bindings_ff.f90's
  # Green Context calls -- distinct from the runtime API libs nvfortran's
  # -cuda flag already links. See CUDA_DRIVER_LIBRARY search above.
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
