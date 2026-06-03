# CNN-based Wall Model via LibTorch

## Context
The NEQWM wall model (`calc_wall_model.f90`, `kind(id_wmles)==4`) is now in place. The goal is to add a CNN variant (`kind(id_wmles)==8`) that replaces the bisection+ODE solver with a 2D convolutional neural network over the x-z wall plane. PyTorch is used both for training (Python) and for runtime inference (LibTorch C++ API). Because the Fortran dispatch already uses a kind-flag for NEQWM, adding kind=8 requires only small, targeted changes. Non-wall-model cases (ETGV, KHI, …) are guarded by a C preprocessor flag so they need no changes to their Makefiles.

---

## Architecture

```
Fortran (GPU)          C++ / LibTorch (CUDA)
─────────────────      ─────────────────────────────────────────
calc_and_apply_wm_cnn  cnn_wm_infer_c
  c_devloc(Q) ──────►  from_blob( Q_dev, strides=Fortran )
  c_devloc(T) ──────►  extract [u,w,T,p] slice at jm_wmles
  c_devloc(F) ◄──────  forward pass → tau_wx, tau_wz
                        write to tau arrays + F(2/4,:,1,:)
```

CNN shape: input `(1, 4, Ni, Nk)` → Conv2d stack → output `(1, 2, Ni, Nk)`.  
Inference runs entirely on GPU; no host-device transfers.

---

## Files to Create

### `3D_solver/src/calc_wall_model_cnn.f90`
New standalone Fortran module `calc_wall_model_cnn`.
- `use iso_c_binding`, `use cudafor`
- C interfaces: `cnn_wm_init_c`, `cnn_wm_infer_c` (bind(C))
- Private device arrays: `tau_wx_cnn(nx,nz)`, `tau_wz_cnn(nx,nz)` 
- `init_wall_model_cnn(nx, nz)` — allocates device arrays, calls `cnn_wm_init_c(trim(cnn_model_path)//c_null_char, nx, nz)`
- `calc_and_apply_wm_cnn(nx, ny, nz, offset, inv_dy, Q, T, mu, F)` — identical signature to `calc_and_apply_wm`; uses `c_devloc` to get raw device pointers and calls `cnn_wm_infer_c`

### `3D_solver/src/cnn_wall_model.cpp`
LibTorch C wrapper (~80 lines).
- Static `torch::jit::script::Module g_model` loaded once
- `cnn_wm_init_c`: loads TorchScript `.pt`, moves to CUDA, sets `eval()`
- `cnn_wm_infer_c`: wraps `Q_dev` and `T_dev` via `torch::from_blob` with Fortran-order strides `{1, nx, nx*5, nx*5*ny}` and `{1, nx, nx*ny}`; extracts interior slice at `jm-1`; stacks `[u, w, T, p]` into `(1, 4, Ni, Nk)` float32 tensor; calls `g_model.forward()`; writes `tau_wx/tau_wz` back and updates `F(comp=1/3, x_face, j_face=0, z_face)` using correct Fortran-order strides `{1, 5, 5*(nx-2), 5*(nx-2)*(ny-1)}`

### `train_cnn_wall_model.py`
PyTorch training script (repository root or `scripts/`).
- `WallModelCNN`: `Conv2d(4→16, 3×3, pad=1) → ReLU → Conv2d(16→32, 3×3, pad=1) → ReLU → Conv2d(32→2, 1×1)`
- Data generation: Python re-implementation of the NEQWM bisection formula — sample random `(u_e, w_e, T_e, p_e, y_jm)` tuples, compute `(tau_wx, tau_wz)` analytically; no solver run required
- Normalise inputs by reference values `(u0, u0, T0, p0)` and output by `rho_ref * u0^2`
- `torch.jit.script` + `model.save("cnn_wall_model.pt")`

---

## Files to Modify

### `3D_solver/src/calc_flux_base.f90`
Two small additions (guarded by `#ifdef USE_CNN_WM`):

1. After the existing `use calc_wall_model` line (inside `calc_EFG_LES`):
```fortran
#ifdef USE_CNN_WM
    use calc_wall_model_cnn, only : calc_and_apply_wm_cnn
#endif
```

2. After the existing `kind(id_wmles)==4` block (≈line 338):
```fortran
#ifdef USE_CNN_WM
    elseif (kind(id_wmles) == 8 .and. id_bc_y) then
      block
        integer :: wm_offset
        if (kind(id_accuracy) == 8) then
          wm_offset = 3
        elseif (kind(id_accuracy) == 4) then
          wm_offset = 2
        else
          wm_offset = 1
        endif
        call calc_and_apply_wm_cnn(nx, ny, nz, wm_offset, inv_dy, Q, T, mu, F)
      end block
#endif
```

### `3D_solver/src/preprocess.f90`
Inside `allocate_device_mem`, after the existing `kind(id_wmles)==4` block:
```fortran
#ifdef USE_CNN_WM
    if (kind(id_wmles) == 8) call init_wall_model_cnn(nx, nz)
#endif
```
Also add `#ifdef USE_CNN_WM / use calc_wall_model_cnn, only: init_wall_model_cnn / #endif` to the module-level use section of `allocate_device_mem`.

### `3D_solver/TBL/mod_globals.f90` and `3D_solver/SBLI/mod_globals.f90`
Two additions:
1. Update the `id_wmles` comment block to include `! kind8 CNN`
2. Add model path parameter after `jm_wmles`:
```fortran
character(len=256), parameter :: cnn_model_path = './cnn_wall_model.pt'
```

### `3D_solver/TBL/Makefile` and `3D_solver/SBLI/Makefile`
Add near the top (after `FC = mpif90`):
```makefile
LIBTORCH_DIR  = $(shell python3 -c "import torch; print(torch.__file__[:-12])" 2>/dev/null)
LIBTORCH_INC  = -I$(LIBTORCH_DIR)/include \
                -I$(LIBTORCH_DIR)/include/torch/csrc/api/include
LIBTORCH_LIBS = -L$(LIBTORCH_DIR)/lib \
                -ltorch -ltorch_cpu -ltorch_cuda -lc10 -lc10_cuda \
                -Wl,-rpath,$(LIBTORCH_DIR)/lib
```
Add to `OUTPUT_OPTIN`: `-DUSE_CNN_WM`

Add compilation rule:
```makefile
cnn_wall_model.o: ../../src/cnn_wall_model.cpp
	g++ -O2 -fPIC -std=c++17 $(LIBTORCH_INC) -c $<
```

Add `calc_wall_model_cnn.o` to `OBJ` (after `calc_wall_model.o`).

Change `a.out` link line to include `cnn_wall_model.o` and `$(LIBTORCH_LIBS)`:
```makefile
a.out: $(OBJ) cnn_wall_model.o
	$(FC) -cuda -acc -fast -Minfo -gpu=rdc,lto -lnvhpcwrapnvtx -I./ $^ $(LIBTORCH_LIBS)
```

Add dependency:
```makefile
calc_wall_model_cnn.o: mod_globals.mod
```

---

## Enabling CNN Wall Model

In `3D_solver/TBL/mod_globals.f90`:
```fortran
integer(kind=8), parameter :: id_wmles = 0   ! kind=8 → CNN wall model
```
Place `cnn_wall_model.pt` in the case directory before running.

---

## Verification

1. Train model and export: `python3 train_cnn_wall_model.py` → `cnn_wall_model.pt`
2. Copy `.pt` to `3D_solver/TBL/`
3. Set `id_wmles = integer(8)` and `id_visc = integer(8)` (LES) in TBL `mod_globals.f90`
4. `cd 3D_solver/TBL && make clean && make` — should compile without errors
5. `bash calc.sh` — simulation runs; check that wall-stress values are non-zero and physically reasonable
6. Compare VTK output wall-plane fields between NEQWM (kind=4) and CNN (kind=8) runs
7. Non-CNN cases: `cd 3D_solver/ETGV && make clean && make` — should compile unchanged (no `-DUSE_CNN_WM`, no `cnn_wall_model.o`)
