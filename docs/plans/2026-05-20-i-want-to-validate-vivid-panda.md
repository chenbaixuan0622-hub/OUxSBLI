# CNN Wall Model Validation Pipeline

## Context

The existing `train_cnn_wall_model.py` generates **synthetic** training data from a Python
NEQWM reference and exports a TorchScript model. The CNN wall model is also embedded in
the Fortran solver via LibTorch (`calc_wall_model_cnn.f90` + `cnn_wall_model.cpp`).

The goal is to create a **self-contained Python validation pipeline** that:
1. Runs the TBL case via the existing Python `Case` API (existing `ouxsbli/` package)
2. Extracts wall model inputs from VTK snapshots and computes NEQWM labels in Python
3. Trains the CNN on real CFD-extracted data (replacing/augmenting synthetic data)
4. Validates CNN predictions against the NEQWM reference with quantitative metrics and plots

This decouples CNN training/validation from the Fortran solver — the solver provides
only flow snapshots; the CNN pipeline is entirely Python.

---

## Critical Files

| File | Action |
|------|--------|
| `run_tbl.py` | **Create** — runs TBL via Python `Case` API |
| `train_cnn_wall_model.py` | **Modify** — add CFD data loading path |
| `validate_cnn_wall_model.py` | **Create** — CNN vs NEQWM comparison |
| `ouxsbli/case.py` | Read-only reference (no changes) |
| `ouxsbli/tests/utils/vtk_reader.py` | Read-only reference (no changes) |
| `3D_solver/TBL/data/` | Existing snapshots `Q00000.vtr`–`Q00010.vtr` |

---

## Step-by-Step Plan

### Step 1 — `run_tbl.py` (new)

Thin wrapper around the existing `Case` API. Demonstrates how to launch a TBL run
programmatically; existing data in `3D_solver/TBL/data/` can be used directly if
a fresh run is not needed.

```
from ouxsbli import Case

Case(
    source="3D_solver/TBL",
    workdir="./output/tbl_run",
    visc="ns",      # Navier-Stokes (integer(4)) — full-resolution for data collection
    scheme="slau",  # SLAU flux
    accuracy=2,
).build().run(nranks=2)
```

Prints the output data directory on completion.

---

### Step 2 — Shared data-extraction helper (added to `train_cnn_wall_model.py`)

New function `load_cfd_data(data_dir, jm_idx=2, offset=1)`:

1. Glob all `Q*.vtr` files in `data_dir`, sorted by numeric suffix.
2. For each snapshot call `vtk_reader.getQ(path, nx, ny, nz)` and
   `vtk_reader.getGrid(path)`.
3. Extract at matching-height index `jm_idx` (Python 0-based; `jm_wmles=3` → `j=2`):
   - `u_jm = data['u'][:, jm_idx, :]` — shape `(Nz, Nx)`
   - `w_jm = data['w'][:, jm_idx, :]`
   - `T_jm = data['p'][:, jm_idx, :] / (data['rho'][:, jm_idx, :] * R_GAS)`
   - `p_jm = data['p'][:, jm_idx, :]`
   - `y_jm` read from `getGrid()` Y-coordinates at index `jm_idx`
4. Apply `neqwm_scalar()` element-wise (already implemented) → labels
   `(tau_wx[i,k], tau_wz[i,k])`.
5. Flatten spatial dims: X shape `(N_pts, 4)`, Y shape `(N_pts, 2)`.
6. Normalise identically to synthetic path:
   - X: `[u/u0, w/u0, T/T0, p/p0]`
   - Y: `[tau_wx, tau_wz] / tau_scale`  where `tau_scale = p0/(R_GAS*T0)*u0**2`
7. Return `X, Y` as `float64` numpy arrays.

**Note on grid size:** TBL has `nx=257, nz=33`; with `offset=1` (2nd order) the
interior region used by the wall model is `Ni=255, Nk=31`.

---

### Step 3 — Modify `train_cnn_wall_model.py`

Add CLI arguments (via `argparse`):

| Argument | Default | Purpose |
|----------|---------|---------|
| `--data-source {synthetic,cfd}` | `synthetic` | Choose data origin |
| `--data-dir PATH` | — | VTK directory (required when `cfd`) |
| `--test-split FLOAT` | `0.2` | Fraction held out for evaluation |
| `--epochs INT` | `100` | Training epochs |
| `--lr FLOAT` | `1e-3` | Adam learning rate |
| `--model-out PATH` | `cnn_wall_model.pt` | Output TorchScript path |

When `--data-source cfd`:
- Call `load_cfd_data(data_dir)` → get full X, Y
- Split indices: last `test_split` fraction → test set; remainder → train set
- Reshape train set into `(N_patches, 4, grid_ni, grid_nk)` tiles for CNN input
  (same reshape logic as the existing synthetic path; grid_ni=16, grid_nk=16 default)

After training, evaluate on the **held-out test set** and print:
- MAE for `tau_wx` and `tau_wz` (in physical units and normalised)
- RMSE for `tau_wx` and `tau_wz`
- R² for each component

Export TorchScript model as before.

---

### Step 4 — `validate_cnn_wall_model.py` (new)

Standalone script; takes `--data-dir` and `--model-path` as arguments.

**Algorithm:**

1. Load VTK snapshots via `load_cfd_data()` (reused from step 2).
2. Reshape X to `(N_patches, 4, Ni, Nk)` matching the TBL interior grid `(255, 31)`.
3. Load TorchScript model: `torch.jit.load(model_path).eval()`.
4. Run CNN forward pass on each patch → `Y_cnn` (normalised).
5. Re-apply `neqwm_scalar()` on same inputs → `Y_neqwm` (as ground truth).
6. Denormalise both sets → physical units `[Pa]`.
7. **Metrics** (printed to stdout):
   - MAE, RMSE, max-abs-error for `tau_wx` and `tau_wz`
   - R² and Pearson correlation for each component
8. **Plots** (saved as PNG):
   - `scatter_tau_wx.png` / `scatter_tau_wz.png`: CNN vs NEQWM scatter with
     1-to-1 reference line and R² annotation
   - `wall_stress_map.png`: side-by-side 2D maps (NEQWM | CNN | relative error)
     for the first test snapshot
   - `error_histogram.png`: distribution of relative errors for both components

---

## Data & Architecture Notes

- **Existing data** `3D_solver/TBL/data/Q00000.vtr`–`Q00010.vtr` (11 snapshots)
  gives ~86,500 wall-model input points — sufficient for a meaningful offline test.
- **CNN architecture** (unchanged): Conv2d 4→16→32→2, 3×3 kernels, padding=1;
  5×5 effective receptive field; input normalised by `(u0, u0, T0, p0)`.
- **`y_jm`**: the physical height of `jm_wmles=3` (j=2 in Python 0-based) is read
  from the VTK Y-coordinates so the NEQWM bisection uses the correct grid spacing.

---

## Verification

```bash
# 1. (Optional) Run fresh TBL simulation
python run_tbl.py

# 2. Train CNN on CFD data
python train_cnn_wall_model.py \
    --data-source cfd \
    --data-dir 3D_solver/TBL/data \
    --test-split 0.2 \
    --epochs 100

# 3. Validate
python validate_cnn_wall_model.py \
    --data-dir 3D_solver/TBL/data \
    --model-path cnn_wall_model.pt

# Expected: R² > 0.95, relative MAE < 5% on held-out snapshots
```

No Fortran recompilation required; the pipeline is self-contained Python.
