# Plan: Load NPZ files to train CNN wall model

## Context

The CNN wall model training script (`cnnwm/train_cnn_wall_model.py`) currently supports two data sources: synthetic (random NEQWM samples) and CFD (VTK `.vtr` snapshots). The user has prepared two large NPZ archives:

- `cnnwm/Qwall.npz` — primitive flow variables `[ρ, u, v, w, p]` at the wall-adjacent cell (j=1), shape `[5, Time, Nz, Nx]`. Used to compute **ground-truth wall shear stress labels** via viscous approximation: `τ_wx = μ(T_wall) × u_j1 / Δy`, `τ_wz = μ(T_wall) × w_j1 / Δy`.
- `cnnwm/Q10pct.npz` — primitive flow variables `[ρ, u, v, w, p]` at y = 10 % of the boundary layer thickness (the CNN matching height), shape `[5, Time, Nz, Nx]`. Used as **CNN inputs** `[u/u0, w/u0, T/T0, p/p0]`.

Files are ~17.3 GB each. They must not be loaded simultaneously. Only the first `n_time` time steps (default 1000) are used.

A new `cnnwm/config.yaml` file will hold all run parameters; the training script reads it at startup, with any CLI flags overriding individual values.

---

## Critical files

| File | Change |
|------|--------|
| `cnnwm/config.yaml` | **Create new** — all parameters live here |
| `cnnwm/train_cnn_wall_model.py` | Add config loading, `load_npz_data()`, and `npz` data-source branch |

---

## 1 — `cnnwm/config.yaml` (new file)

```yaml
# ── Data source ───────────────────────────────────────────
data_source: npz          # synthetic | cfd | npz
data_dir: null            # VTK snapshot dir (cfd mode only)
q10pct: cnnwm/Q10pct.npz  # matching-height NPZ (npz mode)
qwall:  cnnwm/Qwall.npz   # wall-adjacent NPZ   (npz mode)
n_time: 1000              # time steps to load  (npz mode)
dy_wall: 1.0e-4           # Δy from wall to first interior cell [m] (npz mode)

# ── Reference values (match TBL mod_globals) ──────────────
u0: 509.7     # m/s
T0: 171.31    # K
p0: 14924.0   # Pa

# ── Training ──────────────────────────────────────────────
n_epochs:   100
lr:         1.0e-3
batch:      64
grid_ni:    16
grid_nk:    16
test_split: 0.2
n_samples:  200000   # synthetic mode only

# ── Output ────────────────────────────────────────────────
model_out: cnn_wall_model.pt
```

---

## 2 — Changes to `cnnwm/train_cnn_wall_model.py`

### 2a — Add `PyYAML` import and `GAMMA` constant

```python
import yaml        # add to imports
GAMMA = 1.4        # ratio of specific heats, near existing physical constants
```

### 2b — Add `load_config()` helper (top of file, before `main`)

```python
def load_config(path: str) -> dict:
    """Load config.yaml; missing keys fall back to argparse defaults."""
    with open(path) as f:
        return yaml.safe_load(f)
```

### 2c — Add `load_npz_data()` function (after `load_cfd_data`)

```python
def load_npz_data(q10pct_path, qwall_path, u0, T0, p0, dy_wall, n_time=1000):
```

**Step 1 — load Q10pct (inputs), then release:**
- `np.load(q10pct_path, mmap_mode='r')` — only the requested slice is paged in
- Slice key: try `'arr_0'` first, fall back to first key in file
- Slice: `Q10[:, :n_time, :, :]` → shape `[5, T, Nz, Nx]`
- Primitives: `rho=Q[0]`, `u=Q[1]`, `w=Q[3]`, `p=Q[4]`, `T = p / (rho * R_GAS)`
- Build `X` shape `[T*Nz*Nx, 4]` = `[u/u0, w/u0, T/T0, p/p0]`, `float32`
- `del Q10` to free mmap before opening Qwall

**Step 2 — load Qwall (labels):**
- Same mmap load; slice `Qw[:, :n_time, :, :]`
- `rho_w=Qw[0]`, `u_w=Qw[1]`, `w_w=Qw[3]`, `p_w=Qw[4]`
- `T_w = p_w / (rho_w * R_GAS)`
- `mu_w = sutherland(T_w)` — reuses existing `sutherland()` function
- `tau_wx = mu_w * u_w / dy_wall`  (all numpy, no Python loop)
- `tau_wz = mu_w * w_w / dy_wall`
- `tau_scale = (p0 / (R_GAS * T0)) * u0**2`
- Build `Y` shape `[N, 2]` = `[tau_wx/tau_scale, tau_wz/tau_scale]`, `float32`

Print shape and sample count before returning.

### 2d — Update `argparse` in `main()`

- `--config` argument (default `cnnwm/config.yaml`) — path to YAML
- Load YAML first; CLI arguments override individual keys if provided (use `argparse` defaults to detect "not set by user")
- Extend `--data-source` choices to include `"npz"`
- Add `--q10pct`, `--qwall`, `--dy-wall`, `--n-time` arguments (all sourced from config by default)

### 2e — Update dispatch in `main()`

```python
elif args.data_source == "npz":
    X, Y = load_npz_data(args.q10pct, args.qwall,
                          args.u0, args.T0, args.p0,
                          args.dy_wall, args.n_time)
```

The rest of the pipeline (train/test split → `make_patch_dataset` → `train()` → export) is **unchanged**.

---

## Key design choices

- **`config.yaml` is authoritative**: all parameters live there; CLI overrides are optional. This avoids long command lines and makes runs reproducible.
- **`mmap_mode='r'`**: numpy reads only the sliced region from the zip archive, keeping peak RAM well below the 17 GB file size.
- **Sequential loading**: `X` is fully built and the mmap closed before Qwall is opened, so the two 17 GB arrays are never simultaneously in RAM.
- **Pure numpy labels**: `τ = μ × u / Δy` is fully vectorized — no per-point Python loop — and uses actual LES/DNS velocity for higher-quality training targets than NEQWM.

---

## Verification

```bash
cd /home/jhatayama/ouxsbli/wall/OUxSBLI
# Edit config.yaml: set n_time: 10, n_epochs: 2
python cnnwm/train_cnn_wall_model.py --config cnnwm/config.yaml
# Expect: shape printout, 2 epoch losses, test metrics, model saved
```
