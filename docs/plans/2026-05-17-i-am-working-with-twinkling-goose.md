# Plan: Complete test_os.py (2D oblique shock integration test)

## Context

`ouxsbli/tests/test_os.py` tests the 2D_solver/OS case: an M=2 oblique shock (β=37.2°) that reflects off the bottom wall (β_r=37.7°). The test builds and runs the solver, reads the VTK output, then is supposed to assert that the simulated pre-shock and post-shock states match the Rankine-Hugoniot analytical solution. The assertions are missing entirely (file ends at line 74 with no assertions), and there are several bugs in the shared utility `oblique_shock.py`.

---

## Files to modify

| File | Status |
|------|--------|
| `ouxsbli/tests/utils/oblique_shock.py` | Fix 2 bugs + add 2 missing functions |
| `ouxsbli/tests/test_os.py` | Fix signature + add assertions |

---

## 1. Fix `ouxsbli/tests/utils/oblique_shock.py`

### Bug A — `oblique_shock()` uses `gamma` and `R` as free variables
They are not module-level globals and not in the parameter list → `NameError` at runtime.

**Fix:** add `gamma` and `R` as parameters:
```python
def oblique_shock(M0, p0, T0, beta, gamma, R):
```

### Bug B — `reslected_shock()` uses `Ms` and `Ms2` before defining them
Line 24 references `Ms2`; line 28 references `Ms`. Neither is assigned in this function.

**Fix:** add after line 22 (`rho2, p2, T2 = oblique_shock(...)`):
```python
Ms  = M0 * np.sin(beta)
Ms2 = Ms**2
```
Also update both `oblique_shock(...)` calls inside `reslected_shock` to pass `gamma, R`.

### Missing function A — `beta_from_theta(M, theta_deg, gamma=1.4)`
Required by `post_shock_state` (see below) and documented in CLAUDE.md.
Bisects the θ-β-M relation `tan(θ) = 2·cot(β)·(M²sin²β-1)/(M²(γ+cos2β)+2)`
over `β ∈ [arcsin(1/M), π/2]` to find the weak-shock solution.

### Missing function B — `post_shock_state(M, theta_deg, rho1=1.0, p1=1.0, gamma=1.4)`
Required by `test_corn.py` (already written, already imports this function).
1. Call `beta_from_theta` to get β.
2. Ms = M·sin(β)
3. `p_ratio  = 1 + 2γ(Ms²-1)/(γ+1)`
4. `rho_ratio = (γ+1)Ms² / (2+(γ-1)Ms²)`
5. Return `{"p_ratio": p_ratio, "rho_ratio": rho_ratio}`.

---

## 2. Complete `ouxsbli/tests/test_os.py`

### Fix A — invalid fixture `params`
`params` is not a pytest fixture anywhere in conftest.py.  
Change signature to `tmp_path` (consistent with `test_st.py`, `test_corn.py`).

### Fix B — add `free_stream` import
```python
from .utils.oblique_shock import reslected_shock, free_stream
```

### Fix C — add assertions after line 73

```python
# free-stream (pre-shock) analytical state
u0_fs, p0_fs, T0_fs = free_stream(M0, gamma, R, p_tot, T_tot)
rho0_fs = p0_fs / (R * T0_fs)

# Pre-shock region: left 10 % of domain (i < nx//10), top half (j >= ny//2).
# This is purely upstream of the incident shock for a β≈37° shock in a 5×2 mm domain.
rho_pre = rho[0, ny // 2 :, : nx // 10].astype(float)
u_pre   = u[0,   ny // 2 :, : nx // 10].astype(float)
p_pre   = p[0,   ny // 2 :, : nx // 10].astype(float)

# Post-reflected-shock region: right 20 % of domain (i >= 4*nx//5), bottom quarter (j < ny//4).
# At x = 4Lx/5 = 4 mm the reflected shock (from x_hit≈2.6 mm) is at y≈1.1 mm,
# so j < ny//4 ≈ 32 is well below the reflected shock.
rho_post = rho[0, : ny // 4, 4 * nx // 5 :].astype(float)
u_post   = u[0,   : ny // 4, 4 * nx // 5 :].astype(float)
v_post   = v[0,   : ny // 4, 4 * nx // 5 :].astype(float)
p_post   = p[0,   : ny // 4, 4 * nx // 5 :].astype(float)

# pre-shock assertions (free stream)
assert abs(rho_pre.mean() - rho0_fs) / rho0_fs < PRE_RTOL, (
    f"Pre-shock rho: {rho_pre.mean():.4f} vs {rho0_fs:.4f}")
assert abs(u_pre.mean() - u0_fs) / u0_fs < PRE_RTOL, (
    f"Pre-shock u: {u_pre.mean():.4f} vs {u0_fs:.4f}")
assert abs(p_pre.mean() - p0_fs) / p0_fs < PRE_RTOL, (
    f"Pre-shock p: {p_pre.mean():.4f} vs {p0_fs:.4f}")

# post-shock assertions (after reflected shock)
assert abs(rho_post.mean() - rho3) / rho3 < POST_RTOL, (
    f"Post-shock rho: {rho_post.mean():.4f} vs {rho3:.4f}")
assert abs(u_post.mean() - ux3) / abs(ux3) < POST_RTOL, (
    f"Post-shock u: {u_post.mean():.4f} vs {ux3:.4f}")
assert abs(v_post.mean() - uy3) / abs(ux3) < POST_RTOL, (
    f"Post-shock v: {v_post.mean():.4f} vs {uy3:.4f}")
assert abs(p_post.mean() - p3) / p3 < POST_RTOL, (
    f"Post-shock p: {p_post.mean():.4f} vs {p3:.4f}")
```

Key design choices:
- `v` normalized by `abs(ux3)` (not `abs(uy3)`) because after a perfect wall reflection `uy3 ≈ 0`, making a relative tolerance meaningless.
- Spatial regions are chosen to avoid boundary layers, ghost cells, and the shock itself.

---

## Verification

```bash
# Unit-level: smoke-test the analytical functions without GPU
python -c "
from ouxsbli.tests.utils.oblique_shock import free_stream, reslected_shock, beta_from_theta, post_shock_state
import numpy as np
u0,p0,T0 = free_stream(2.0,1.4,287.03,1e5,295)
print('free_stream:', u0, p0, T0)
r3,ux3,uy3,p3 = reslected_shock(2.0,1.4,287.03,1e5,295,np.pi*37.2/180,np.pi*37.7/180)
print('reslected_shock rho3 ux3 uy3 p3:', r3, ux3, uy3, p3)
ref = post_shock_state(2.0, 8.0, gamma=1.4)
print('post_shock_state:', ref)
"

# Full integration test (requires GPU)
pytest ouxsbli/tests/test_os.py -v -m integration
```
