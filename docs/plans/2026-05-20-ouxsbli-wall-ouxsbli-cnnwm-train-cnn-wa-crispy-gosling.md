# Plan: Fix NaN Training and Upgrade to ResNet Architecture

## Context

Training diverges to NaN at epoch ~30.  Root cause: the per-channel inverse-variance
weighting added to the tau loss computes `w = 1 / var(τ_wz_norm)`.  For a TBL,
τ_wz ≈ 0 everywhere, so `var(τ_wz_norm) ~ 1e-10`, giving `w_τ_wz ~ 1e10`.  The `.clamp(min=1e-8)`
guard still allows weights up to `1e8`.  Any nonzero τ_wz prediction in early training
produces a gradient of O(1e8) → overflow → NaN.

Additionally, the flat sequential architecture (no skip connections) is prone to
gradient vanishing/exploding in deeper networks.  Replacing it with a small ResNet
(stem + 2 residual blocks + head) fixes both gradient flow and overall accuracy.

## Files to Modify

| File | Changes |
|------|---------|
| `cnnwm/train_cnn_wall_model.py` | Replace architecture with ResNet, remove weighted loss, add gradient clipping |

`config.yaml`, `dataloader.py`, `plot_cnn_wall_model.py` are unchanged.

---

## Step 1 — Replace `WallModelCNN` with a ResNet

Add a private `_ResBlock` module and rewrite `WallModelCNN` to use it.
Keep the same public interface (`n_out`, `dropout`).

```python
import torch.nn.functional as F

class _ResBlock(nn.Module):
    def __init__(self, ch: int) -> None:
        super().__init__()
        self.conv1 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.bn1   = nn.BatchNorm2d(ch)
        self.conv2 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.bn2   = nn.BatchNorm2d(ch)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        out = F.gelu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        return F.gelu(out + residual)


class WallModelCNN(nn.Module):
    """
    2-D ResNet wall model.

    Stem (3×3) + two residual blocks (each 3×3+3×3) + 1×1 head.
    Effective receptive field: 11×11.
    Parameters: ~150 k for n_out=2.
    """
    def __init__(self, n_out: int, dropout: float = 0.1) -> None:
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv2d(4, 64, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(64),
            nn.GELU(),
        )
        self.body = nn.Sequential(
            _ResBlock(64),
            _ResBlock(64),
        )
        self.head = nn.Sequential(
            nn.Dropout2d(p=dropout),
            nn.Conv2d(64, n_out, kernel_size=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.body(self.stem(x)))
```

Add `import torch.nn.functional as F` to the imports section.

---

## Step 2 — Remove per-channel weighting, add gradient clipping

### 2a. `train()` — drop `ch_weights`, add `clip_grad_norm_`

```python
def train(model: nn.Module, loader: DataLoader,
          n_epochs: int, lr: float, device: torch.device) -> list[float]:
    model.to(device)
    optimiser = torch.optim.RAdam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, T_max=n_epochs)
    losses: list[float] = []
    for epoch in range(1, n_epochs + 1):
        model.train()
        total = 0.0
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimiser.zero_grad()
            loss = F.mse_loss(model(xb), yb)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimiser.step()
            total += loss.item() * xb.size(0)
        scheduler.step()
        epoch_loss = total / len(loader.dataset)
        losses.append(epoch_loss)
        if epoch % max(1, n_epochs // 10) == 0:
            print(f"  epoch {epoch:4d}/{n_epochs}  loss={epoch_loss:.6e}"
                  f"  lr={scheduler.get_last_lr()[0]:.2e}")
    return losses
```

### 2b. `main()` — remove weight computation and passing

Delete the block:
```python
var_tau = torch.from_numpy(Y_tau_train.var(axis=0)).float().clamp(min=1e-8)
w_tau   = (1.0 / var_tau).view(1, 2, 1, 1)
var_qw  = torch.from_numpy(Y_qw_train.var(axis=0)).float().clamp(min=1e-8)
w_qw    = (1.0 / var_qw).view(1, 1, 1, 1)
```

Change train calls to drop `ch_weights`:
```python
losses_tau = train(model_tau, loader_tau, args.n_epochs, args.lr, device)
losses_qw  = train(model_qw,  loader_qw,  args.n_epochs, args.lr, device)
```

---

## What stays unchanged

- `make_patch_dataset` with z-flip augmentation
- `CosineAnnealingLR` scheduler
- Patch-based evaluation in `main()`
- `compute_metrics`, `compute_metrics` call sites
- All CLI args and config values (`n_epochs=300`, `dropout=0.1`)
- `torch.jit.trace` export

---

## Why this works

| Problem | Fix |
|---------|-----|
| NaN at epoch 30 | Removing per-channel `1/var` weighting eliminates O(1e8) gradient spikes |
| Gradient instability in deeper net | Residual skip connections provide a stable gradient path back |
| Spurious large-negative τ_wz | ResNet + standard MSE + dropout → model regularised enough to output ≈ 0 for near-zero target |
| q_w underfitting | Larger capacity (150k vs 57k), 11×11 RF, 300 epochs |

---

## Verification

```bash
cd cnnwm
python train_cnn_wall_model.py --config config.yaml
# Confirm: NO NaN loss at any epoch
# Confirm: tau loss descends smoothly, no spikes
# Confirm: q_w loss continues declining past epoch 100

python plot_cnn_wall_model.py --results results.npz --out-dir ./plots
# Confirm: scatter τ_wx R² meaningfully positive
# Confirm: scatter τ_wz R² > 0 (no longer catastrophically negative)
# Confirm: wall map τ_wz near-zero matching reference
```
