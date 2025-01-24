import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mod_AI.rc import ReservoirComputing


# 1. Lorenzモデルを生成
def lorenz_system(t, state, sigma=10, beta=8/3, rho=28):
  x, y, z = state
  dx = sigma * (y - x)
  dy = x * (rho - z) - y
  dz = x * y - beta * z
  return [dx, dy, dz]


def generate_lorenz_data(initial_state, t_span, t_eval):
  sol = solve_ivp(lorenz_system, t_span, initial_state, t_eval=t_eval, method="RK45")
  return sol.y.T  # shape: (timesteps, 3)


# 初期条件と時間設定
initial_state = [1.0, 1.0, 1.0]
t_span = (0, 40)  # 時間範囲
t_eval = np.linspace(0, 40, 4000)  # 時間ステップ
data = generate_lorenz_data(initial_state, t_span, t_eval)

# 正規化
data = (data - np.mean(data, axis=0)) / np.std(data, axis=0)

# トレーニングとテストデータに分割
train_data = data[:3000]
test_data = data[3000:]

tau = 1

# 3. モデルのトレーニング
# データ準備
train_input = torch.tensor(train_data[:-tau, 0:1], dtype=torch.float32).unsqueeze(0)  # shape: (1, timesteps, 3)
train_target = torch.tensor(train_data[tau:, 1:2], dtype=torch.float32).unsqueeze(0)  # shape: (1, timesteps, 3)

test_input = torch.tensor(test_data[:-tau, 0:1], dtype=torch.float32).unsqueeze(0)
test_target = torch.tensor(test_data[tau:, 1:2], dtype=torch.float32).unsqueeze(0)

# モデル初期化
reservoir_size = 30
ridge_alpha    = 1e-6
model = ReservoirComputing(
                           input_dim=train_input.shape[-1],
                           reservoir_size=reservoir_size,
                           output_dim=train_target.shape[-1],
                           ridge_alpha=ridge_alpha
                           )


# トレーニングループ
model.train_readout(train_input, train_target, epochs=2)

# 4. テスト
with torch.no_grad():
  predictions = model(test_input).squeeze(0).numpy()

# 5. 結果プロット
plt.figure(figsize=(12, 6))
plt.plot(test_data[100:, 0], label="True X")
plt.plot(predictions[100:, 0], label="Predicted X", linestyle="dashed")
plt.legend()
plt.title("Lorenz System Prediction (X) - Simple RC")
plt.show()

