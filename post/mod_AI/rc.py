import torch
import torch.nn as nn
import optuna


class RidgeRegression:
  def __init__(self, alpha=1.e-6):
    self.alpha   = alpha
    self.weights = None


  def fit(self, R, Y):
    """
    R: reservoir (timesteps, reservoir_size)
    Y: target    (timesteps, output_dim)
    """
    I = torch.eye(R.shape[1], device=R.device)
    self.weights = torch.linalg.solve(
      R.T @ R + self.alpha * I, R.T @ Y
    )


  def predict(self, R):
    if self.weights is None:
      raise ValueError("Model is not trained yet. Call 'fit' first")
    return R @ self.weights



class ESN(nn.Module):
  def __init__(self, input_dim, reservoir_size, output_dim, \
               spectral_radius=0.9e0, sparsity=0.1e0, ridge_alpha=1e-6):
    super(ESN, self).__init__()
    self.input_dim      = input_dim
    self.reservoir_size = reservoir_size
    self.output_dim     = output_dim
    # parameters
    self.spectral_radius = spectral_radius
    self.input_scaling   = 0.1e0
    self.ridge_alpha     = ridge_alpha
    self.gamma           = 0.9e0
    # input weight
    self.W_in = torch.randn(reservoir_size, input_dim)
    # reservoir weight
    W_re = torch.randn(reservoir_size, reservoir_size) * 0.1e0
    mask = torch.rand(reservoir_size, reservoir_size) > sparsity
    W_re[mask] = 0.e0
    eigvals    = torch.linalg.eigvals(W_re).abs().max()
    self.W_re  = W_re / eigvals
    # output weight
    self.ridge = RidgeRegression(alpha=self.ridge_alpha)


  def compute_reservoir_states(self, inputs):
    '''
    inputs: (batch_size, seq_len, input_dim)
    return: (total_timesteps, reservoir_size)
    '''
    batch_size, seq_len, _ = inputs.shape
    r = torch.zeros(batch_size, self.reservoir_size)
    states = []
    # scaling input weight
    W_in     = self.W_in * self.input_scaling
    W_scaled = self.W_re * self.spectral_radius
    for t in range(seq_len):
      u_t   = inputs[:,t,:]
      r_new = torch.tanh(W_in @ u_t.T + W_scaled @ r.T).T
      r     = (1.e0 - self.gamma) * r + self.gamma * r_new
      states.append(r.clone())
    return torch.stack(states, dim=1).reshape(-1, self.reservoir_size)


  def train_readout(self, train_inputs, train_targets, epochs=1, lr=1e-3):
    '''
    train weight of readout layer
    train_inputs : (batch_size, seq_len, input_dim)
    train_targets: (batch_size, seq_len, output_dim)
    '''
    def calc_RC(spectral_radius, input_scaling, ridge_alpha, gamma):
      self.spectral_radius = spectral_radius
      self.input_scaling   = input_scaling
      self.ridge_alpha     = ridge_alpha
      self.gamma           = gamma

      R = self.compute_reservoir_states(train_inputs)
      Y = train_targets.reshape(-1, self.output_dim)

      self.ridge.alpha = self.ridge_alpha
      self.ridge.fit(R, Y)

      predictions = self.ridge.predict(R)
      loss = torch.mean((predictions - Y)**2)
      return loss


    def optuna_objective(trial):
      spectral_radius = trial.suggest_float('spectral_radius', 0.9e0, 0.95e0)
      input_scaling   = trial.suggest_float('input_scaling', 0.09e0, 0.11e0)
      ridge_alpha     = trial.suggest_float('ridge_alpha', 1e-7, 1e-6)
      gamma           = trial.suggest_float('gamma', 0.85e0, 1.e0)

      loss = calc_RC(spectral_radius, input_scaling, ridge_alpha, gamma)
      return loss


    study = optuna.create_study(direction='minimize')
    study.optimize(optuna_objective, n_trials=100)
    spectral_radius, input_scaling, ridge_alpha, gamma = study.best_params.items()
    self.spectral_radius = spectral_radius[1]
    self.input_scaling   = input_scaling[1]
    self.ridge_alpha     = ridge_alpha[1]
    self.gamma           = gamma[1]


  def forward(self, inputs):
    '''
    inputs: (batch_size, seq_len, input_dim)
    '''
    R = self.compute_reservoir_states(inputs)
    return self.ridge.predict(R).reshape(inputs.shape[0], inputs.shape[1], self.output_dim)

