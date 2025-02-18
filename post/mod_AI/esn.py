import numpy as np
import optuna


class RidgeRegression:
  def __init__(self, tikh=1.e-5):
    self.tikh = tikh
    self.W    = None


  def fit(self, R, Y):
    LHS = R[1:].T @ R[1:] + self.tikh * np.eye(R[1:].shape[1])
    RHS = R[1:].T @ Y
    self.Wout = np.linalg.solve(LHS, RHS)



# ESN with bias architecture
class ESN:
  def __init__(self, n_units, in_dim, out_dim, sigma_in=0.5, rho=0.9e0, tikh=1e-5):
    self.n_units  = n_units # reservoir size
    self.in_dim   = in_dim
    self.out_dim  = out_dim
    
    # hyper parameters
    self.bias_in  = 1.0 # input bias
    self.bias_out = 1.0 # output bias 
    connectivity  = 2
    self.sigma_in = sigma_in # input scaling
    self.rho      = rho      # spectral radius
    sparseness    = 1 - connectivity/(n_units-1) 
    self.tikh     = tikh     # Tikhonov factor
    self.gamma    = 0.9e0
    
    # reservoir weight
    self.Win = np.zeros((in_dim+1, n_units))
    for j in range(self.n_units):
      self.Win[np.random.randint(0, in_dim+1),j] = np.random.uniform(-1, 1) #only one element different from zero per row
    # practical way to set the sparseness
    self.W = np.random.uniform(-1, 1, (n_units, n_units)) * (np.random.rand(n_units, n_units) < (1-sparseness))
    spectral_radius = np.max(np.abs(np.linalg.eigvals(self.W)))
    self.W /= spectral_radius #scaled to have unitary spec radius

    # output weight
    self.ridge = RidgeRegression(tikh=self.tikh)


  def step(self, r_pre, x):
    """ Advances one ESN time step.
      Args:
        x_pre: reservoir state(reservoir size)
        x:     input(dimension)
      Returns:
        new augmented state (new state with bias_out appended)
    """
    x_augmented = np.hstack([x, self.bias_in]) # input bias added
    # scaling input and weight
    W_in     = self.Win * self.sigma_in
    W_scaled = self.W   * self.rho
    r_new    = np.tanh(x_augmented @ W_in + W_scaled @ r_pre) 
    r_post   = (1.e0 - self.gamma) * r_pre + self.gamma * r_new
    return np.hstack([r_post, self.bias_out]) # output bias added


  def open_loop(self, X, r0):
    """ Advances ESN in open-loop / idle iteration.
      Args:
        X: input time series
        r0: initial reservoir state
      Returns:
        time series of augmented reservoir states
        X[0] to X[N-1] -> R[0] to R[N]
    """
    N = X.shape[0]
    R = np.empty((N + 1, self.n_units+1))
    R[0] = np.hstack([r0, self.bias_out])
    for i in 1 + np.arange(N):
      R[i] = self.step(R[i-1,:self.n_units], X[i-1])
    return R


  def train(self, X_washout, X_train, Y_train):
    """ Train the ESN.
      Args:
        X_washout: idle iteration input time series
        X_train: training input time series
    """
    ## washout phase / idle iteration
    rf_washout = self.open_loop(X_washout, np.zeros(self.n_units))[-1,:self.n_units]
    ## open-loop train phase
    R = self.open_loop(X_train, rf_washout)
    self.ridge.tikh = self.tikh
    self.ridge.fit(R, Y_train)
  

  def train_optim(self, X, N_washout, N_train, N_test, N_evo, N_tstart, offset):
    def predict(sigma_in, rho, tikh, gamma):
      self.sigma_in = sigma_in
      self.rho      = rho
      self.tikh     = tikh
      self.gamma    = gamma
      
      X_washout = X[:N_washout]
      X_train   = X[N_washout:N_washout+N_train-1]
      Y_train   = X[N_washout+1:N_washout+N_train]
      self.train(X_washout, X_train, Y_train)

      X_test_washout = X[N_tstart - N_washout + offset:N_tstart + offset]
      Y    = X[N_tstart + offset : N_tstart + offset + N_evo + 1]
      r0   = self.open_loop(X_test_washout, np.zeros(self.n_units))[-1]
      Yh   = self.evolve(r0, N_evo)
      loss = np.mean((Y - Yh)**2)
      return loss
    

    def optuna_objective(trial):
      sigma_in = trial.suggest_float('sigma_in', 0.49e0, 0.51e0)
      rho      = trial.suggest_float('rho',      0.89e0, 0.91e0)
      tikh     = trial.suggest_float('tikh',       1e-6, 1e-4)
      gamma    = trial.suggest_float('gamma',     0.890, 0.91e0)
      loss = predict(sigma_in, rho, tikh, gamma)
      return loss

    study = optuna.create_study(direction='minimize')
    study.optimize(optuna_objective, n_trials=2)
    sigma_in, rho, tikh, gamma = study.best_params.items()
    self.sigma_in = sigma_in[1]
    self.rho      = rho[1]
    self.tikh     = tikh[1]
    self.gamma    = gamma[1]



  def evolve(self, r0, N_evo):
    """Advances ESN in closed-loop
      Args:
        r0: initial reservoir state
        N_evo: Number of iterations
      Returns:
        time series of x_CL
        R[0] -> X[0] to X[N_evo]
    """
    R = np.empty((N_evo + 1, self.n_units+1))
    R[0] = r0
     
    Yh = np.zeros((N_evo+1, self.out_dim))
    Yh[0] = np.dot(R[0], self.ridge.Wout)
     
    for i in 1 + np.arange(N_evo):
      R[i] = self.step(R[i-1,:self.n_units], Yh[i-1])
      Yh[i] = np.dot(R[i], self.ridge.Wout) 
    return Yh


  def evolve_edge_removal(self, r0, j, i, N_evo):
    """Advances ESN in Intervened-loop
      Args:
        r0: initial reservoir state
        N_evo: Number of iterations
        j: dim index of the cause
        i: dim index of the effect
      Returns:
        time series of x_j->i
        R[0] -> X[0] to X[N_evo]
    """
    R = np.empty((N_evo + 1, self.n_units+1))
    Rprime = np.empty((N_evo + 1, self.n_units+1))
    R[0] = r0
    Rprime[0] = r0

    Yh = np.zeros((N_evo+1, self.out_dim))
    Yh[0] = np.dot(r0, self.ridge.Wout)

    Ain = np.ones(self.in_dim) # 0_j
    Ain[j] = 0

    Aout = np.ones(self.out_dim) # 0_i
    Aout[i] = 0

    Aoutprime = np.zeros(self.out_dim) # 1_i
    Aoutprime[i] = 1

    for i in 1 + np.arange(N_evo):
      R     [i] = self.step(R     [i-1,:self.n_units], Yh[i-1])
      Rprime[i] = self.step(Rprime[i-1,:self.n_units], Yh[i-1] * Ain)
      Yh[i] = np.dot(R[i], self.ridge.Wout) * Aout + np.dot(Rprime[i], self.ridge.Wout) * Aoutprime
    return Yh

