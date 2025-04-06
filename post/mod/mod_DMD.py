import numpy as np
import matplotlib.pyplot as plt


class DMD:
  def __init__(self, X, dt):
    # X[space, time]
    self.X1 = X[:,:-1]
    self.X2 = X[:,1:]
    self.dt = dt


  def exec(self, rank):
    # X1, X2[space, time]
    U, S, Vh = np.linalg.svd(self.X1, full_matrices=False)
    Ur  = U[:,:rank]
    Sr  = np.diag(S[:rank])
    Vhr = Vh[:rank,:]
    Atilde = Ur.conj().T @ self.X2 @ Vhr.conj().T @ np.linalg.inv(Sr)
    self.Lambda, w = np.linalg.eig(Atilde)
    self.Omega = np.log(self.Lambda) / dt
    self.Phi = self.X2 @ Vhr.conj().T @ np.linalg.inv(Sr) @ w
    self.b = np.linalg.pinv(self.Phi) @ self.X1[:,0]
 
  
  def reconst(self, ranks, Nt):
    omega = np.imag(self.Omega)
    freq  = omega / (2.e0 * np.pi)
    time  = np.linspace(0.e0, 4.e0 * np.pi, Nt)
    size  = self.X1.shape[0]
    Xr    = np.zeros((size, Nt), dtype=np.float32)
    for i, t in enumerate(time):
      x = np.zeros((rows, 1))
      for rank in ranks:
        if freq[rank] == 0.e0:
          x = np.reshape(self.Phi[:,rank] * b[rank], (size, 1)) + x
        else:
          x = np.reshape(self.Phi[:,rank] * np.exp(complex(0, t)) * b[rank], (size, 1)) + x
      Xr[:,i:i+1] = x
    return Xr


  def plot_energy(self, file_path=None, figsize=(18,6)):
    _, S, _ = np.linalg.svd(self.X1, full_matrices=False)
    rank = []
    for i in range(len(S)):
      rank.append(i+1)
    cumlative_cont = (S / np.sum(S)) * 100
    fig, ax = plt.subplots(2, 1, figsize=figsize)
    ax[0].semilogy(rank, S)
    ax[0].set_xlabel("Rank")
    ax[0].set_ylabel("Singular value")
    ax[1].plot(rank, cumlative_cont)
    ax[0].set_xlabel("Rank")
    ax[0].set_ylabel("Cumulative Contribution [%]")
    if file_path is not None:
      plt.savefig(file_path)
    plt.show()


  def plot_eig(self, file_path=None, figsize=(8,6)):
    if self.Lambda is None:
      raise ValueError("call 'exec' first")
    theta = np.linspace(0.e0, 2.e0 * np.pi, 150)
    x = np.cos(theta)
    y = np.sin(theta)
    fig, ax = plt.subplots(figsize=figsize, dpi=100)
    ax.scatter(np.real(self.Lambda), np.imag(self.Lambda))
    ax.plot(x, y, color='black', ls='-')
    ax.set_xlim(-1.1e0, 1.1e0)
    ax.set_ylim(-1.1e0, 1.1e0)
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    ax.set_aspect('equal')
    if file_path is not None:
      plt.savefig(file_path)
    plt.show()

