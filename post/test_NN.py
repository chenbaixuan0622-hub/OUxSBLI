import numpy as np
from sklearn.neighbors import NearestNeighbors as NN

X = np.array([[-1,-1], [-2,-1], [-3,-2], [1,1], [2,1], [3,2]])
nrbs = NN(n_neighbors=3, algorithm='ball_tree').fit(X)
distances, indices = nrbs.kneighbors(X)

print(indices)
print(X[indices])
print(distances)

