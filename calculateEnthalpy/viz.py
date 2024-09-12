from matplotlib import pyplot as plt
import numpy as np

n = 2
density = 10
temp_list = []
for i in range(n):
	temp_list.append(np.linspace(0, 1, density))


g = np.meshgrid(*temp_list)
for idx, i in enumerate(g):
	g[idx] = i.reshape(-1,1)

coords = np.hstack(g)

filtered_coords = coords[np.sum(coords, axis=1) == 1]

print(coords.shape, filtered_coords.shape)

# Plot the 2D projection
plt.scatter(coords[:,0], coords[:,1], s=5, cmap='Spectral')
plt.scatter(filtered_coords[:,0], filtered_coords[:,1], s=5, c = 'black')
plt.title('UMAP projection of the 3D coordinates into 2D')
plt.show()