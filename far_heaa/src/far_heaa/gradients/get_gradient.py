import itertools

from far_heaa.grids_and_combinations import grid_creation
import numpy as np


mol_grid = grid_creation.CompositionGrid.create_mol_grid(4, 20)

y_grid = np.sin(mol_grid[:,  0])

print(mol_grid.shape, y_grid.shape)

def get_partial_gradients(mol_grid,
                          y_grid):

    gradients = []
    dx_n = np.gradient(mol_grid[:, -1])
    dx_n += 1e-9
    dy_dx_n = np.gradient(y_grid) / dx_n
    dy_dx_n += 1e-9
    for i in range(mol_grid.shape[-1]):

        dx = np.gradient(mol_grid[:, i])
        dx += 1e-9
        dy_dx = np.gradient(y_grid)/dx
        dy_dx += 1e-9
        gradients.append(dy_dx - dy_dx_n)

    gradients = np.array(gradients)

    return gradients

def create_direction_vector(k,
                            n,
                            counter):

    assert len(counter) == k
    assert k < n

    vec = np.zeros((n , 1))
    vec[counter] = (1/k)

    return vec

def compile_direction_vectors(
    n,
    k
):

    ind = list(range(n))

    direction_vector_list = []
    combination = np.array(list(itertools.combinations(ind, r = k)))
    for i in combination:
        direction_vector_list.append(create_direction_vector(k, n, i))

    direction_vector_list = np.array(direction_vector_list)
    return direction_vector_list

def get_sphere_indices(mol_grid, n, p):
    indices = []
    sphere_condition = np.sum(mol_grid**2, 1) - 1/n
    indices = np.where(sphere_condition <= p**2)[0]
    return indices


indices = get_sphere_indices(mol_grid, 4, 0.1)

gradients = get_partial_gradients(mol_grid, y_grid).T
print(gradients.shape)
constrained_gradients = gradients[indices]

def mega_compile_direction_vectors(n):
    direction_vectors = []
    for i in range(2, n):
        direction_vectors.append(compile_direction_vectors(n, i))

    # direction_vectors = np.array(direction_vectors)

    return direction_vectors


direc = mega_compile_direction_vectors(4)

def get_directional_derivative(gradient, direction):
    return np.dot(gradient, direction)



for i in direc:
    for j in i:
        print(get_directional_derivative(constrained_gradients[-1], j))
