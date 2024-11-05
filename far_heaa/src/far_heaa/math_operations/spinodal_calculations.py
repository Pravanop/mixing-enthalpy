import numpy as np
from scipy.spatial.distance import pdist


# Step 1: Calculate a dynamic tolerance based on grid spacing
def calculate_dynamic_tolerance(mol_grid, factor=1.0):
    distances = pdist(mol_grid)
    median_distance = np.median(distances)
    tolerance = factor * median_distance
    return tolerance


# Step 2: Find the closest neighbors for each mol fraction direction and return their indices
def find_closest_neighbor_indices(mol_fractions, mol_grid, factor=0.5):
    tolerance = calculate_dynamic_tolerance(mol_grid, factor)
    n = len(mol_fractions)
    closest_neighbors = [{} for _ in range(n)]

    for i in range(n):
        # Create a mask to find points that match on all fractions except the ith one
        mask = np.all(np.isclose(mol_grid[:, np.arange(n) != i], mol_fractions[np.arange(n) != i], atol=tolerance),
                      axis=1)
        potential_neighbors = mol_grid[mask]
        potential_indices = np.where(mask)[0]  # Get indices of potential neighbors

        # Identify forward and backward neighbors for the ith fraction
        forward_index = None
        backward_index = None
        min_forward_distance = float('inf')
        min_backward_distance = float('inf')

        for idx, point in zip(potential_indices, potential_neighbors):
            delta = point[i] - mol_fractions[i]
            distance = np.linalg.norm(point - mol_fractions)

            if delta > 0 and distance < min_forward_distance:
                min_forward_distance = distance
                forward_index = idx
            elif delta < 0 and distance < min_backward_distance:
                min_backward_distance = distance
                backward_index = idx

        closest_neighbors[i]["forward_index"] = forward_index
        closest_neighbors[i]["backward_index"] = backward_index
        closest_neighbors[i]["delta"] = min_forward_distance if forward_index is not None else min_backward_distance

    return closest_neighbors


def calculate_hessian(mol_fractions, gibbs_value, mol_grid, gibbs, tolerance_factor=0.5):
    n = len(mol_fractions)
    hessian = np.zeros((n, n))

    # Find the closest neighbors' indices for each direction
    closest_neighbors = find_closest_neighbor_indices(mol_fractions, mol_grid, factor=tolerance_factor)

    # Diagonal elements: Second partial derivatives with respect to each mol fraction
    for i in range(n):
        forward_index = closest_neighbors[i]["forward_index"]
        backward_index = closest_neighbors[i]["backward_index"]

        if forward_index is not None and backward_index is not None:
            forward_gibbs = gibbs[forward_index]
            backward_gibbs = gibbs[backward_index]
            delta = closest_neighbors[i]["delta"]
            hessian[i, i] = (forward_gibbs - 2 * gibbs_value + backward_gibbs) / (delta ** 2)

    # Off-diagonal elements: Mixed partial derivatives
    for i in range(n):
        for j in range(i + 1, n):
            # Find neighbors that vary both i and j simultaneously in all combinations
            f_pp_index = closest_neighbors[i].get("forward_index")
            f_pm_index = closest_neighbors[j].get("backward_index")
            f_mp_index = closest_neighbors[i].get("backward_index")
            f_mm_index = closest_neighbors[j].get("backward_index")

            if all(index is not None for index in [f_pp_index, f_pm_index, f_mp_index, f_mm_index]):
                f_pp = gibbs[f_pp_index]
                f_pm = gibbs[f_pm_index]
                f_mp = gibbs[f_mp_index]
                f_mm = gibbs[f_mm_index]

                delta_i = closest_neighbors[i]["delta"]
                delta_j = closest_neighbors[j]["delta"]
                hessian[i, j] = hessian[j, i] = (f_pp - f_pm - f_mp + f_mm) / (4 * delta_i * delta_j)

    return hessian

def calculate_eigenvalues(hessian):

    eigenvalues = np.linalg.eigvals(hessian)
    return eigenvalues