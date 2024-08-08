def get_nearest_neighbors_periodic_3d(arr, x, y, z):
    max_x, max_y, max_z = arr.shape
    neighbors = []

    # Iterate over all possible neighbors
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                # Skip the center point itself
                if dx == 0 and dy == 0 and dz == 0:
                    continue

                # Calculate new indices with modulo for periodic boundary conditions
                nx, ny, nz = (x + dx) % max_x, (y + dy) % max_y, (z + dz) % max_z

                if arr[nx, ny, nz] != 0:
                    # Append the indices of the neighbor
                    neighbors.append((nx, ny, nz))

    return neighbors


def get_second_nearest_neighbors_periodic(arr, x, y, z):
    max_x, max_y, max_z = arr.shape
    neighbors = []

    # List all combinations that sum up to 2 steps away including diagonals
    # We consider steps like (2,0,0), (0,2,0), (0,0,2), (-2,0,0), (0,-2,0), (0,0,-2)
    # and all one-step diagonals like (1,1,0), (1,-1,0), etc.
    second_steps = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (-2, 0, 0), (0, -2, 0), (0, 0, -2),
                    (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
                    (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
                    (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1)]

    for dx, dy, dz in second_steps:
        # Calculate new indices with modulo for periodic boundary conditions
        nx, ny, nz = (x + dx) % max_x, (y + dy) % max_y, (z + dz) % max_z

        # Append the indices of the neighbor
        if arr[nx, ny, nz] != 0:
            # Append the indices of the neighbor
            neighbors.append((nx, ny, nz))

    return neighbors


def create_neighbor_list(arr, flag):
    max_x, max_y, max_z = arr.shape
    neighbor_list = []

    # Loop through each index in the array
    for x in range(max_x):
        for y in range(max_y):
            for z in range(max_z):
                # Get neighbors for the current index
                if arr[x][y][z] != 0:
                    if flag == 1:
                        neighbors = get_nearest_neighbors_periodic_3d(arr, x, y, z)
                    elif flag == 2:
                        neighbors = get_second_nearest_neighbors_periodic(arr, x, y, z)
                    # Append the current index and its neighbors
                    neighbor_list.append([(x, y, z), neighbors])

    return neighbor_list


