import sys
import os
import numpy as np

def update_curvature_history(curvature_folder, iteration, angles, num_history=3):
    """
    Saves the current curvature matrix and maintains a history of the last k curvature files.
    Each file is stored as a .npy containing the angles at iteration. 
    """
    os.makedirs(curvature_folder, exist_ok=True)
    hist_path = os.path.join(curvature_folder, f"curv_{iteration}.npy")

    # Save the current angles as a .npy
    np.save(hist_path, angles)

    # Cleanup older files beyond the window
    for old_iter in range(iteration - num_history):
        old_file = os.path.join(curvature_folder, f"curv_{old_iter}.npy")
        if os.path.exists(old_file):
            os.remove(old_file)

def load_memory_matrix(curvature_folder, iteration, theta_rec, num_history=3):
    """
    Constructs a matrix (same shape as angles) counting how many times
    each angle was below the threshold over the last `num_history` iterations.
    """
    memory_sum = None
    for i in range(iteration - num_history + 1, iteration + 1):
        hist_file = os.path.join(curvature_folder, f"curv_{i}.npy")
        if not os.path.exists(hist_file):
            continue
        hist_angles = np.load(hist_file)
        below_thresh = (np.abs(hist_angles) < theta_rec).astype(float)
        if memory_sum is None:
            memory_sum = below_thresh
        else:
            memory_sum += below_thresh
    return memory_sum

def apply_network_coupling(b_stiffness, kappa=0.02):
    """
    Applies a stiffness-based coupling damage, where each node
    is reduced if it is stiffer than its neighbors.
    
    b_stiffness: (m, n) array
    kappa: coupling coefficient
    Returns: coupling_damage (same shape as b_stiffness)
    """
    m, n = b_stiffness.shape
    # Prepare an array to accumulate coupling damage
    coupling_damage = np.zeros((m, n), dtype=float)
    
    def neighbors(i, j):
        # Example adjacency:
        #   up/down in i-direction => (i-1, j), (i+1, j)
        #   left/right in j-direction => (i, j-1), (i, j+1)
        #   diagonal neighbors => (i-1, j-1), (i-1, j+1), (i+1, j-1), (i+1, j+1)
        # Adjust or reduce as needed.
        nb_list = []
        for di, dj in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]:
            ip = i + di
            jp = j + dj
            if (0 <= ip < m) and (0 <= jp < n):
                nb_list.append((ip, jp))
        return nb_list

    # For each node, compute how much damage is introduced by neighbors
    for i in range(m):
        for j in range(n):
            # difference in stiffness from neighbors
            neighbors_list = neighbors(i, j)
            diff_sum = 0.0
            for (ip, jp) in neighbors_list:
                diff_sum += (b_stiffness[i,j] - b_stiffness[ip,jp])
            # total coupling damage is scaled by kappa
            coupling_damage[i,j] = kappa * diff_sum

    return coupling_damage

def compute_damage_with_memory(
    angles, prev_stiffness, memory_matrix, k=2,
    alpha=0.03, delta=0.0001, beta=0.0005,
    theta_rec=np.pi/6, B_max=0.01, max_angle=np.pi/4,
    kappa=0.02
):
    """
    Damage + conditional memory-based recovery + network coupling.
    1) Local curvature-based damage
    2) Memory-based partial recovery
    3) Network coupling damage
    """
    num_filaments, num_segments = angles.shape
    b_array = np.copy(prev_stiffness)  # Start from previous

    # 1) Local curvature-based damage + partial uniform degrade
    for i in range(num_filaments):
        for j in range(num_segments):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0, 1)  # normalised angle [0..1]

            # If memory_matrix[i,j] == k => all last k iters below threshold => allow recovery
            allow_recovery = (memory_matrix[i, j] == k)
            # Recovery factor: let's do a simple approach with an exponential or linear term
            R_ij = beta * (B_max - b_array[i, j]) if allow_recovery else 0.0

            # local damage
            local_damage = alpha * gamma + delta

            # update
            b_array[i, j] = b_array[i, j] - local_damage + R_ij
            if b_array[i, j] < 0.001:
                b_array[i, j] = 0.001
            if b_array[i, j] > B_max:
                b_array[i, j] = B_max

    # 2) Now apply network-based coupling damage
    #    1) compute how much damage from neighbors
    coupling_damage = apply_network_coupling(b_array, kappa=kappa)
    #    2) subtract it from b_array
    b_array = b_array - coupling_damage
    # clamp again
    b_array = np.clip(b_array, 0.001, B_max)

    return b_array

def adjust_bending_stiffness(
    iteration, num_filaments, N_SEG,
    initial_stiffness=0.01, lower_rat=1.0,
    alpha=0.05, delta=0.0001, beta=0.001,
    theta_rec=np.pi/6, B_max=0.01, max_angle=np.pi/4,
    memory_length=3, curvature_dir="curvature_history",
    kappa=0.02  # network coupling parameter
):
    """
    Main function to read angles from previous iteration, load previous b_array, 
    and update to new stiffness using local+memory+network coupling damage.
    """
    num_bend = N_SEG - 1

    if iteration == 1:
        # Initialize uniform stiffness
        b_array = np.full((num_filaments, num_bend), initial_stiffness)
        for i in range(num_filaments):
            if i % 2 != 0:
                b_array[i, :] *= lower_rat

    else:
        last_iter_folder = f"Iteration_{iteration-1}"
        lowest_file = os.path.join(last_iter_folder, "lowest")
        if not os.path.isfile(lowest_file):
            print(f"[WARN] Could not find {lowest_file}, using default.")
            b_array = np.full((num_filaments, num_bend), initial_stiffness)
            for i in range(num_filaments):
                if i % 2 != 0:
                    b_array[i, :] *= lower_rat
        else:
            angles_all = np.loadtxt(lowest_file, skiprows=4)
            try:
                angles = angles_all.reshape((num_filaments, num_bend))
            except Exception as e:
                print("Error reshaping angles:", e)
                angles = np.zeros((num_filaments, num_bend))
            
            # Read previous stiffness
            prev_b_file = os.path.join(last_iter_folder, "b_array.in")
            if os.path.isfile(prev_b_file):
                prev_stiffness_1d = np.loadtxt(prev_b_file)
                if prev_stiffness_1d.ndim == 1:
                    prev_stiffness = prev_stiffness_1d.reshape((num_filaments, num_bend))
                else:
                    prev_stiffness = prev_stiffness_1d
            else:
                prev_stiffness = np.full((num_filaments, num_bend), initial_stiffness)
                for i in range(num_filaments):
                    if i % 2 != 0:
                        prev_stiffness[i, :] *= lower_rat

            # Update curvature history
            update_curvature_history(curvature_dir, iteration, angles, num_history=memory_length)

            # Build memory matrix
            memory_matrix = load_memory_matrix(curvature_dir, iteration, theta_rec, num_history=memory_length)
            if memory_matrix is None:
                memory_matrix = np.zeros((num_filaments, num_bend))

            # Compute new b_array with local + memory + network damage
            b_array = compute_damage_with_memory(
                angles, 
                prev_stiffness, 
                memory_matrix, 
                k=memory_length,
                alpha=alpha, 
                delta=delta,
                beta=beta, 
                theta_rec=theta_rec, 
                B_max=B_max, 
                max_angle=max_angle,
                kappa=kappa
            )

    # Save new b_array
    with open("b_array.in", "w") as f:
        np.savetxt(f, b_array, fmt="%.6f")

    print(f"[INFO] iteration={iteration} updated B_ARRAY:\n", b_array)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 adjust_damage.py <iteration> <num_filaments> <num_segments>")
        sys.exit(1)
    
    iteration = int(sys.argv[1])
    num_filaments = int(sys.argv[2])
    N_SEG = int(sys.argv[3])

    adjust_bending_stiffness(iteration, num_filaments, N_SEG)
