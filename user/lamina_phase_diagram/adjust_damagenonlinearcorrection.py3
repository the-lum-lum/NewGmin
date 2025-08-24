import sys
import os
import numpy as np


def update_curvature_history(curvature_folder, iteration, angles, num_history=3):
    """
    Saves the current curvature matrix and maintains a history of the last k curvature files.
    Each file is stored as .npy containing the angles.
    """
    os.makedirs(curvature_folder, exist_ok=True)
    hist_path = os.path.join(curvature_folder, f"curv_{iteration}.npy")
    # Save the current curvature matrix
    np.save(hist_path, angles)
    
    # Clean up old history files beyond the window
    for old_iter in range(iteration - num_history):
        old_file = os.path.join(curvature_folder, f"curv_{old_iter}.npy")
        if os.path.exists(old_file):
            os.remove(old_file)


def load_memory_matrix(curvature_folder, iteration, theta_rec, num_history=3):
    """
    Constructs a memory matrix (same shape as angles) indicating how many times
    each angle was below the 'recovery threshold' over the last `num_history` iterations.
    """
    memory_sum = None
    for i in range(iteration - num_history + 1, iteration + 1):
        hist_file = os.path.join(curvature_folder, f"curv_{i}.npy")
        if not os.path.exists(hist_file):
            continue
        angles = np.load(hist_file)
        below_thresh = (np.abs(angles) < theta_rec).astype(int)
        if memory_sum is None:
            memory_sum = below_thresh
        else:
            memory_sum += below_thresh
    return memory_sum


def compute_damage_with_memory(
    angles: np.ndarray,
    prev_stiffness: np.ndarray,
    memory_matrix: np.ndarray,
    *,
    k: int = 3,
    alpha: float = 0.05,
    alpha_nc: float = 1.0,
    delta: float = 2e-5,
    beta: float = 0.01,
    theta_rec: float = np.pi / 6,
    B_max: float = 0.015,
    max_angle: float = np.pi / 4
):
    """
    Curvature-only damage and memory-based recovery with nonlinear correction (Van Paepegem style).
    No network coupling or damage spikes.
    """
    n_fil, n_seg = angles.shape
    new_B = np.empty_like(prev_stiffness)
    tot_loc = 0.0
    for i in range(n_fil):
        for j in range(n_seg):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0.0, 1.0)

            # Local damage with nonlinear correction
            D_loc = alpha * gamma
            B_tmp = (prev_stiffness[i, j] * (1.0 - D_loc)) / (1.0 + alpha_nc * D_loc)

            # Memory-based recovery
            if memory_matrix[i, j] >= k:
                recovery = beta * (B_max - B_tmp) * np.exp(-theta / theta_rec)
                tot_loc += recovery
            else:
                recovery = 0.0

            new_B[i, j] = np.clip(B_tmp - delta + recovery, 1e-3, B_max)
    print(f"[recovery] local: {tot_loc:.4e}")
    return new_B


def adjust_bending_stiffness(
    iteration, num_filaments, N_SEG, 
    initial_stiffness=0.01,
    lower_rat=1.0, 
    alpha=0.05,         # used in compute_damage_with_memory => local damage factor
    alpha_nc=1.0,       # new param for the "nonlinear correction" in the denominator
    delta=0.0001, 
    beta=0.001,
    theta_rec=np.pi/6, 
    B_max=0.01, 
    max_angle=np.pi/4,
    memory_length=3, 
    curvature_dir="curvature_history"
):
    """
    Main pipeline for memory-based stiffness update with nonlinear correction.
    """
    num_bend = N_SEG - 1

    if iteration == 1:
        # First iteration => uniform stiffness
        b_array = np.full((num_filaments, num_bend), initial_stiffness)
        for i in range(num_filaments):
            # Example: scale every odd rod
            if i % 2 != 0:
                b_array[i, :] *= lower_rat

    else:
        # Load angles from previous iteration
        last_iteration_folder = f"Iteration_{iteration - 1}"
        lowest_file = os.path.join(last_iteration_folder, "lowest")
        if not os.path.isfile(lowest_file):
            print(f"Error: Could not find {lowest_file}. Using default stiffness.")
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

            # Load previous stiffness
            prev_b_file = os.path.join(last_iteration_folder, "b_array.in")
            if os.path.isfile(prev_b_file):
                with open(prev_b_file, 'r') as f:
                    prev_stiffness = np.loadtxt(f)
                    if prev_stiffness.ndim == 1:
                        prev_stiffness = prev_stiffness.reshape((num_filaments, num_bend))
            else:
                print(f"Warning: No b_array.in from iteration {iteration-1}, default used.")
                prev_stiffness = np.full((num_filaments, num_bend), initial_stiffness)
                for i in range(num_filaments):
                    if i % 2 != 0:
                        prev_stiffness[i, :] *= lower_rat

            # Update curvature history
            update_curvature_history(curvature_dir, iteration, angles, num_history=memory_length)

            # Load memory matrix from the last k steps
            memory_matrix = load_memory_matrix(curvature_dir, iteration, theta_rec, num_history=memory_length)
            if memory_matrix is None:
                memory_matrix = np.zeros_like(angles, dtype=int)

            # Now do the memory-based approach with a NONLINEAR correction
            b_array = compute_damage_with_memory(
                angles=angles,
                prev_stiffness=prev_stiffness,
                memory_matrix=memory_matrix,
                k=memory_length,
                alpha=alpha,            # damage factor
                alpha_nc=alpha_nc,      # non-linear correction param
                delta=delta,
                beta=beta,
                theta_rec=theta_rec,
                B_max=B_max,
                max_angle=max_angle
            )

    # Write out
    with open("b_array.in", "w") as file:
        np.savetxt(file, b_array, fmt="%.6f")

    print(f"Updated B_ARRAY for iteration {iteration} (with nonlinear correction).")
    print(b_array)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 adjust_damage.py <iteration> <num_filaments> <num_segments>")
        sys.exit(1)
    iteration = int(sys.argv[1])
    num_filaments = int(sys.argv[2])
    N_SEG = int(sys.argv[3])

    # Example usage with some defaults:
    adjust_bending_stiffness(
        iteration, 
        num_filaments, 
        N_SEG,
        initial_stiffness=0.01,
        lower_rat=1.0,
        alpha=0.05,            # local damage factor
        alpha_nc=1.0,          # new param for the rational denominator
        delta=5e-5, 
        beta=0.01,
        theta_rec=np.pi/6, 
        B_max=0.01, 
        max_angle=np.pi/4,
        memory_length=3, 
        curvature_dir="curvature_history"
    )
