import sys
import os
import numpy as np





def update_curvature_history(curvature_folder, iteration, angles, num_history=3):
    """
    Saves the current curvature matrix and maintains a history of the last k curvature files.
    Each file is stored as a binary indicator of whether curvature is below the threshold.
    """
    import os
    import numpy as np

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
    each angle was below the recovery threshold over the last `num_history` iterations.
    """
    import numpy as np
    import os

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
    delta: float = 2e-5,
    beta: float = 0.01,
    theta_rec: float = np.pi / 6,
    B_max: float = 0.015,
    max_angle: float = np.pi / 4
):
    """
    Curvature-only damage and memory-based recovery without network coupling or nonlinear correction.
    """
    n_fil, n_seg = angles.shape
    new_B = np.empty_like(prev_stiffness)

    for i in range(n_fil):
        for j in range(n_seg):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0.0, 1.0)

            # Damage term
            D_loc = alpha * gamma
            B_tmp = prev_stiffness[i, j] * (1.0 - D_loc)

            # Memory-based recovery
            if memory_matrix[i, j] >= k:
                recovery = beta * (B_max - B_tmp) * np.exp(-theta / theta_rec)
            else:
                recovery = 0.0

            # Apply degradation, correction offset, and recovery
            new_B[i, j] = np.clip(B_tmp - delta + recovery, 1e-3, B_max)

    return new_B




def adjust_bending_stiffness(iteration, num_filaments, N_SEG, initial_stiffness=0.01,
                              lower_rat=1.0, alpha=0.05, delta=0.0001, beta=0.001,
                              theta_rec=np.pi/6, B_max=0.01, max_angle=np.pi/4,
                              memory_length=3, curvature_dir="curvature_history"):

    num_bend = N_SEG - 1

    if iteration == 1:
        # Initialisation: assign uniform bending stiffness
        b_array = np.full((num_filaments, num_bend), initial_stiffness)
        for i in range(num_filaments):
            if i % 2 != 0:
                b_array[i, :] *= lower_rat

    else:
        # Load angles from previous run
        last_iteration_folder = f"Iteration_{iteration - 1}"
        lowest_file = os.path.join(last_iteration_folder, "lowest")
        if not os.path.isfile(lowest_file):
            print(f"Error: Could not find {lowest_file}. Using default stiffness values.")
            b_array = np.full((num_filaments, num_bend), initial_stiffness)
            for i in range(num_filaments):
                if i % 2 != 0:
                    b_array[i, :] *= lower_rat
        else:
            with open(lowest_file, 'r') as f:
                angles_all = np.loadtxt(f, skiprows=4)
            try:
                angles = angles_all.reshape((num_filaments, num_bend))
            except Exception as e:
                print("Error reshaping angles:", e)
                angles = np.zeros((num_filaments, num_bend))

            # Load previous bending stiffness
            prev_b_file = os.path.join(last_iteration_folder, "b_array.in")
            if os.path.isfile(prev_b_file):
                with open(prev_b_file, 'r') as f:
                    prev_stiffness = np.loadtxt(f)
                    if prev_stiffness.ndim == 1:
                        prev_stiffness = prev_stiffness.reshape((num_filaments, num_bend))
            else:
                print(f"Warning: Previous b_array.in not found. Using default stiffness.")
                prev_stiffness = np.full((num_filaments, num_bend), initial_stiffness)
                for i in range(num_filaments):
                    if i % 2 != 0:
                        prev_stiffness[i, :] *= lower_rat

            # Update curvature history files
            update_curvature_history(curvature_dir, iteration, angles, num_history=memory_length)

            # Load memory matrix from last k curvature frames
            memory_matrix = load_memory_matrix(curvature_dir, iteration, theta_rec, num_history=memory_length)

            # Compute updated bending stiffness using cumulative damage memory rule
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
                max_angle=max_angle
            )

    # Save updated stiffness array
    with open("b_array.in", "w") as file:
        np.savetxt(file, b_array, fmt="%.6f")

    print(f"Updated B_ARRAY for iteration {iteration}:")
    print(b_array)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 adjust_damage.py <iteration> <num_filaments> <num_segments>")
        sys.exit(1)
    iteration = int(sys.argv[1])
    num_filaments = int(sys.argv[2])
    N_SEG = int(sys.argv[3])
    adjust_bending_stiffness(
        iteration, 
        num_filaments, 
        N_SEG,
        initial_stiffness=0.01,
        lower_rat=1.0,
        alpha=0.05,            # local damage factor        
        delta= 5e-5, 
        beta=0.01,
        theta_rec=np.pi/6, 
        B_max=0.01, 
        max_angle=np.pi/4,
        memory_length=3, 
        curvature_dir="curvature_history"
    )
