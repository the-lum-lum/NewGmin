import sys
import os
import numpy as np
from numpy.random import default_rng         # already used below


def compute_damage(angles, prev_stiffness, damage_factor, linear_decrement, max_angle=np.pi/4, alpha_nc=1.0, Bmin=0.001):
    """
    Nonlinear damage update rule:
    For each angle θ_ij, compute a damage variable D_ij = damage_factor * |θ_ij| / max_angle
    Then apply nonlinear correction as:
        B_ij = ((1 - D_ij) * B_prev) / (1 + alpha_nc * D_ij) - delta
    """
    num_filaments, num_segments = angles.shape
    b_array = np.zeros((num_filaments, num_segments))
    for i in range(num_filaments):
        for j in range(num_segments):
            angle_val = abs(angles[i, j])
            gamma_ij = min(angle_val / max_angle, 1.0)
            D_ij = damage_factor * gamma_ij
            numerator = (1 - D_ij) * prev_stiffness[i, j]
            denominator = 1 + alpha_nc * D_ij
            new_val = numerator / denominator - linear_decrement
            b_array[i, j] = max(new_val, Bmin)
    return b_array



def adjust_bending_stiffness(iteration, num_filaments, N_SEG, initial_stiffness=0.01, lower_rat=1.0, damage_factor=0.05, linear_decrement=5e-5):
    # Number of bending segments per filament
    # Number of bending segments per filament
    num_bend = N_SEG - 1

    if iteration == 1:
        # --- truncated‑normal initial stiffness -----------------------------
        mean_B  = initial_stiffness           # 0.01
        spread  = 0.15                        # ±10 %
        sigma   = 0.05 * mean_B              # 3 % s.d. before truncation
        low, hi = mean_B * (1 - spread), mean_B * (1 + spread)

        rng = default_rng()                  # remove argument for full randomness
        draws = rng.normal(loc=mean_B, scale=sigma, size=(num_filaments, num_bend))
        random_base = np.clip(draws, low, hi)   # hard truncation

        

        b_array = random_base.copy()
# -------------------------------------------------------------------------

    else:
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
            # Read previous bending stiffness (now a 2D array) from previous iteration's b_array.in
            prev_b_file = os.path.join(last_iteration_folder, "b_array.in")
            if os.path.isfile(prev_b_file):
                with open(prev_b_file, 'r') as f:
                    prev_stiffness = np.loadtxt(f)
                    if prev_stiffness.ndim == 1:
                        # If loaded as 1D, reshape it appropriately
                        prev_stiffness = prev_stiffness.reshape((num_filaments, num_bend))
            else:
                print(f"Warning: Previous b_array.in not found in {last_iteration_folder}. Using default stiffness.")
                prev_stiffness = np.full((num_filaments, num_bend), initial_stiffness)
                for i in range(num_filaments):
                    if i % 2 != 0:
                        prev_stiffness[i, :] *= lower_rat
            b_array = compute_damage(angles, prev_stiffness,
                         damage_factor=damage_factor,
                         linear_decrement=linear_decrement,alpha_nc=1.0)


    # Write the updated 2D bending stiffness to file.
    with open('b_array.in', 'w') as file:
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
    adjust_bending_stiffness(iteration, num_filaments, N_SEG)