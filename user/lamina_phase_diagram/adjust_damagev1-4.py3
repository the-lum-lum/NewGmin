import sys
import os
import numpy as np

def compute_damage(angles, prev_stiffness, damage_factor, linear_decrement, max_angle=np.pi/4):
    """
    Computes updated bending stiffness for each filament and each segment based on:
      - the local absolute angle (normalized by max_angle) to compute a multiplicative damage factor,
      - and then subtracts a fixed linear decrement.
    angles and prev_stiffness are arrays of shape (num_filaments, num_segments).
    """
    num_filaments, num_segments = angles.shape
    b_array = np.zeros((num_filaments, num_segments))
    for i in range(num_filaments):
        for j in range(num_segments):
            angle_val = abs(angles[i, j])
            angle_effect = np.clip(angle_val / max_angle, 0, 1)
            damage_multiplier = 1 - angle_effect * damage_factor
            damage_multiplier = max(damage_multiplier, 0.05)
            new_val = prev_stiffness[i, j] * damage_multiplier - linear_decrement
            b_array[i, j] = max(new_val, 0.001)
    return b_array


def compute_damage_v2(angles, prev_stiffness, a=10, theta_0=np.pi/6, linear_decrement=0.0001):
    num_filaments, num_segments = angles.shape
    b_array = np.zeros((num_filaments, num_segments))
    for i in range(num_filaments):
        for j in range(num_segments):
            theta = abs(angles[i, j])
            gamma = 1 / (1 + np.exp(-a * (theta - theta_0)))  # sigmoid function
            new_val = prev_stiffness[i, j] * (1 - gamma) - linear_decrement
            b_array[i, j] = max(new_val, 0.001)
    return b_array


def compute_damage_v3(angles, prev_stiffness, damage_factor, delta_d=0.0001, delta_r=0.00005, max_angle=np.pi/4):
    num_filaments, num_segments = angles.shape
    b_array = np.zeros((num_filaments, num_segments))
    
    for i in range(num_filaments):
        for j in range(num_segments):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0, 1)  # curvature effect
            damage_multiplier = 1 - damage_factor * gamma
            
            # Apply damage and recovery
            new_val = prev_stiffness[i, j] * damage_multiplier - delta_d + delta_r
            b_array[i, j] = max(new_val, 0.001)
    
    return b_array

def compute_damage_v4(angles, prev_stiffness, iteration,
                      alpha=0.05, delta=0.0001, beta=0.0005,
                      theta_rec=np.pi/6, B_max=0.01, max_angle=np.pi/4):
    """
    Damage + damage-weighted conditional recovery model.
    Recovery is proportional to damage and suppressed at high curvature.
    """
    num_filaments, num_segments = angles.shape
    b_array = np.zeros((num_filaments, num_segments))

    for i in range(num_filaments):
        for j in range(num_segments):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0, 1)

            # Damage
            damaged_val = (1 - alpha * gamma) * prev_stiffness[i, j] - delta

            # Recovery only if curvature is low and some damage has occurred
            if iteration >= 3 and theta < theta_rec:
                Rij = np.exp(-theta / theta_rec)
                recovery = beta * (B_max - prev_stiffness[i, j]) * Rij
            else:
                recovery = 0.0

            updated_val = damaged_val + recovery
            b_array[i, j] = min(B_max, max(updated_val, 0.001))

    return b_array




def adjust_bending_stiffness(iteration, num_filaments, N_SEG, initial_stiffness=0.01, lower_rat=1.0, damage_factor=0.05, linear_decrement=0.0001):
    # Number of bending segments per filament
    num_bend = N_SEG - 1
    if iteration == 1:
        # Initialize a 2D array: each filament gets the same stiffness along all bending segments,
        # with secondary rods (odd-indexed) scaled by lower_rat.
        b_array = np.full((num_filaments, num_bend), initial_stiffness)
        for i in range(num_filaments):
            if i % 2 != 0:
                b_array[i, :] *= lower_rat
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
            b_array = compute_damage_v4(angles, prev_stiffness,iteration)

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
