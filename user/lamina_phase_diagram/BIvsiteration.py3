import os
import numpy as np
import matplotlib.pyplot as plt

# === USER SETTINGS ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\lamina_NROD10_NSEG_200_RAT5_initial"  # Change this to your simulation directory
num_iterations = 20  # Or however many iterations you ran

# === Collect mean B values over iterations ===
mean_B_values = []
available_iters = []

for i in range(1, num_iterations + 1):
    folder = f"Iteration_{i}"
    file_path = os.path.join(base_path, folder, "b_array.in")

    if os.path.exists(file_path):
        b_array = np.loadtxt(file_path)
        if b_array.ndim == 1:
            print(f"Warning: {file_path} appears to be 1D. Skipping.")
            continue
        mean_B = np.mean(b_array)
        mean_B_values.append(mean_B)
        available_iters.append(i)
    else:
        print(f"Missing: {file_path}")

# === Plot ===
plt.figure(figsize=(8, 5))
plt.plot(available_iters, mean_B_values, marker='o', linestyle='-', linewidth=2)
plt.title("Mean Bending Stiffness $\\bar{B}^{(t)}$ Across Iterations")
plt.xlabel("Iteration")
plt.ylabel("Mean $\\bar{B}$")
plt.grid(True)
plt.tight_layout()
plt.show()
