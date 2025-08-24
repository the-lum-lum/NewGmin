import numpy as np
import matplotlib.pyplot as plt
import os

# === USER SETTINGS ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 spike beta0.01, delta 2e-5, k=2"
iter_prev = 25  # previous iteration
iter_curr = 35 # current iteration

# === Load Files ===
b_prev_path = os.path.join(base_path, f"Iteration_{iter_prev}", "b_array.in")
b_curr_path = os.path.join(base_path, f"Iteration_{iter_curr}", "b_array.in")

if not os.path.exists(b_prev_path) or not os.path.exists(b_curr_path):
    print("Error: One of the files does not exist.")
else:
    B_prev = np.loadtxt(b_prev_path)
    B_curr = np.loadtxt(b_curr_path)

    if B_prev.ndim == 1 or B_curr.ndim == 1:
        print("Warning: One of the arrays is 1D. Reshape may be required.")
    elif B_prev.shape != B_curr.shape:
        print("Mismatch in array shapes.")
    else:
        delta_B = B_curr - B_prev

        plt.figure(figsize=(10, 6))
        im = plt.imshow(delta_B, cmap='coolwarm', vmin=-0.001, vmax=0.001, aspect='auto')
        plt.title(f"Bending Stiffness Change (Iteration {iter_prev} → {iter_curr})")
        plt.xlabel("Segment Index $j$")
        plt.ylabel("Filament Index $i$")
        cbar = plt.colorbar(im, label="Change in $B_{ij}$")
        plt.tight_layout()
        plt.show()
