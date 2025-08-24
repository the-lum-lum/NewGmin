import numpy as np
import matplotlib.pyplot as plt
import os

# === USER SETTINGS ===
iteration = 70  # Change this to the iteration you want to visualise
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.0001 100iter damage spikes"
b_path = os.path.join(base_path, f"Iteration_{iteration}", "b_array.in")
spike_path = os.path.join(base_path, f"Iteration_{iteration}", "spike_mask.txt")

# === Load Data ===
if not os.path.exists(b_path) or not os.path.exists(spike_path):
    print("Missing files. Check paths.")
else:
    B = np.loadtxt(b_path)
    if B.ndim == 1:
        raise ValueError("B array is 1D. Check reshape logic.")

    spike_mask = np.loadtxt(spike_path)
    if spike_mask.shape != B.shape:
        raise ValueError("spike_mask and b_array shape mismatch")

    # === Plot ===
    plt.figure(figsize=(10, 6))
    plt.imshow(B, cmap='viridis', aspect='auto')
    plt.colorbar(label="Bending Stiffness $B_{ij}$")
    plt.title(f"Damage Heatmap with Random Spike Overlay — Iteration {iteration}")
    plt.xlabel("Segment Index $j$")
    plt.ylabel("Filament Index $i$")

    # === Overlay spikes as red X markers ===
    spike_indices = np.argwhere(spike_mask == 1)
    for (i, j) in spike_indices:
        plt.plot(j, i, 'rx', markersize=6, markeredgewidth=1.5)

    plt.tight_layout()
    plt.show()
