import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === USER SETTINGS (match your other script) ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 spike beta0.01, delta 2e-5, k=2"
iteration = 50  # ⬅️ Change this for other iterations
filename = "b_array.in"

# === BUILD FILE PATH ===
iteration_folder = f"Iteration_{iteration}"
b_array_path = os.path.join(base_path, iteration_folder, filename)

# === LOAD AND PLOT ===
if not os.path.exists(b_array_path):
    print(f"Could not find file: {b_array_path}")
else:
    b_array = np.loadtxt(b_array_path)

    # Reshape check (optional)
    if b_array.ndim == 1:
        raise ValueError("Loaded b_array is 1D — reshape manually if needed.")

    # Plot heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(
    b_array,
    cmap='viridis',
    vmin=0.001,  # lower bound of stiffness scale (B_min)
    vmax=0.01,   # upper bound of initial stiffness (B_0)
    cbar_kws={'label': 'Bending Stiffness $B_{ij}$'}
)
    plt.title(f"Bending Stiffness Heatmap at Iteration {iteration}")
    plt.xlabel("Segment Index $j$")
    plt.ylabel("Filament Index $i$")
    plt.tight_layout()

    # Optional: save to file
    save_path = os.path.join(base_path, f"Iteration_{iteration}_heatmap.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved heatmap to {save_path}")

    plt.show()
