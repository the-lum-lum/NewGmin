import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === USER SETTINGS ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 network recovrandomBI"
filename = "b_array.in"
iteration_range = range(1,50,5)  # e.g., plot every 4 iterations up to 20

# === Fixed color scale across plots ===
vmin = 0.001  # B_min
vmax = 0.01   # B_0 or B_max (depending on your setup)

# === Plot grid setup ===
n_cols = 5
n_rows = int(np.ceil(len(iteration_range) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows), constrained_layout=True)

axes = axes.flatten()

for idx, iteration in enumerate(iteration_range):
    ax = axes[idx]
    iteration_folder = f"Iteration_{iteration}"
    file_path = os.path.join(base_path, iteration_folder, filename)

    if os.path.exists(file_path):
        b_array = np.loadtxt(file_path)
        if b_array.ndim == 1:
            print(f"Warning: {file_path} is 1D. Skipping.")
            continue

        sns.heatmap(b_array, ax=ax, cmap='viridis', vmin=vmin, vmax=vmax, cbar=False)
        ax.set_title(f"Iteration {iteration}", fontsize=14)
        ax.set_xlabel("Segment Index $j$", fontsize=12)
        ax.set_ylabel("Filament Index $i$", fontsize=12)
        ax.tick_params(axis='both', labelsize=10)

        # Optional: reduce number of x-tick labels
        num_xticks = 6
        xtick_locs = np.linspace(0, b_array.shape[1] - 1, num_xticks, dtype=int)
        ax.set_xticks(xtick_locs)
        ax.set_xticklabels(xtick_locs, fontsize=10)

    else:
        ax.set_visible(False)
        print(f"Missing: {file_path}")

# Add a single shared colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
norm = plt.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
#fig.colorbar(sm, ax=axes, orientation='vertical', shrink=0.9, label='Bending Stiffness $B_{ij}$')
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', shrink=0.9)
cbar.set_label('Bending Stiffness $B_{ij}$', fontsize=12, labelpad=10)
cbar.ax.tick_params(labelsize=10)

#plt.suptitle("Bending Stiffness Heatmaps Across Iterations", fontsize=16)
plt.show()
