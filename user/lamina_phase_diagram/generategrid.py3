import numpy as np
import matplotlib.pyplot as plt
import os

# User inputs and settings
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 memory recov beta0.01, k=3"
N_ROD = 10
N_SEG = 200
L = 5
nheader = 4
nheaderinit = 0

# Which iterations you want to plot (choose manually)
chosen_iterations = [1,3,5,7,9,11,13,15,17,19,20,24,28,32,36,38,40,41,43,47,50]

# Grid layout settings
rows = 7
cols = 3  # adjust depending how many you have

# Initialize figure
fig, axes = plt.subplots(rows, cols, figsize=(15, 10))
axes = axes.flatten()

# Calculate segment length
LS = L / N_SEG

for idx, iter in enumerate(chosen_iterations):
    if idx >= len(axes):
        break  # Prevent index error if more iterations than axes
    
    iteration_folder = f"Iteration_{iter}"
    coords_file = os.path.join(base_path, iteration_folder, "lowest")
    initpos_file = os.path.join(base_path, iteration_folder, "initpos")
    
    # Check if files exist
    if not os.path.isfile(coords_file) or not os.path.isfile(initpos_file):
        print(f"Files not found for {iteration_folder}. Skipping iteration {iter}.")
        continue
    
    # Load data
    with open(coords_file, 'r') as fid:
        angles = np.loadtxt(fid, skiprows=nheader)
    with open(initpos_file, 'r') as fid:
        inits = np.loadtxt(fid, skiprows=nheaderinit)
    
    print(f"Iteration {iter}: {len(angles)} angles loaded.")
    
    # Build coordinates
    coords = np.zeros([N_ROD, N_SEG, 2])
    coords[:, 0, 0] = inits[:, 0]
    coords[:, 0, 1] = inits[:, 1]
    s = 0
    for i in range(N_ROD):
        for j in range(1, N_SEG):
            coords[i, j, 0] = coords[i, j-1, 0] + LS * np.cos(angles[s])
            coords[i, j, 1] = coords[i, j-1, 1] + LS * np.sin(angles[s])
            s += 1

    # Plot
    ax = axes[idx]
    for i in range(N_ROD):
        ax.plot(coords[i, :, 0], coords[i, :, 1])
    ax.set_title(f"Iter {iter}", fontsize=8)
    ax.axis('equal')
    ax.axis('off')

# Remove any empty subplots
for k in range(len(chosen_iterations), len(axes)):
    fig.delaxes(axes[k])

# Adjust layout
plt.tight_layout()

# Save
save_path = os.path.join(base_path, "chosen_grid_plot.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0)
print(f"Saved selected grid plot to {save_path}")

# Show
plt.show()
