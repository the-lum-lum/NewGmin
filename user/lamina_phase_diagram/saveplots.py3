import numpy as np
import matplotlib.pyplot as plt
import os

# User inputs and settings
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 spike beta0.01, delta 2e-5, k=2"
N_ROD = 10
N_SEG = 200
L = 5
nheader = 4
nheaderinit = 0

# Which iterations you want to plot
chosen_iterations = [1,4,7,10,13,16,19,21,24,25,27,30,32,35,38,39,40,43,44,47,50]

# Calculate segment length
LS = L / N_SEG

for iter in chosen_iterations:
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

    # Create a figure for this iteration
    fig, ax = plt.subplots(figsize=(8, 4))  # Nice big individual figure

    for i in range(N_ROD):
        ax.plot(coords[i, :, 0], coords[i, :, 1])

    ax.set_title(f"Iteration {iter}", fontsize=20)
    ax.axis('equal')
    ax.axis('off')

    # Save
    save_path = os.path.join(base_path, f"iteration_{iter}.png")
    fig.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)  # Important to close to avoid memory overload

    print(f"Saved {save_path}")
