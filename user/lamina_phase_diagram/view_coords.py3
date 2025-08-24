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

# Set the number of iterations (adjust as needed)
num_iterations = 50

# List of iterations for which to save the plot
#save_iterations = []  # adjust as needed

# Loop over each iteration folder
for iter in range(1, num_iterations+1):
    iteration_folder = f"Iteration_{iter}"
    coords_file = os.path.join(base_path, iteration_folder, "lowest")
    initpos_file = os.path.join(base_path, iteration_folder, "initpos")
    
    # Check if the files exist
    if not os.path.isfile(coords_file) or not os.path.isfile(initpos_file):
        print(f"Files not found for {iteration_folder}. Skipping iteration {iter}.")
        continue

    # Read the angles and initial positions
    with open(coords_file, 'r') as fid:
        angles = np.loadtxt(fid, skiprows=nheader)
    with open(initpos_file, 'r') as fid:
        inits = np.loadtxt(fid, skiprows=nheaderinit)
    
    print(f"Iteration {iter}: {len(angles)} angles loaded.")
    
    # Restructure the data and convert to xy coordinates
    LS = L / N_SEG
    coords = np.zeros([N_ROD, N_SEG, 2])
    coords[:, 0, 0] = inits[:, 0]
    coords[:, 0, 1] = inits[:, 1]
    s = 0
    for i in range(N_ROD):
        for j in range(1, N_SEG):
            coords[i, j, 0] = coords[i, j-1, 0] + LS * np.cos(angles[s])
            coords[i, j, 1] = coords[i, j-1, 1] + LS * np.sin(angles[s])
            s += 1

    # Create the plot for the current iteration
    plt.figure()
    #plt.title(f"Iteration {iter}")
    for i in range(N_ROD):
        plt.plot(coords[i, :, 0], coords[i, :, 1])
    #plt.xlabel("X")
    #plt.ylabel("Y")
    #plt.axis('equal')
    #plt.tight_layout()
    plt.axis('equal')
    plt.axis('off')  # Turn off axes completely

    # Save the figure if the iteration is in our save list
    #if iter in save_iterations:
    save_path = os.path.join(base_path, f"Iteration_{iter}_plot.png")
    #plt.savefig(save_path)
    plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
    print(f"Saved plot for iteration {iter} to {save_path}")

    # Display the plot briefly (or wait for a key press)
    plt.show(block=False)
    plt.pause(1.5)  # Pause for 1.5 seconds; adjust as needed
    plt.close()
