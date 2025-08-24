import numpy as np
import matplotlib.pyplot as plt
import os

# === USER SETTINGS ===
iteration = 30  # Change this to any iteration you want to visualise
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#damage spikes"  # Your simulation folder
file_path = os.path.join(base_path, f"Iteration_{iteration}", "spike_mask.txt")

# === Load and Plot ===
if os.path.exists(file_path):
    spike_mask = np.loadtxt(file_path)
    plt.figure(figsize=(10, 6))
    plt.imshow(spike_mask, cmap='Reds', interpolation='none', aspect='auto')
    plt.title(f"Random Damage Spike Map — Iteration {iteration}")
    plt.xlabel("Segment Index $j$")
    plt.ylabel("Filament Index $i$")
    plt.colorbar(label="Spike (1 = Yes)")
    plt.tight_layout()
    plt.show()
else:
    print(f"No spike_mask.txt found for iteration {iteration}")
