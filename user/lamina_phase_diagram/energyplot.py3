import os
import numpy as np
import matplotlib.pyplot as plt

# === USER SETTINGS ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 spike beta0.01, delta 2e-5, k=2"
filename = "GMIN_OUT"
num_iterations = 50  # Adjust as needed

# === Extract final energy per iteration ===
energies = []
iters_found = []

for i in range(1, num_iterations + 1):
    folder = f"Iteration_{i}"
    file_path = os.path.join(base_path, folder, filename)

    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            lines = [line for line in f if "Energy" in line]
            if lines:
                # Take energy from the last line containing "Energy"
                last_line = lines[-1]
                try:
                    energy_str = last_line.split()[1]
                    energy_val = float(energy_str)
                    energies.append(energy_val)
                    iters_found.append(i)
                except Exception as e:
                    print(f"Could not parse energy in {file_path}: {e}")
    else:
        print(f"Missing file: {file_path}")

# === Plot energy evolution ===
plt.figure(figsize=(8, 5))
plt.plot(iters_found, energies, marker='o', linestyle='-', linewidth=2)
plt.title("Energy Evolution Across Iterations")
plt.xlabel("Iteration")
plt.ylabel("Total Energy $\\mathcal{E}$")
plt.grid(True)
plt.tight_layout()
plt.show()
