import os
import numpy as np
import pandas as pd

# === USER SETTINGS ===
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#FI_0.1_FS_0.00005 spike beta0.01, delta 2e-5, k=2"
filename = "b_array.in"
num_iterations = 50
B_min = 0.001  # Fully damaged threshold

# === Collect stats per iteration ===
data = {
    "Iteration": [],
    "Mean_B": [],
    "Std_B": [],
    "Min_B": [],
    "# Fully Damaged": []
}

for i in range(1, num_iterations + 1):
    folder = f"Iteration_{i}"
    file_path = os.path.join(base_path, folder, filename)

    if os.path.exists(file_path):
        b_array = np.loadtxt(file_path)
        if b_array.ndim == 1:
            print(f"Warning: {file_path} is 1D — skipping.")
            continue

        data["Iteration"].append(i)
        data["Mean_B"].append(np.mean(b_array))
        data["Std_B"].append(np.std(b_array))
        data["Min_B"].append(np.min(b_array))
        data["# Fully Damaged"].append(np.sum(b_array <= B_min))
    else:
        print(f"Missing: {file_path}")

# === Create DataFrame and show ===
df = pd.DataFrame(data)
print(df)


latex_table = df.to_latex(index=False, float_format="%.4f", caption="Bending stiffness statistics over iterations.", label="tab:bending_stats")
print(latex_table)
# === Optionally save to CSV ===
#save_csv_path = os.path.join(base_path, "bending_stiffness_summary.csv")
#df.to_csv(save_csv_path, index=False)
#print(f"\nSaved summary to: {save_csv_path}")
