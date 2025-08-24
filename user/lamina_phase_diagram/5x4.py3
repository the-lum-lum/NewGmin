import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# Set path and filenames
base_path = r"C:\Users\callu\Documents\NewGmin\user\lamina_phase_diagram\#memory'"
num_plots = 20
rows, cols = 5, 4

# Set figure size (adjust as needed)
fig, axes = plt.subplots(rows, cols, figsize=(12, 15))

# Flatten axes for easier iteration
axes = axes.flatten()

# Loop through plots
for i in range(num_plots):
    img_path = os.path.join(base_path, f"Iteration_{i+1}_plot.png")
    img = mpimg.imread(img_path)
    axes[i].imshow(img)
    axes[i].axis('off')  # Turn off axis
    # Optional: remove borders
    for spine in axes[i].spines.values():
        spine.set_visible(False)

# Remove spacing between subplots
plt.subplots_adjust(wspace=0, hspace=0)

# Save the whole figure
save_path = os.path.join(base_path, "combined_5x4_grid.png")
plt.savefig(save_path, bbox_inches='tight', pad_inches=0, dpi=300)
plt.close()

print(f"Saved 5x4 grid to {save_path}")
