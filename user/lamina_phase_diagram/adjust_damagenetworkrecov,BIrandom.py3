
import sys
import os
import numpy as np
from numpy.random import default_rng 
# -----------------------------------------------------------------------------
# adjust_damage_network.py  --  bending‑stiffness update with
#   • progressive curvature‑based damage
#   • curvature‑memory healing
#   • non‑linear neutral‑axis shift (Van Paepegem‑style)
#   • NEW: network‑coupled recovery that diffuses between neighbouring nodes
# -----------------------------------------------------------------------------
# 2025‑04‑19  –  draft integration for Filament‑Bundle Energy‑Minimisation project
# -----------------------------------------------------------------------------

# --------------------------- helper: history I/O --------------------------------

def update_curvature_history(folder: str, iteration: int, angles: np.ndarray, *, num_history: int = 3):
    """Save local curvature matrix for *iteration* and keep only last *num_history* files."""
    os.makedirs(folder, exist_ok=True)
    np.save(os.path.join(folder, f"curv_{iteration}.npy"), angles)

    # purge files older than the history window
    for old_iter in range(max(0, iteration - num_history)):
        old_file = os.path.join(folder, f"curv_{old_iter}.npy")
        if os.path.exists(old_file):
            os.remove(old_file)


def load_memory_matrix(folder: str, iteration: int, theta_rec: float, *, num_history: int = 3):
    """Return an int matrix whose entries give the number of times each node was below
    *theta_rec* in the last *num_history* iterations (inclusive). If history is missing
    fall back to zeros.
    """
    memory_sum = None
    for step in range(iteration - num_history + 1, iteration + 1):
        f = os.path.join(folder, f"curv_{step}.npy")
        if not os.path.exists(f):
            continue
        angles = np.load(f)
        below = (np.abs(angles) < theta_rec).astype(int)
        memory_sum = below if memory_sum is None else memory_sum + below
    return memory_sum if memory_sum is not None else np.zeros_like(angles, dtype=int)




def compute_damage_with_memory_network(
    angles: np.ndarray,
    prev_stiffness: np.ndarray,
    memory_matrix: np.ndarray,
    *,
    k: int = 3,
    alpha: float = 0.05,
    alpha_nc: float = 1.0,
    delta: float = 1e-4,
    beta: float = 1e-3,
    theta_rec: float = np.pi / 6,
    B_max: float = 0.015,
    max_angle: float = np.pi / 4,
    lambda_net: float = 1e-3,
    epsilon_cross: float = 0.4,
    seed_contrast: bool = False,
    iteration=None,
):
    """
    Update bending stiffness with:
      • curvature-based degradation
      • calm-memory local recovery
      • network-coupled recovery
      • NO random spike events
    """
    B_min=1e-3
    n_fil, n_seg = angles.shape
    rho = np.zeros_like(prev_stiffness)

    # Compute local recovery weight matrix (rho)
    calm = memory_matrix >= k
    rho[calm] = beta * (B_max - prev_stiffness[calm]) * \
                np.exp(-np.abs(angles[calm]) / theta_rec)

    new_B = np.empty_like(prev_stiffness)

    # Track local and network contributions
    local_recovery = np.zeros_like(prev_stiffness)
    network_recovery = np.zeros_like(prev_stiffness)
    tot_loc, tot_net = 0.0, 0.0

    for i in range(n_fil):
        for j in range(n_seg):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0.0, 1.0)
            D_loc = alpha * gamma

            # Curvature-based degradation
            B_tmp = (prev_stiffness[i, j] * (1.0 - D_loc)) / (1.0 + alpha_nc * D_loc)

            # Calm memory-based local recovery
            rec_local = rho[i, j]
            

            # Network-coupled healing from neighbours
            rec_net = 0.0
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < n_fil and 0 <= nj < n_seg:
                        weight = rho[ni, nj]
                        if ni != i:
                            weight *= epsilon_cross
                        diff = prev_stiffness[ni, nj] - prev_stiffness[i, j]
                        if diff > 0.0:
                            rec_net += lambda_net * weight * diff
            
            if prev_stiffness[i, j] <= B_min * 1.000001:
                rec_local = 0.0
                rec_net = 0.0

            local_recovery[i, j] = rec_local
            tot_loc += rec_local
            network_recovery[i, j] = rec_net
            tot_net += rec_net

            # Final update
            new_B[i, j] = np.clip(B_tmp - delta + rec_local + rec_net, B_min, B_max)

    print(f"[recovery] local: {tot_loc:.4e}  network: {tot_net:.4e}")

    # Optionally save matrices
    if iteration is not None:
        out_dir = f"Iteration_{iteration}"
        os.makedirs(out_dir, exist_ok=True)
        np.savetxt(os.path.join(out_dir, "local_recovery.txt"), local_recovery, fmt="%.6e")
        np.savetxt(os.path.join(out_dir, "network_recovery.txt"), network_recovery, fmt="%.6e")

    return new_B




def adjust_bending_stiffness(
    iteration: int,
    num_filaments: int,
    N_SEG: int,
    *,
    initial_stiffness: float = 1e-2,
    lower_rat: float = 1.0,
    memory_length: int = 3,
    curvature_dir: str = "curvature_history",
    **kwargs,  # pass‑through to compute_damage_with_memory_network
):
    """Pipeline wrapper that writes the file *b_array.in* for the Fortran core."""
    num_bend = N_SEG - 1

    # --------------- first iteration (uniform B) ------------------------------
    if iteration == 1:
        # --- truncated‑normal initial stiffness -----------------------------
        mean_B  = initial_stiffness           # 0.01
        spread  = 0.15                        # ±10 %
        sigma   = 0.05 * mean_B              # 3 % s.d. before truncation
        low, hi = mean_B * (1 - spread), mean_B * (1 + spread)

        rng = default_rng()                  # remove argument for full randomness
        draws = rng.normal(loc=mean_B, scale=sigma, size=(num_filaments, num_bend))
        random_base = np.clip(draws, low, hi)   # hard truncation

        

        b_array = random_base.copy()
    else:
        # ----- load angles -----
        last_folder = f"Iteration_{iteration - 1}"
        lowest_file = os.path.join(last_folder, "lowest")
        if not os.path.isfile(lowest_file):
            print(f"Warning: {lowest_file} not found – using uniform stiffness.")
            angles = np.zeros((num_filaments, num_bend))
        else:
            angles_all = np.loadtxt(lowest_file, skiprows=4)
            angles = angles_all.reshape((num_filaments, num_bend))

        # ----- load previous B -----
        prev_file = os.path.join(last_folder, "b_array.in")
        if os.path.isfile(prev_file):
            prev_B_flat = np.loadtxt(prev_file)
            prev_B = prev_B_flat.reshape((num_filaments, num_bend))
        else:
            prev_B = np.full((num_filaments, num_bend), initial_stiffness)

        # ----- curvature history + memory matrix -----
        update_curvature_history(curvature_dir, iteration - 1, angles, num_history=memory_length)
        memory_mat = load_memory_matrix(curvature_dir, iteration - 1, kwargs.get("theta_rec", np.pi/6),
                                         num_history=memory_length)

        # ----- compute new B -----
        b_array = compute_damage_with_memory_network(
            angles, prev_B, memory_mat,
            k=memory_length,
            seed_contrast=(iteration == 1),
            iteration=iteration,
            **kwargs)


    # ---------- write to disk ----------
    with open('b_array.in', 'w') as file:
        np.savetxt(file, b_array, fmt="%.6f")
    print(f"Updated B_ARRAY for iteration {iteration}:")
    print(b_array)
    


# ----------------------------- CLI --------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 adjust_damage_network.py <iteration> <num_filaments> <num_segments>")
        sys.exit(1)

    iter_no = int(sys.argv[1])
    n_fil   = int(sys.argv[2])
    n_seg   = int(sys.argv[3])

    adjust_bending_stiffness(
        iter_no,
        n_fil,
        n_seg,
        # --- tweakables ------------------------------------------------------
        initial_stiffness = 1e-2,
        lower_rat         = 1.0,
        alpha             = 0.04,
        alpha_nc          = 1.0,
        delta             = 5e-5,
        beta              = 0.015,
        theta_rec         = np.pi/6,
        B_max             = 0.015,
        max_angle         = np.pi/4,
        memory_length     = 3,
        lambda_net        = 5e-3,
        epsilon_cross     = 0.2,
    )
