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

# ----------------------------- core update -------------------------------------

def compute_damage_with_memory_network(
    angles: np.ndarray,
    prev_stiffness: np.ndarray,
    memory_matrix: np.ndarray,
    *,
    k: int = 3,               # memory length required to trigger recovery
    alpha: float = 0.05,      # damage factor (linear part)
    alpha_nc: float = 1.0,    # non‑linear correction coefficient
    delta: float = 1e-4,      # constant degradation per step
    beta: float = 1e-3,       # local recovery magnitude
    theta_rec: float = np.pi/6,
    B_max: float = 0.015,
    max_angle: float = np.pi/4,
    # --- network‑coupled parameters ---
    lambda_net: float = 2e-4, # coupling strength λ
    epsilon_cross: float = 0.4, # attenuation across different filaments ε
    iteration=None
):
    """Return updated *b_array* with network‑coupled recovery.

    Parameters
    ----------
    angles : ndarray (n_filaments, n_seg−1)
        Local bend angles |θ| from previous iteration.
    prev_stiffness : ndarray (same shape)
        Bending stiffness from previous iteration.
    memory_matrix : ndarray (same shape)
        Counts how many consecutive steps each node had |θ| < theta_rec.
    k : int
        Recovery requires exactly *k* consecutive below‑threshold steps.
    ... (see defaults above) ...
    """
    n_fil, n_seg = angles.shape

    # ---------------- ρ‑matrix (local recovery weight) -------------------------
    rho = np.zeros_like(prev_stiffness)
    for i in range(n_fil):
        for j in range(n_seg):
            if memory_matrix[i, j] == k:
                rho[i, j] = beta * (B_max - prev_stiffness[i, j]) * \
                             np.exp(-abs(angles[i, j]) / theta_rec)

    # ---------------- stiffness update ----------------------------------------
    new_B = np.zeros_like(prev_stiffness)
    # Tracking matrix for where spikes occurred
    spike_mask = np.zeros_like(prev_stiffness, dtype=bool)
    rng = np.random.default_rng()  # moved outside loop for consistency


    total_rec_local = 0.0
    total_rec_net = 0.0

    for i in range(n_fil):
        for j in range(n_seg):
            theta = abs(angles[i, j])
            gamma = np.clip(theta / max_angle, 0.0, 1.0)
            D_loc = np.clip(alpha * gamma, 0.0, 1.0)

           # ---------------- Random spike addition ----------------
            p_spike = 0.01         # 1% chance of random damage per segment per iteration
            D_spike = 0.2          # Spike magnitude added to D_loc (tune as needed)
            rng = np.random.default_rng()

            if rng.random() < p_spike:
                D_loc += D_spike
                spike_mask[i, j] = True
                # Optional: log spike
                # print(f"Random spike at ({i},{j}) — D_loc = {D_loc:.3f}")

            D_loc = np.clip(D_loc, 0.0, 1.0)  # ensure valid range



            # ξ = damage with non‑linear denominator
            B_temp = (prev_stiffness[i, j] * (1.0 - D_loc)) / (1.0 + alpha_nc * D_loc)

            # local recovery contribution
            rec_local = rho[i, j]
            total_rec_local += rec_local

            # network‑coupled recovery
            rec_net = 0.0
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < n_fil and 0 <= nj < n_seg:
                        weight = rho[ni, nj]
                        if ni != i:  # hopping to a different filament ⇒ attenuate
                            weight *= epsilon_cross
                        diff = prev_stiffness[ni, nj] - prev_stiffness[i, j]
                        if diff > -1e-4:
                            rec_net += lambda_net * weight * diff
            
            total_rec_net += rec_net

            B_new = B_temp - delta + rec_local + rec_net
            new_B[i, j] = np.clip(B_new, 1e-3, B_max)
    total_spikes = np.sum(spike_mask)
    print(f"[random damage] total spikes this iteration: {total_spikes}")
    print(f"[recovery] Total local recovery added: {total_rec_local:.6f}")
    print(f"[recovery] Total network-coupled recovery added: {total_rec_net:.6f}")
    if iteration is not None:
        out_path = os.path.join(f"Iteration_{iteration}", "spike_mask.txt")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        np.savetxt(out_path, spike_mask.astype(int), fmt="%d")

    return new_B

# ------------------------------ driver -----------------------------------------

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
                k=memory_length,iteration=iteration,
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
        alpha             = 0.05,
        alpha_nc          = 1.0,
        delta             = 2e-5,
        beta              = 0.01,
        theta_rec         = np.pi/6,
        B_max             = 0.015,
        max_angle         = np.pi/4,
        memory_length     = 1,
        lambda_net        = 2e-4,
        epsilon_cross     = 0.4,
    )
