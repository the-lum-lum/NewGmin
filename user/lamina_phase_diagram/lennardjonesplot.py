import numpy as np
import matplotlib.pyplot as plt

# Constants (you can tweak these)
C1 = 1.0
C2 = 5.0

# Normalised area A/A0 from near 0 to 2
A_over_A0 = np.linspace(0.01, 2, 500)

# Repulsion energy
E_rep = C1 * np.exp(-C2 * A_over_A0)

# Plot
plt.figure(figsize=(7, 5))
plt.plot(A_over_A0, E_rep, label=r'$E_{\mathrm{rep}} = C_1 \exp\left(-C_2 \frac{A}{A_0} \right)$', color='darkorange')

plt.title("Exponential Area-Based Repulsion Energy", fontsize=14)
plt.xlabel(r'Normalised Area $A / A_0$', fontsize=12)
plt.ylabel(r'Energy $E_{\mathrm{rep}}$', fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
