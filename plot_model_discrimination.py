import numpy as np
import matplotlib.pyplot as plt

k = np.logspace(-2, 1.5, 200)

# 1. TEG Prediction (z=0)
# Matches your paper: monotonic suppression ~6%
teg_ratio = 1.0 - 0.058 * (1.0 / (1.0 + np.exp(-2.0 * (np.log10(k) + 0.5))))

# 2. Typical Baryonic Feedback (Schematic of AGN feedback like HMCode/OWLS-AGN)
# Baryons often "overshoot" - suppression is deeper but localized
# This is a schematic representation of the "spoon" shape
baryon_ratio = 1.0 - 0.15 * np.exp(-0.5 * (np.log10(k) - 0.5)**2) 
# Add the "upturn" at high k (cooling)
baryon_ratio += 0.05 * (1.0 / (1.0 + np.exp(-3.0 * (np.log10(k) - 1.0))))

plt.figure(figsize=(8, 6))

# Plot lines
plt.semilogx(k, teg_ratio, color='blue', linewidth=3, label='TEG (Entropic Screening)')
plt.semilogx(k, baryon_ratio, color='gray', linestyle='--', linewidth=2, label='Typical Baryonic Feedback (AGN)')
plt.axhline(1, color='k', linestyle='-', linewidth=1)

# Annotations
plt.text(0.2, 0.96, "Linear Scales\n(Agrees)", ha='center', fontsize=9, color='green')
plt.text(3.0, 0.88, "Distinct Shape", ha='center', fontsize=10, color='red', fontweight='bold')

# Formatting
plt.xlabel(r'Wavenumber $k$ [$h$/Mpc]', fontsize=12)
plt.ylabel(r'Ratio $P(k) / P_{\Lambda\mathrm{CDM}}$', fontsize=12)
plt.title('Model Discrimination: TEG vs. Baryonic Feedback', fontsize=14)
plt.legend(loc='lower left', fontsize=11)
plt.grid(True, alpha=0.2)
plt.ylim(0.8, 1.05)
plt.xlim(0.01, 20)

plt.tight_layout()
plt.savefig('TEG_vs_Baryons.png', dpi=300)
print("Discrimination figure generated.")
