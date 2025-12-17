import numpy as np
from scipy.integrate import simps

# --- Constants & Setup ---
h = 0.674
Om0 = 0.315
Ob0 = 0.049
ns = 0.965
sigma8_target_planck = 0.811

# --- Core Physics Functions ---
def eisenstein_hu_transfer(k, h, Om0, Ob0):
    """Eisenstein & Hu (1998) transfer function (no-wiggles)."""
    theta_cmb = 2.7255 / 2.7
    s = 44.5 * np.log(9.83 / (Om0 * h**2)) / np.sqrt(1 + 10 * (Ob0 * h**2)**0.75)
    alpha_gamma = 1 - 0.328 * np.log(431 * Om0 * h**2) * (Ob0/Om0) + \
                  0.38 * np.log(22.3 * Om0 * h**2) * (Ob0/Om0)**2
    gamma_eff = Om0 * h * (alpha_gamma + (1 - alpha_gamma) / (1 + (0.43 * k * s)**4))
    q = k / (gamma_eff) * theta_cmb**2
    L = np.log(2 * np.e + 1.8 * q)
    C = 14.2 + 731.0 / (1 + 62.5 * q)
    return L / (L + C * q**2)

def top_hat_window(k, R):
    """Fourier transform of spherical top-hat window."""
    x = k * R
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3

def calculate_sigma8(k, Pk, R=8.0):
    """Compute sigma8 via integration."""
    W = top_hat_window(k, R)
    integrand = Pk * W**2 * k**2 / (2 * np.pi**2)
    return np.sqrt(simps(integrand * k, np.log(k))) # Integrate in log space

def teg_suppression(k, kappa, Gamma):
    """
    TEG suppression factor D(k). 
    P_TEG = P_LCDM * D(k)^2
    """
    k_nl = 0.13
    # Max suppression scales with kappa. Baseline 0.088 at kappa=0.0012
    max_suppression = 0.088 * (kappa / 0.0012)
    sharpness = 1.5 * Gamma
    # The factor D(k)
    D_k = 1.0 - max_suppression / (1.0 + (k_nl / k)**sharpness)
    return D_k

# --- Sensitivity Analysis Run ---
def run_sensitivity_analysis():
    # K-space array
    k = np.logspace(-4, 2, 1000)
    
    # 1. Generate Baseline LCDM Spectrum
    Tk = eisenstein_hu_transfer(k, h, Om0, Ob0)
    Pk_raw = k**ns * Tk**2
    
    # Normalize to Planck sigma8
    s8_raw = calculate_sigma8(k, Pk_raw)
    norm = (sigma8_target_planck / s8_raw)**2
    Pk_lcdm = Pk_raw * norm
    
    print(f"{'Kappa':<10} | {'Gamma':<10} | {'Sigma8_TEG':<12} | {'Suppression (%)':<15}")
    print("-" * 55)
    
    # 2. Parameter Sweep
    kappas = [0.0005, 0.0008, 0.0010, 0.0012, 0.0015, 0.0020]
    gamma_fixed = 5.0/3.0
    
    for kap in kappas:
        # Calculate Suppression Factor D(k)
        D_k = teg_suppression(k, kap, gamma_fixed)
        
        # Apply TEG: P_TEG = P_LCDM * D(k)^2
        Pk_teg = Pk_lcdm * (D_k**2)
        
        # Calculate Sigma8
        s8_teg = calculate_sigma8(k, Pk_teg)
        suppression_pct = (1 - s8_teg / sigma8_target_planck) * 100
        
        print(f"{kap:<10.4f} | {gamma_fixed:<10.2f} | {s8_teg:<12.3f} | {suppression_pct:<15.2f}")

    print("\n--- Gamma Sensitivity (at Kappa=0.0012) ---")
    gammas = [1.0, 1.33, 1.67, 2.0]
    kappa_fixed = 0.0012
    
    for gam in gammas:
        D_k = teg_suppression(k, kappa_fixed, gam)
        Pk_teg = Pk_lcdm * (D_k**2)
        s8_teg = calculate_sigma8(k, Pk_teg)
        suppression_pct = (1 - s8_teg / sigma8_target_planck) * 100
        print(f"{kappa_fixed:<10.4f} | {gam:<10.2f} | {s8_teg:<12.3f} | {suppression_pct:<15.2f}")

if __name__ == "__main__":
    run_sensitivity_analysis()
