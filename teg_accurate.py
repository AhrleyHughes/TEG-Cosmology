import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

# --- 1. COSMOLOGICAL PARAMETERS (Planck 2018) ---
params = {
    'h': 0.674,
    'Om0': 0.315,
    'Ob0': 0.049,
    'ns': 0.965,
    'sigma8_lcdm': 0.811,
    'T_cmb': 2.725
}

# --- 2. TEG THEORETICAL PARAMETERS ---
# Derived from QED and Topology - NOT fitted
ALPHA_QED = 1.0 / 137.036
KAPPA_TEG = ALPHA_QED / (2 * np.pi)  # ~ 0.00116
GAMMA_TEG = 5.0 / 3.0                # nu = 5/3 FQHE State

# --- 3. PHYSICS MODULES ---

def eisenstein_hu_transfer(k, p):
    """Approximate Transfer Function (BBKS style for dependency-free use)"""
    # This approximates the shape of the power spectrum without needing external libraries
    q = k / (p['Om0'] * p['h']**2 * np.exp(-p['Ob0'] - np.sqrt(2*p['h'])*p['Ob0']/p['Om0']))
    return np.log(1 + 2.34*q) / (2.34*q) * (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**-0.25

def linear_power_spectrum(k, p):
    """ P_lin(k) = A * k^ns * T(k)^2 """
    Tk = eisenstein_hu_transfer(k, p)
    return k**p['ns'] * Tk**2

def get_sigma8_normalization(p):
    """Finds amplitude A to match Planck sigma8"""
    R = 8.0 / p['h'] 
    def integrand(k):
        x = k * R
        W = 3 * (np.sin(x) - x*np.cos(x)) / x**3
        return linear_power_spectrum(k, p) * W**2 * k**2
    integral, _ = quad(integrand, 1e-4, 100.0)
    return p['sigma8_lcdm']**2 / ((1.0 / (2 * np.pi**2)) * integral)

# --- 4. TEG SUPPRESSION LOGIC ---

def teg_suppression_factor(k, kappa):
    """D_TEG(k) = 1 - delta_max * f(k)"""
    # Linear scaling from sensitivity analysis
    delta_max = 75.0 * kappa 
    
    # Sigmoid activation at non-linear scale (k > 0.1)
    k_pivot = 0.1 
    log_k = np.log10(k)
    log_p = np.log10(k_pivot)
    
    # Smooth step function
    fk = 1.0 / (1.0 + np.exp(-2.0 * (log_k - log_p)/0.4))
    
    return 1.0 - delta_max * fk

# --- 5. MAIN CALCULATION ---

def main():
    print(f"--- TEG CALCULATION STARTED ---")
    print(f"Theory Kappa (QED): {KAPPA_TEG:.6f}")
    
    k = np.logspace(-3, 2, 500)
    
    # 1. Normalize Linear P(k)
    A_norm = get_sigma8_normalization(params)
    P_lin = A_norm * linear_power_spectrum(k, params)
    
    # 2. Apply TEG Suppression
    suppression = teg_suppression_factor(k, KAPPA_TEG)
    P_teg = P_lin * suppression 
    
    # 3. Calculate New Sigma8
    R = 8.0 / params['h']
    def teg_integrand_func(k_val):
        p_val = np.interp(k_val, k, P_teg)
        x = k_val * R
        W = 3 * (np.sin(x) - x*np.cos(x)) / (x**3 + 1e-10)
        return p_val * W**2 * k_val**2
        
    integral_teg, _ = quad(teg_integrand_func, 1e-4, 50.0, limit=100)
    sigma8_teg = np.sqrt(integral_teg / (2 * np.pi**2))
    
    suppression_pct = (1 - sigma8_teg/params['sigma8_lcdm']) * 100
    
    print(f"LambdaCDM Sigma8: {params['sigma8_lcdm']:.3f}")
    print(f"TEG Sigma8:       {sigma8_teg:.3f}")
    print(f"Suppression:      {suppression_pct:.2f}%")
    
    # --- 6. PLOTTING FIGURES ---
    plt.figure(figsize=(10, 5))
    
    # FIG 1: Power Spectrum Ratio
    plt.subplot(1, 2, 1)
    ratio = P_teg / P_lin
    plt.semilogx(k, ratio, 'r-', linewidth=3, label=f'TEG ($\\kappa={KAPPA_TEG:.5f}$)')
    plt.axhline(1, color='k', linestyle='--')
    plt.axvspan(0.1, 10, color='yellow', alpha=0.2, label='Lensing Sensitive')
    plt.xlim(0.01, 10)
    plt.ylim(0.85, 1.02)
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('Ratio $P_{TEG} / P_{\\Lambda CDM}$')
    plt.title('Figure 1: Power Suppression')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # FIG 2: Cusp-Core (Analytic Approximation)
    plt.subplot(1, 2, 2)
    M = np.logspace(10, 15, 50)
    # Standard Duffy et al c(M)
    c_lcdm = 6.71 * (M / 2e12)**(-0.091)
    
    # TEG Modification: Virial expansion factor 85 * kappa
    virial_factor = 1.0 + 85.0 * KAPPA_TEG * (c_lcdm/5.0)**1.5
    c_teg = c_lcdm / virial_factor
    
    plt.loglog(M, c_lcdm, 'k--', label='LambdaCDM')
    plt.loglog(M, c_teg, 'b-', linewidth=3, label='TEG Prediction')
    plt.xlabel('Halo Mass $M_{200} [M_{\\odot}]$')
    plt.ylabel('Concentration c')
    plt.title('Figure 2: Cusp-Core Solution')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('TEG_Figures_1_and_2.png', dpi=150)
    print("Figures saved to TEG_Figures_1_and_2.png")
    
    # Save Data
    np.savez('teg_data.npz', k=k, ratio=ratio, M=M, c_lcdm=c_lcdm, c_teg=c_teg)

if __name__ == "__main__":
    main()
