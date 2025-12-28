# Topological Entropic Gravity (TEG)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18051561.svg)](https://doi.org/10.5281/zenodo.18051561)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains code and supporting material for:

**Hughes, A. (2025)**  
*Topological Entropic Gravity: Unifying the Quantum Hall Vacuum with Cosmic Structure Formation to Resolve the S₈ Tension*  
Preprint: https://doi.org/10.5281/zenodo.18051561

---

## Scope and Motivation

Weak lensing surveys (KiDS, DES, HSC) consistently infer a lower clustering amplitude than Planck CMB constraints under ΛCDM, producing the well-known **S₈ tension** at the ~3σ level. Empirically, resolving the tension requires a **5–7% suppression of the matter power spectrum** at non-linear scales while preserving linear-scale growth.

Baryonic feedback models can produce suppression, but generally require strong, scale-localized AGN prescriptions and often introduce redshift-dependent behavior in conflict with CMB and Lyman-α constraints.

**Topological Entropic Gravity (TEG)** introduces an alternative mechanism: a scale-dependent, non-dissipative pressure sourced by the topology of the vacuum itself. The framework is intentionally minimal and designed to be falsifiable at the level of the matter power spectrum and halo structure.

---
 ![TEG-HUGHES-2025-visual_abstract](TEG-HUGHES-2025-visual_abstract.png)
## Core Assumptions

1. The vacuum admits an effective description as an incompressible topological fluid.
2. The equation of state is stiff, with adiabatic index Γ = 5/3, motivated by stability of the ν = 5/3 fractional quantum Hall state.
3. Gravitational compression locally increases the effective filling factor, inducing a repulsive **entropic pressure**.
4. The coupling strength κ is small and fixed by fundamental physics, not fitted to galaxy data.

---

## Parameter Determination

TEG introduces a single dimensionless coupling κ. Importantly, the value required to resolve the S₈ tension coincides with independent theoretical estimates.

| Origin | Estimate | κ |
|--------|----------|---|
| Vacuum polarization (QED) | κ = α / 2π | 0.00116 |
| Thermodynamic scaling (BBN) | κ ~ η⁻¹ᐟ³ | 0.00085 |
| Required suppression (S₈) | — | 0.00120 |

Agreement is within O(10%), consistent with expected uncertainties in early-universe thermodynamic estimates.

No additional parameters are tuned.

---

## Observable Consequences

### 1. Matter Power Spectrum

- Suppression activates only for k ≳ 0.1 h/Mpc
- Linear scales remain unchanged
- Net suppression at z ≈ 0 is ~5–6%, consistent with weak lensing data

### 2. Halo Structure

- Reduced virial overdensity modifies halo concentrations
- For M ≲ 10¹¹ M☉, NFW-like cusps are softened into approximately constant-density cores
- Effect arises from modified collapse dynamics, not feedback energy injection

---

## Theoretical Formalism (Sketch)

Spherical collapse with entropic pressure:

```
d²R/dt² = -GM/R² - (1/ρ)∇P_ent + ΛR/3
```

with

```
P_ent = κ c² ρ_crit (ρ/ρ̄)^(5/3)
```

At turnaround, the pressure term reduces the maximum overdensity.  
At virialization, the modified virial theorem

```
2T + W - 3∫P dV = 0
```

leads to a larger equilibrium radius for fixed mass, reducing concentration.

This mechanism operates without violating linear growth constraints.

---

## Code Structure

### Installation

```bash
git clone https://github.com/AhrleyHughes/TEG-Cosmology.git
cd TEG-Cosmology
pip install -r requirements.txt
```

### Main Calculation

```bash
python teg_accurate.py
```

Outputs:

* Matter power spectrum ratios
* σ₈ comparison with Planck and lensing values
* Saved numerical arrays for external analysis

### ΛCDM Recovery Test
```bash
python test_lcdm_limit.py
```

Critical credibility check: Demonstrates that TEG cleanly reduces to ΛCDM when κ → 0.

![ΛCDM Recovery Test](TEG_LCDM_Recovery.png)

**Figure:** Exact recovery of ΛCDM in the κ → 0 limit. P(k) and c(M) ratios deviate from unity by < 10⁻⁸%, confirming the absence of numerical artifacts or hidden assumptions. The model reduces cleanly to standard cosmology when the topological coupling is disabled.

### Sensitivity Analysis

```bash
python sensitivity_analysis.py
```

Explores robustness against κ variation and verifies absence of fine-tuning.

---

## Additional Diagnostics

### Redshift Evolution

Script: `plot_redshift_evolution.py`

![Redshift Evolution](TEG_Redshift_Evolution.png)

**Figure: Predicted evolution of suppression.** The effect activates at late times (z < 1), consistent with thawing dark energy constraints.

### Model Discrimination

Script: `plot_model_discrimination.py`

![Model Discrimination](TEG_vs_Baryons.png)

**Figure: TEG produces a stable 'shelf' of suppression, distinct from the characteristic 'spoon' shape of AGN feedback, offering a clear discriminant for future observations.**

---
## Robustness Tests

TEG is intentionally constrained. These tests demonstrate that its predictions are
**monotonic, shape-invariant, and non-degenerate** with baryonic feedback models.

### 1. Merciless Parameter Sweep

The response of σ₈ to the coupling constant κ is smooth and monotonic.
There are no instabilities, jumps, or tuning artifacts.

![σ8 sweep](kappa_sweep.png)

Key result:
- Increasing κ produces controlled suppression
- The weak-lensing target range is reached naturally
- No parameter fine-tuning is required

---
## Distinguishing TEG from Baryonic Feedback

A common objection to modifications of the matter power spectrum on non-linear scales is the degeneracy with baryonic feedback processes (e.g., AGN feedback, supernovae winds). However, TEG is phenomenologically distinct from baryonic mechanisms in three testable ways:

1. **Derivative Signature (The "Fingerprint"):** Baryonic feedback models (such as those in BAHAMAS or IllustrisTNG) inject energy at specific scales, creating a characteristic 'spoon' feature in the power spectrum ratio—a suppression followed by a high-$k$ rebound or inflection. TEG, by contrast, produces a monotonic, scale-dependent suppression with a smooth gradient (see Figure `TEG_vs_Baryons_Derivative.png`). There is no 'rebound' because there is no energy injection, only a modification of the collapse threshold.

2. **Shape Invariance:** The TEG suppression profile is structurally rigid. Varying the coupling $\kappa$ changes the amplitude but preserves the functional form (see Figure `shape_invariance.png`). Baryonic models are effective field theories with multiple free parameters that can arbitrarily alter the shape of the suppression. TEG cannot be 'tuned' to mimic arbitrary shapes; if the data shows a rebound, TEG is falsified.

3. **Redshift Dependence:** Baryonic effects are cumulative and strongly redshift-dependent, tracking the history of star formation and AGN activity. TEG suppression is tied to the virialization of structure and the effective density contrast, predicting a late-time activation that distinctively tracks the growth factor, independent of local thermodynamic history.

---

### 2. Shape Invariance Test

Changing κ alters only the **amplitude**, not the **shape**, of suppression.

![Shape invariance](shape_invariance.png)

All curves collapse to a single universal profile when normalized.
This distinguishes TEG from baryonic feedback models, which freely alter shape.

---

### 3. Model Discrimination via Derivatives

Even when power-spectrum ratios appear similar, their derivatives reveal
distinct physical mechanisms.

![Derivative comparison](TEG_vs_Baryons_Derivative.png)

TEG exhibits:
- Smooth
- Monotonic
- Single-scale behavior

Baryonic feedback exhibits:
- Inflection points
- Rebound features
- Multi-scale tuning

---
## Falsifiability and Future Tests

### Robust, Near-Term Predictions
The TEG framework makes concrete, falsifiable predictions that stem directly from the modified power spectrum logic. These are not fitted parameters; they are structural consequences of the theory.

* **Scale-Dependent Suppression:** A specific reduction in power on non-linear scales with **zero** impact on linear (CMB) scales (see `kappa_sweep.png`).
* **Shape Invariance:** The suppression profile maintains a universal shape across different coupling strengths, distinguishing it from arbitrary curve-fitting (see `shape_invariance.png`).
* **Monotonic High-k Behavior:** Unlike baryonic feedback models (AGN), TEG predicts no high-$k$ power recovery or "rebound" (see `TEG_vs_Baryons_Derivative.png`).

These signatures are testable *now* with current weak-lensing data.

### Optional Phenomenological Extensions (Distinct from Core Theory)
The core mechanism resolves the S₈ tension without these additions. However, if the TEG coupling arises from deeper microphysical structure, secondary effects may appear:

* **Mild Halo Oblateness:** Potential signature from chiral edge modes, detectable by Gaia/Euclid.
* **Temperature-Dependent κ Scaling:** Possible variations in high-z lensing observable by JWST.
* **Discrete κ Transitions:** Theoretical possibility of vacuum topology phase changes.

These effects are **avenues for investigation**, not load-bearing pillars of the primary result.

---

## Citation

```bibtex
@article{Hughes2025TEG,
  title   = {Topological Entropic Gravity: Unifying the Quantum Hall Vacuum with Cosmic Structure Formation to Resolve the S8 Tension},
  author  = {Hughes, Ahrley},
  year    = {2025},
  doi     = {10.5281/zenodo.18051561},
  url     = {https://doi.org/10.5281/zenodo.18051561}
}
```

---

## Contact

Ahrley Hughes  
Independent Researcher  
arhleyhughes@proton.me

---

## License

MIT License — see `LICENSE`.
