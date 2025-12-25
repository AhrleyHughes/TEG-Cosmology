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
- Matter power spectrum ratios
- σ₈ comparison with Planck and lensing values
- Saved numerical arrays for external analysis

### Sensitivity Analysis

```bash
python sensitivity_analysis.py
```

Explores robustness against κ variation and verifies absence of fine-tuning.

---

## Additional Diagnostics

### Redshift Evolution

Script: `plot_redshift_evolution.py`

Implements the ansatz:
```
κ(z) ∝ (1+z)^(-n)
```

Used to test whether suppression naturally activates at late times (z ≲ 1), consistent with:
- Absence of tension at recombination
- Preference for thawing dark energy in DESI DR2

### Model Discrimination

Script: `plot_model_discrimination.py`

Compares TEG suppression to schematic baryonic feedback:
- Baryons produce a localized "spoon" feature with high-k upturn
- TEG predicts a monotonic, saturating shelf

This distinction is observable with next-generation lensing surveys.

---

## Falsifiability and Future Tests

### Robust Predictions

1. Scale-dependent suppression without linear-scale impact
2. Systematic reduction in dwarf-halo concentrations
3. Absence of high-k power recovery

### Speculative Extensions

1. Mild halo oblateness from chiral edge modes (Gaia / Euclid)
2. Temperature-dependent κ scaling (JWST high-z lensing)
3. Possible discrete κ transitions if vacuum topology changes phase

Speculative elements are not required for resolving the S₈ tension.

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
