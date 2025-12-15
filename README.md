
markdown# TEG Cosmology: Thermodynamic Entropic Gravity

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Python implementation of Thermodynamic Entropic Gravity (TEG) — a phenomenological framework for resolving the S₈ tension through local entropic screening during structure formation.

**Paper:** *"Thermodynamic Self-Regulation of Cosmic Structure: Resolving the S₈ Tension via Local Entropic Screening"* — Ahrley Hughes (2025)

---

## Overview

The standard ΛCDM model faces two major tensions:

1. **S₈ Tension**: Weak lensing surveys measure S₈ ≈ 0.76, while Planck CMB infers S₈ ≈ 0.83 (approximately 3.5σ discrepancy)
2. **Dark Energy Evolution**: DESI DR2 favors dynamical "thawing" dark energy over a cosmological constant

TEG proposes that irreversible information processing during halo collapse generates a local repulsive pressure (P_ent), amplified by the photon-to-baryon ratio η ≈ 10⁹. This mechanism:

- Suppresses the matter power spectrum by approximately 6-8% in the non-linear regime (k > 0.1 h/Mpc)
- Resolves the S₈ tension without altering linear growth
- Predicts a steeper concentration-mass relation at low masses, addressing the Cusp-Core problem

---

## Key Results

| Observable | ΛCDM (Planck) | TEG Prediction | Weak Lensing |
|------------|---------------|----------------|--------------|
| σ₈ | 0.811 | 0.764 | ~0.76 |
| Suppression | — | 5.8% | Match |
| c(M) slope | Shallow | Steep (M < 10¹¹ M☉) | Testable |
---
![TEG Visual Abstract](visual_abstract_2x.png)
---

## Installation
```bash
git clone https://github.com/ahrleyhughes/TEG-Cosmology.git
cd TEG-Cosmology
pip install -r requirements.txt
```

## Usage
```bash
python teg_accurate.py
```

**Output:**
- Prints σ₈ values and suppression percentage to console
- Generates `TEG_Figures_1_and_2.png` reproducing paper figures
- Saves data arrays to `teg_data.npz` for further analysis

---

## Figures

The script generates two figures matching the paper:

**Figure 1: Power Spectrum Suppression**
- Shows P_TEG / P_ΛCDM ratio versus wavenumber k
- Suppression activates at k ~ 0.1 h/Mpc in the non-linear regime
- Linear scales (k < 0.1) remain unaffected, preserving CMB constraints

**Figure 2: Concentration-Mass Relation**
- Compares ΛCDM versus TEG halo concentrations
- TEG predicts "puffier" dwarf halos (M < 10¹¹ M☉)
- Falsifiable with Euclid weak lensing and Gaia kinematics

---

## Theory Summary

### Core Mechanism: Entropic Pressure
```
P_ent ≈ κ c² ρ_crit (ρ_b / ρ̄_b)^Γ
```

**Parameters:**
- κ ≈ 0.0012: Coupling constant derived from η ≈ 10⁹
- Γ = 5/3: Polytropic index (phenomenological stiffness)

**Physical origin:** The vacuum acts as a thermodynamic reservoir, resisting local entropy reduction during collapse.

### Modified Power Spectrum
```
P_TEG(k) = P_ΛCDM(k) × D²_TEG(k)
```

where D_TEG(k) is the suppression factor encoding density-dependent screening.

---

## Code Structure

Main components of `teg_accurate.py`:

1. `eisenstein_hu_transfer()` - Linear transfer function (BAO-accurate)
2. `linear_power_spectrum()` - P(k) ∝ k^n_s T(k)²
3. `calculate_sigma8()` - Integrate with top-hat filter
4. `teg_suppression()` - Core TEG modification
5. `concentration_mass_relation()` - Halo structure prediction

---

## Observational Tests

| Test | Dataset | Status |
|------|---------|--------|
| S₈ suppression | KiDS-1000, DES Y3, HSC | Consistent |
| c(M) steepening | Euclid (2025+) | Pending |
| Dwarf kinematics | Gaia DR4 | Pending |
| Thawing DE | DESI DR2 | Consistent |

---

## Citation

If you use this code, please cite:
```bibtex
@article{Hughes2025TEG,
  title={Thermodynamic Self-Regulation of Cosmic Structure: Resolving the $S_8$ Tension via Local Entropic Screening},
  author={Hughes, Ahrley},
  journal={arXiv preprint arXiv:XXXX.XXXXX},
  year={2025}
}
```

---

## Customization

Modify TEG parameters in the script:
```python
kappa = 1.2e-3    # Coupling strength
Gamma = 5.0/3.0   # Polytropic index
```

Then re-run to see impact on σ₈ and c(M).

---

## Contact

**Ahrley Hughes**  
ArhleyHughes@proton.me
Cincinnati, OH  
*Seeking arXiv endorsement in astro-ph.CO*

---

## License

MIT License - see `LICENSE` file for details.

---

## Acknowledgments

- Eisenstein & Hu (1998) for transfer function formalism
- Planck, KiDS, DES, HSC collaborations for observational data
- DESI collaboration for BAO constraints
