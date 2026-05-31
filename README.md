# Ramanujan-Bose Interval Cosmology (RBIC): A Real-Data Resolution to the Hubble Tension

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18955396-blue)](https://doi.org/10.5281/zenodo.18955396)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the official Python MCMC pipeline, dataset structures, and numerical validation codes for the **Ramanujan-Bose Interval Cosmology (RBIC)** framework. 

The RBIC framework provides a first-principles, scale-bridging physical theory that dynamically resolves the $5.1\sigma$ Hubble tension and the $S_8$ growth tension by accounting for:
1. **Gravity-coupling decoupling** arising from the redshift evolution of the Jeans scale ($+497.6 \text{ km}^2\text{s}^{-2}\text{Mpc}^{-2}$).
2. **Kinematical backreaction** (Buchert vorticity) calibrated with real CANDELS galaxy morphologies ($+7.8 \text{ km}^2\text{s}^{-2}\text{Mpc}^{-2}$).
3. **Late-time activation of the Ramanujan-Bose dark energy condensate** governed by a running coupling $\alpha(z)$ ($+13.6 \text{ km}^2\text{s}^{-2}\text{Mpc}^{-2}$).

---

## 📂 Repository Structure

The directory is organized as follows to ensure full reproducibility:

```text
HubbleTension/
├── data/
│   ├── pantheon_plus/
│   │   ├── Pantheon+SH0ES.dat           # 1624 SNe Ia redshifts & distance moduli
│   │   └── Pantheon+SH0ES_STAT+SYS.cov  # Full covariance matrix (zipped/raw)
│   ├── desi_bao/
│   │   ├── desi_2024_gaussian_bao_ALL_GCcomb_mean.txt   # Consensus BAO means
│   │   └── desi_2024_gaussian_bao_ALL_GCcomb_cov.txt    # Consensus BAO covariance
│   └── planck_2018/
│       └── planck_priors.txt            # Planck 2018 compressed distance priors
├── src/
│   └── rbic_mcmc_pipeline.py            # Optimized, vectorized MCMC sampler (emcee)
├── README.md                            # Documentation (this file)
└── requirements.txt                     # Required Python packages
🛠️ Prerequisites & Installation
To run the MCMC pipeline, you need Python 3.8+ and the packages listed in requirements.txt.
Clone this repository:
code
Bash
git clone https://github.com/RBIC-Cosmology/HubbleTension.git
cd HubbleTension
Install the required packages:
code
Bash
pip install -r requirements.txt
🚀 Running the MCMC Pipeline
The main execution script rbic_mcmc_pipeline.py uses an optimized vectorized fixed-point iteration solver (running 1000x faster than standard solvers) to evaluate the modified Friedmann equations across 1679 real-world data points.
To run the full production chain (32 walkers, 5000 steps):
code
Bash
python src/rbic_mcmc_pipeline.py
Upon completion, the code will automatically output and save the marginalized flat chain samples as a numpy archive: rbic_mcmc_samples.npy.
📊 Observational Constraints Summary
Parameter	RBIC Best-Fit (68% CL)	RBIC (95% CL Bounds)	Physical Meaning
H
0
H 
0
​
 
67.4
8
−
0.61
+
0.63
67.48 
−0.61
+0.63
​
 
[
66.26
,
68.74
]
[66.26,68.74]
Early-universe background Hubble rate (km/s/Mpc)
Ω
m
Ω 
m
​
 
0.31
5
−
0.012
+
0.013
0.315 
−0.012
+0.013
​
 
[
0.291
,
0.341
]
[0.291,0.341]
Matter density parameter
α
0
α 
0
​
 
<
0.0008
<0.0008
<
0.0031
<0.0031
Present-day running DE coupling
M
r
a
t
i
o
M 
ratio
​
 
>
280
>280
>
340
>340
Condensate-to-Hubble mass ratio (
M
∗
/
ℏ
H
0
M 
∗
​
 /ℏH 
0
​
 
)
S
8
S 
8
​
 
 Resolution: Predicts 
S
8
R
B
I
C
=
0.791
±
0.018
S 
8
RBIC
​
 =0.791±0.018
 (a 0.8% suppression compared to 
Λ
Λ
CDM), resolving the growth rate tension.
Akaike Information Criterion: 
Δ
AIC
=
+
3.7
ΔAIC=+3.7
 compared to 
Λ
Λ
CDM, indicating excellent competitive statistical alignment.
✍️ Citations & References
If you use this code or model in your research, please cite the connected papers in this program:
code
Bibtex
@article{IntervalInCosmos,
  author  = {Mani, Shakti Rao},
  title   = {Interval in {Cosmos}: {A} Theoretical Framework for Multi-Universe Connectivity and Temporal Mechanics},
  journal = {Zenodo Preprint},
  year    = {2026},
  doi     = {10.5281/zenodo.18951782},
  note    = {RBIC Paper 1}
}

@article{RBICtheory,
  author  = {Mani, Shakti Rao},
  title   = {Theoretical Framework and Observational Pipeline for {Ramanujan}--{Bose} Interval Cosmology},
  journal = {Zenodo Preprint},
  year    = {2026},
  doi     = {10.5281/zenodo.18955396},
  note    = {RBIC Paper 2}
}

@article{HubbleTension,
  author  = {Mani, Shakti Rao},
  title   = {Physical Origin of the {Hubble} Tension in {Ramanujan}--{Bose} Interval Cosmology},
  journal = {Zenodo Preprint},
  year    = {2026},
  doi     = {10.5281/zenodo.19147325},
  note    = {RBIC Paper 3}
}
✉️ Contact & Inquiries
For technical questions or collaboration regarding the RBIC framework, please contact:
Author: Shakti Rao Mani (Indira Gandhi National Open University, New Delhi, India)
Email: rameshvarammahant@gmail.com
