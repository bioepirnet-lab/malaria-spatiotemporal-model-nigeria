## Description: malaria-spatiotemporal-model-nigeria
Code for the study “Spatio-Temporal and Mathematical Modeling of Malaria Burden and Intervention Efficacy in Nigeria”, including R (INLA Bayesian spatio-temporal model) and Python (ODE transmission model) implementations.


# Malaria Spatio-Temporal Model – Nigeria

This repository contains the **code used in the study**:

**"Spatio-Temporal and Mathematical Modeling of Malaria Burden and Intervention Efficacy in Nigeria."**

The repository provides reproducible implementations of the **statistical and mathematical models** used to analyze malaria prevalence patterns and intervention impacts in Nigeria.

---

# Overview

The study combines two complementary modeling approaches:

## 1. Bayesian Spatio-Temporal Statistical Model (R – INLA)

Used to analyze malaria prevalence among children under five using **Nigeria Malaria Indicator Survey (NMIS)** data from:

- 2010  
- 2015  
- 2021  

## 2. Mathematical Transmission Model (Python – ODE System)

A compartmental model describing malaria transmission dynamics between humans and mosquitoes.  
The model explores the impact of interventions such as:

- Seasonal Malaria Chemoprevention (SMC)
- Vaccination
- Bed net use
- Environmental factors (e.g., irrigation and mosquito breeding sites)

---

# Repository Structure
malaria-spatiotemporal-model-nigeria/
/code

├── code
│ ├── spatiotemporal_malaria_inla_model.R
│ └── malaria_transmission_ode_model.py
│
├── figures
│
├── outputs
│
├── data_info
│ └── nmis_data_description.txt
│
└── README.md
/code


---

# Statistical Model

The statistical analysis uses a **Bayesian hierarchical spatio-temporal model** where malaria infection status is modeled at the individual child level using a **Bernoulli logistic regression framework**.

### Model components

- Fixed effects for socio-demographic and household variables
- Spatial random effects using **Intrinsic Conditional Autoregressive (ICAR) priors**
- Temporal random effects modeled as a **first-order random walk**
- Spatio-temporal interaction terms
- Penalized Complexity (PC) priors for random-effect precision parameters

Model estimation was performed using **Integrated Nested Laplace Approximation (INLA)**.

---

# Mathematical Model

The mathematical model is a **compartmental system of ordinary differential equations (ODEs)** describing malaria transmission between humans and mosquito vectors.

### Human compartments

- Susceptible
- Vaccinated
- Exposed
- Asymptomatic
- Infectious
- Recovered

### Vector compartments

- Susceptible
- Exposed
- Infectious

### Model features

The model incorporates environmental and intervention factors such as:

- Mosquito breeding sites
- Irrigation effects
- Bed net usage
- Treatment access
- Seasonal malaria chemoprevention
- Vaccination

The model simulations are implemented using **Python and SciPy**.

---

# Data Availability

The analysis uses **Nigeria Malaria Indicator Survey (NMIS)** data.

Due to data access restrictions, the NMIS datasets are **not included in this repository**.

Researchers wishing to reproduce the analysis should obtain the data from the official sources and use the scripts provided here.

---

# Software Requirements

## R

Required packages include:
- INLA
- spdep
- sf
- dplyr
- ggplot2


## Python

Required libraries include:
- numpy
- scipy
- matplotlib
- pandas
---

# Reproducibility

All scripts required to reproduce the statistical analysis and model simulations are provided in this repository.

The code is intended to promote **transparency, reproducibility, and further research on malaria transmission dynamics**.

---

# Citation

If you use this code, please cite:

Ukwajunor, E. E., Nwankwo-Eluwa, A., Egonmwan, A. O., & Anyasodor, A. E.
**Spatio-Temporal and Mathematical Modeling of Malaria Burden and Intervention Efficacy in Nigeria.**
*Communications in Statistics – Case Studies and Data Analysis.*

---

# Contact

**Eunice Ukwajunor, PhD**
Founder & Lead Research Scientist
Biostatistics and Epidemiological Research Network (BioEpiRNet)
bioepirnet@gmail.com | +234 8023437236
