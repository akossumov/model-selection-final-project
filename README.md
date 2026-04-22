# Parametric and Semiparametric Copula Model Selection

**Author:** Aibat Kossumov  
**Course:** MATH 782 - Model Selection and Sparsity in Statistical Learning  

## About This Repository
This repository contains the code, report, and presentation slides for my final project on copula model selection. The project explores how to choose the best copula model to capture the dependence structure between variables. 

It compares model selection criteria across two distinct statistical frameworks:
1. **Fully Parametric Approach:** Using Two-Stage Maximum Likelihood Estimation (2MLE).
2. **Semiparametric Approach:** Using Maximum Pseudo-Likelihood Estimation (MPLE), where the marginal distributions are estimated nonparametrically.

## How to Run the Code
The R code for the simulations is organized into two main folders corresponding to the two approaches.

* **Two-Stage MLE Simulations:** To run the simulations for the fully parametric approach, navigate to the `two_stage_MLE` folder and execute the `simulation_two_stage_MLE.R` script.

* **Pseudo-Likelihood (MPLE) Simulations:** To run the simulations evaluating the semiparametric criteria, navigate to the `MPLE` folder and execute the `main.R` script.

## Files Included
* `Report.pdf`: The complete final project report detailing the methodology and simulation results.
* `Presentation.pdf`: The Beamer slide deck used for the project defense.
* `two_stage_MLE/`: Directory containing the parametric simulation scripts.
* `MPLE/`: Directory containing the semiparametric simulation scripts.
