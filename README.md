# emicmGapTime

# Weighted EM–ICM for Bivariate Interval Censoring

Open‑source reference implementation for the paper  
“Non‑parametric estimation of sequential gap times under
panel censoring” (Azevedo & All, 2025).

## Repository structure

| Folder      | Purpose |
|-------------|---------|
| src/      | Core R functions: Turnbull geometry and weighted EM–ICM. |
| sim/      | Monte‑Carlo driver (run_sim_grid.R) that reproduces all results in Section 4. |
| tandmobiel/ | Notebook for the dental‑caries cohort (Section 4.3). |
| macs/     | Notebook for the MACS HIV/AIDS example (Section 4.4). |

## Requirements

* *R ≥ 4.1*  
* CRAN packages (installed automatically on first run):  
  data.table, Matrix, copula, isotone, ranger,  
  SQUAREM, future, future.apply, ggplot2,  
  plus example‑specific packages  
  * icensBKL for Tandmobiel data  
  * KMsurv   for MACS data

## Quick start

```bash
git clone https://github.com/<user>/emicmGapTime.git
cd emicmGapTime
Rscript sim/run_sim_grid.R            # Monte‑Carlo study
Rscript tandmobiel/analysis_tandmobiel.R
Rscript macs/analysis_macs.R
