# Functional Regression via Set Embedding — Simulation Code

Simulation code accompanying:

> Wong, R. K. W., Cao, G., and Li, Y. "Functional Regression via Set Embedding." 

This repository contains the R code used to run the simulation study in Section 6.1 of the
paper and Section S.4 of the Supplementary Material (Tables 1–2 in the main text and
Tables S.1–S.7 in the supplement).

## What the simulation does

For each configuration in a grid of

- sample size `n` (50, 100),
- within-curve observation count `M` (fixed at 10, 30, or 60; or random on {10,...,20}),
- response link function (`flr_linear`, `flr_square`, `flr_cube` — linear, quadratic, and
  cubic-type nonlinear),
- score distribution for the functional predictor's Karhunen–Loève coefficients
  (`gaussian`, `laplace`, `t5`, `gamma`),
- measurement-error ("X-noise") distribution (`gaussian`, `t3`, `gamma`), and
- estimation method (see mapping below),

the script generates 500 (`N_REPS`) independent Monte Carlo replications, fits the specified
estimator, and records the prediction mean-squared error on an independent test set. Results
are written per configuration and then aggregated into a single `all_method_config_summary.csv`
per run (mean/SD of prediction MSE and other diagnostics), which is what the paper's tables
report.

### Method name mapping

The paper compares the proposed estimator against five benchmarks (Section 6.1). The
`METHOD_GRID` values used in the code map to the paper as follows:

| Code (`method`) | Paper                                                              |
|------------------|--------------------------------------------------------------------|
| `embedding`      | Proposed set-embedding KRR estimator ("KRR" in the tables)         |
| `plugin`         | Plug-in estimator of Hall & Horowitz (2007)                        |
| `in`             | Integral-approximation scores estimator ("IN"), Yao et al. (2005a) |
| `pace`           | Conditional-expectation scores estimator ("PACE"), Yao et al. (2005a) |
| `pace_rkhs`      | RKHS-based nonparametric estimator ("RKHS"), Avery et al. (2014)   |
| `split_s5`       | Sample-splitting estimator ("SS"), Zhou et al. (2023)              |

## Repository structure

```
.
├── code/
│   ├── fdr.R      # main simulation driver (data generation, all six estimators, CV/GCV tuning)
│   └── fdr.sh      # SLURM array-job submission script used to run fdr.R on an HPC cluster
├── results_summary/
│   └── *_summary.csv   # per-run aggregated results (mean/SD prediction MSE per configuration);
│                         # directly reproduces Tables 1–2 (main text) and S.1–S.7 (supplement)
├── LICENSE
└── README.md
```

Column definitions:

| Column | Meaning |
|---|---|
| `run_tag` | Encodes `n`, `obs_count`, `obs_count_mode`, `y_mode`, `score_dist`, `x_noise_dist`, `method` |
| `n_reps` / `n_success` | Number of Monte Carlo replications requested / completed without error |
| `mean_pred_mse`, `sd_pred_mse` | Mean and standard deviation of test-set prediction MSE across replications — these are the numbers reported in the paper's tables |
| `mean_runtime_seconds`, `total_runtime_seconds` | Wall-clock timing, used for the computational-cost discussion in Section 6.1 |

**Note on raw output.** The full per-replication output (every Monte Carlo draw, ~3,300 CSV
files / ~286 MB from the original HPC runs) is not included here to keep the repository light.
It is fully reproducible from `fdr.R` with the seeds recorded in each `results_summary/*.csv`,
and is available from the authors on request.

## Requirements

- R (≥ 4.3 recommended; developed against `R-bundle-CRAN/2023.12-foss-2023a`)
- R packages: `fdapace`, `parallel`, `VGAM`, `splines`

Install the packages with:

```r
install.packages(c("fdapace", "VGAM"))
# 'parallel' and 'splines' ship with base R
```

## Running the simulation

The driver script supports two modes, controlled by the first command-line argument.

**1. Build the configuration manifest** (lists every `n` × `obs_count` × ... × `method`
combination implied by the grid environment variables, without running anything):

```bash
Rscript code/fdr.R manifest
```

**2. Run a single configuration** by its row index in that manifest (this is what each SLURM
array task does):

```bash
Rscript code/fdr.R task <task_id>
```

The grid and run settings are controlled entirely through environment variables (all have
defaults baked into `fdr.R`):

| Variable | Meaning | Example |
|---|---|---|
| `N_GRID` | Sample size(s) | `50,100` |
| `OBS_COUNT_GRID` | Within-curve observation count(s) | `10,30,60` |
| `OBS_COUNT_MODE` | `fixed` (exactly `obs_count` points) or `random` (discrete uniform up to `2*obs_count`) | `fixed` |
| `Y_MODE_GRID` | Response link function(s) | `flr_linear,flr_square,flr_cube` |
| `SCORE_DIST_GRID` | KL-coefficient distribution(s) | `gaussian,laplace,t5,gamma` |
| `X_NOISE_DIST_GRID` | Measurement-error distribution(s) | `gaussian,t3,gamma` |
| `METHOD_GRID` | Estimator(s) to run, see mapping above | `embedding` |
| `SEED_START` | First random seed | `42` |
| `N_REPS` | Monte Carlo replications per configuration | `500` |
| `N_CORES` | Cores for within-task parallelism | `1` |
| `RESULT_DIR` | Output directory | `results_embedding_50` |

Example — reproduce the proposed method's results at `n = 50` (matches
`results_summary/embedding_50_summary.csv`):

```bash
export N_GRID=50
export OBS_COUNT_GRID=10
export METHOD_GRID=embedding
export RESULT_DIR=results_embedding_50
Rscript code/fdr.R manifest
Rscript code/fdr.R task 1   # repeat for each row of the manifest, or loop over all rows
```

### Running on a SLURM cluster

`code/fdr.sh` is the batch script used on the authors' cluster. It sets the same environment
variables, loads the R module, and runs one manifest row per SLURM array task:

```bash
sbatch --array=1-N code/fdr.sh
```

You will need to adjust the `--account`, `module load` line, and `R_LIBS_USER` path for your
own cluster; `N` should match the number of rows printed by `Rscript code/fdr.R manifest` for
your chosen grid.

## Citation

If you use this code, please cite:

```bibtex
@article{wong_cao_li_2026_setembedding,
  title   = {Functional Regression via Set Embedding},
  author  = {Wong, Raymond K. W. and Cao, Guanqun and Li, Yehua},
  journal = {},
  year    = {2026},
  note    = {Submitted}
}
```
*(Update the BibTeX entry with the final volume/issue/page numbers once the paper is published.)*

## Contact

Questions about the code can be directed to Guanqun Cao (caoguanq@msu.edu).
