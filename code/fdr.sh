#!/usr/bin/env bash
#SBATCH --job-name=fdr_2_pace_rkhs
#SBATCH --account=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=1-12
#SBATCH --output=logs/fdr_2-%A_%a.out
#SBATCH --error=logs/fdr_2-%A_%a.err

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR}"
mkdir -p "${WORKDIR}/logs"

export RESULT_DIR="${RESULT_DIR:-${WORKDIR}/results_plugin_n100}"
export N_CORES="${N_CORES:-1}"

export SEED_START="${SEED_START:-42}"
export N_REPS="${N_REPS:-50}"
export N_GRID="${N_GRID:-100}"
export OBS_COUNT_GRID="${OBS_COUNT_GRID:-10}"
#10,30,60
export OBS_COUNT_MODE="${OBS_COUNT_MODE:-random}"
#fixed,random
#random=2*fixed
export Y_MODE_GRID="${Y_MODE_GRID:-flr_linear,flr_square,flr_cube}"
#flr_linear,flr_square,flr_cube
export SCORE_DIST_GRID="${SCORE_DIST_GRID:-gaussian,laplace,t5,gamma}"
#gaussian,laplace,t5,gamma
export X_NOISE_DIST_GRID="${X_NOISE_DIST_GRID:-gaussian}"
#gaussian,t3,gamma
export METHOD_GRID="${METHOD_GRID:-plugin}"
#embedding,split_s5,plugin,pace,in,pace_rkhs

if ! command -v module >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh
fi
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL _CE_CONDA _CE_M || true
export PATH="$(echo "${PATH}" | tr ':' '\n' | grep -v -i 'conda\|anaconda' | paste -sd: -)"
export LD_LIBRARY_PATH="$(echo "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v -i 'conda\|anaconda' | paste -sd: -)"

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
mkdir -p "$R_LIBS_USER" "$RESULT_DIR"

cd "$WORKDIR"

echo "[fdr.sh] WORKDIR=$WORKDIR"
echo "[fdr.sh] RESULT_DIR=$RESULT_DIR"
echo "[fdr.sh] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "[fdr.sh] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-NA}"
echo "[fdr.sh] N_CORES=${N_CORES}"
echo "[fdr.sh] SEED_START=${SEED_START}"
echo "[fdr.sh] N_REPS=${N_REPS}"
echo "[fdr.sh] N_GRID=${N_GRID}"
echo "[fdr.sh] OBS_COUNT_GRID=${OBS_COUNT_GRID}"
echo "[fdr.sh] OBS_COUNT_MODE=${OBS_COUNT_MODE}"
echo "[fdr.sh] Y_MODE_GRID=${Y_MODE_GRID}"
echo "[fdr.sh] SCORE_DIST_GRID=${SCORE_DIST_GRID}"
echo "[fdr.sh] X_NOISE_DIST_GRID=${X_NOISE_DIST_GRID}"
echo "[fdr.sh] METHOD_GRID=${METHOD_GRID}"

Rscript fdr.R manifest >/dev/null 2>&1 || true

srun --export=ALL \
  --hint=nomultithread --cpu-bind=cores \
  Rscript fdr_2.R task "${SLURM_ARRAY_TASK_ID}"
