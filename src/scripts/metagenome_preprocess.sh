#!/bin/bash
#
# Run metagenome preprocessing and collect outputs for downstream stages.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=src/scripts/common.sh
source "${SCRIPT_ROOT}/src/scripts/common.sh"

function load_master_values() {
  local master_output

  master_output="$(
    python3 "${SCRIPT_ROOT}/src/helpers/master_config.py" metagenome \
      --master-config "${MASTER_CONFIG}"
  )"
  mapfile -t MASTER_VALUES <<< "${master_output}"

  ROOT="${MASTER_VALUES[0]}"
  MASTER_MANIFEST="${MASTER_VALUES[1]}"
  MASTER_PREPROCESS_CONFIG="${MASTER_VALUES[2]}"
  MASTER_METAGENOME_PREPROCESS_CONFIG="${MASTER_VALUES[3]}"
  MASTER_CONTIG_ANNOTATION_CONFIG="${MASTER_VALUES[4]}"
  MASTER_LOG="${MASTER_VALUES[5]}"
}

function main() {
  local script_conda_env
  local config_shell

  script_conda_env="${METARIBO_METAGENOME_PREPROCESS_CONDA_ENV:-metaribo_pipeline}"
  metaribo_reexec_in_conda "${script_conda_env}" "$@"
  metaribo_parse_master_config "${SCRIPT_ROOT}" "$@"
  load_master_values

  HELPER="${ROOT}/src/helpers/metagenome_preprocess_helper.py"
  METAGENOME_PREPROCESS_CONFIG="${METAGENOME_PREPROCESS_CONFIG:-${MASTER_METAGENOME_PREPROCESS_CONFIG}}"
  MANIFEST="${METAGENOME_MANIFEST:-${MASTER_MANIFEST}}"
  LOG="${METAGENOME_LOG:-${MASTER_LOG}}"

  metaribo_require_no_args "${PASSTHRU_ARGS[@]}"
  metaribo_start_log "${LOG}"

  cd "${ROOT}"
  unset JAVA_CMD
  unset JAVA_HOME

  config_shell="$(
    python3 "${HELPER}" command-shell \
      --metagenome-preprocess-config "${METAGENOME_PREPROCESS_CONFIG}"
  )"
  eval "${config_shell}"

  export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-}"
  export NXF_APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR}"

  METAGENOME_PREPROCESS_CMD=(
    "${NEXTFLOW}" run nf-core/mag
    -r 5.4.2
    -profile apptainer
    -params-file "${PARAMS_FILE}"
    -c "${NEXTFLOW_CONFIG}"
    -work-dir "${WORK_DIR}"
    --input "${SAMPLESHEET}"
    --outdir "${OUTDIR}"
    -resume
  )
  if [[ -n "${GTDB_DB:-}" ]]; then
    METAGENOME_PREPROCESS_CMD+=(--gtdb_db "${GTDB_DB}")
  fi

  echo "[$(date -Is)] - Starting metagenome preprocessing"

  python3 "${HELPER}" check-paths \
    --metagenome-preprocess-config "${METAGENOME_PREPROCESS_CONFIG}"
  "${METAGENOME_PREPROCESS_CMD[@]}"
  python3 "${HELPER}" collect-outputs \
    --manifest "${MANIFEST}" \
    --samplesheet "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --output-dir "${LOCAL_OUTDIR}" \
    --metagenome-preprocess-config "${METAGENOME_PREPROCESS_CONFIG}"

  echo "[$(date -Is)] - Metagenome preprocessing complete"
}

main "$@"
