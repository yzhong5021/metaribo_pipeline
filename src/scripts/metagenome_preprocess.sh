#!/bin/bash

# Preprocess metagenomic data.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly SCRIPT_ROOT

source "${SCRIPT_ROOT}/src/scripts/common.sh"

function main() {
  local -a master_values
  local -a metagenome_preprocess_cmd
  local config_shell
  local master_output

  metaribo_parse_master_config "${SCRIPT_ROOT}" "$@"
  master_output="$(
    python3 "${SCRIPT_ROOT}/src/helpers/master_config.py" metagenome \
      --master-config "${MASTER_CONFIG}"
  )"
  mapfile -t master_values <<< "${master_output}"

  ROOT="${master_values[0]}"
  MASTER_MANIFEST="${master_values[1]}"
  MASTER_RIBO_PREPROCESS_CONFIG="${master_values[2]}"
  MASTER_METAGENOME_PREPROCESS_CONFIG="${master_values[3]}"
  MASTER_CONTIG_ANNOTATION_CONFIG="${master_values[4]}"
  MASTER_LOG="${master_values[5]}"

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

  metagenome_preprocess_cmd=(
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
    metagenome_preprocess_cmd+=(--gtdb_db "${GTDB_DB}")
  fi

  echo "[$(date -Is)]: Starting metagenome preprocessing"

  python3 "${HELPER}" check-paths \
    --metagenome-preprocess-config "${METAGENOME_PREPROCESS_CONFIG}"
  "${metagenome_preprocess_cmd[@]}"
  python3 "${HELPER}" collect-outputs \
    --manifest "${MANIFEST}" \
    --samplesheet "${SAMPLESHEET}" \
    --metagenome-preprocess-config "${METAGENOME_PREPROCESS_CONFIG}"

  echo "[$(date -Is)]: Metagenome preprocessing complete"
}

main "$@"
