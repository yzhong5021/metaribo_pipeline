#!/bin/bash
#
# Annotate metagenome contigs and update the workflow manifest.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=src/scripts/common.sh
source "${SCRIPT_ROOT}/src/scripts/common.sh"

function load_master_values() {
  local master_output

  master_output="$(
    python3 "${SCRIPT_ROOT}/src/helpers/master_config.py" annotation \
      --master-config "${MASTER_CONFIG}"
  )"
  mapfile -t MASTER_VALUES <<< "${master_output}"

  ROOT="${MASTER_VALUES[0]}"
  MASTER_MANIFEST="${MASTER_VALUES[1]}"
  MASTER_PREPROCESS_CONFIG="${MASTER_VALUES[2]}"
  MASTER_METAGENOME_PREPROCESS_CONFIG="${MASTER_VALUES[3]}"
  MASTER_CONTIG_ANNOTATION_CONFIG="${MASTER_VALUES[4]}"
  MASTER_SAMPLESHEET="${MASTER_VALUES[5]}"
  MASTER_LOG="${MASTER_VALUES[6]}"
}

function run_annotation_for_sample() {
  local sample_id="$1"
  local contigs="$2"
  local bam="$4"
  local contig_to_genus
  local cat_prefix
  local cat_contig_calls
  local cat_translated
  local bakta_run_dir
  local coverage_tsv
  local genome_file
  local cds_bed_sorted
  local kofam_tsv
  local kofam_tmp
  local master_sheet_tsv
  local cat_updated
  local cat_total

  contig_to_genus="${CONTIG_TO_GENUS_DIR}/contig_to_genus_${sample_id}.tsv"
  cat_prefix="${CAT_DIR}/${sample_id}"
  cat_contig_calls="${CAT_DIR}/${sample_id}.contig2classification.txt"
  cat_translated="${CAT_DIR}/${sample_id}_translated.txt"
  bakta_run_dir="${BAKTA_DIR}/${sample_id}"
  coverage_tsv="${BAKTA_COVERAGE_DIR}/${sample_id}_gene_mean_coverage.tsv"
  genome_file="${BAKTA_DIR}/${sample_id}.genome"
  cds_bed_sorted="${BAKTA_DIR}/${sample_id}_cds.sorted.bed"
  kofam_tsv="${KOFAM_DIR}/${sample_id}_kofam.tsv"
  kofam_tmp="${ANNOT_WORKDIR}/kofam/${sample_id}"
  master_sheet_tsv="${MASTER_SHEET_DIR}/${sample_id}_master_sheet.tsv"
  mkdir -p "${kofam_tmp}"

  BAKTA_CMD=(
    "${BAKTA}" --db "${BAKTA_DB}"
    --output "${bakta_run_dir}"
    --prefix "${sample_id}"
    --meta
    --keep-contig-headers
    --threads "${THREADS}"
    --force
    "${BAKTA_PARAM_ARGS[@]}"
    "${contigs}"
  )
  "${BAKTA_CMD[@]}"

  eval "$(
    python3 "${HELPER}" process-bakta \
      --sample "${sample_id}" \
      --bakta-run-dir "${bakta_run_dir}" \
      --bakta-outdir "${BAKTA_DIR}" \
      --local-outdir "${LOCAL_OUTDIR}" \
      --shell
  )"
  echo "${sample_id} Number of CDSs annotated by Bakta: ${NUM_CDS}"
  echo "${sample_id} Number of ncRNAs annotated by Bakta: ${NUM_NCRNA}"

  python3 "${HELPER}" sort-bed-for-coverage \
    --fasta "${contigs}" \
    --bed "${CDS_BED}" \
    --sorted-bed "${cds_bed_sorted}" \
    --genome "${genome_file}"
  "${BEDTOOLS}" coverage \
    -a "${cds_bed_sorted}" \
    -b "${bam}" \
    -mean \
    -sorted \
    -g "${genome_file}" > "${coverage_tsv}"
  "${EXEC_ANNOTATION}" \
    -o "${kofam_tsv}" \
    -p "${KOFAM_PROFILES}" \
    -k "${KOFAM_LIST}" \
    --cpu "${THREADS}" \
    --tmp-dir "${kofam_tmp}" \
    -E 1e-3 \
    -f detail-tsv \
    "${CDS_FAA}"

  MASTER_SHEET_CMD=(
    python3 "${HELPER}" make-master-sheet
    --cds-gff "${CDS_GFF}"
    --coverage "${coverage_tsv}"
    --kofam "${kofam_tsv}"
    --output "${master_sheet_tsv}"
  )
  CAT_CMD=(
    "${CAT_PATH}" contigs
    "${CAT_PARAM_ARGS[@]}"
    -c "${contigs}"
    -d "${CAT_DB}"
    -t "${CAT_TAX}"
    -o "${cat_prefix}"
    -n "${THREADS}"
    --force
  )
  ADD_NAMES_CMD=(
    "${CAT_PATH}" add_names
    -i "${cat_contig_calls}"
    -o "${cat_translated}"
    -t "${CAT_TAX}"
    --force
  )
  "${CAT_CMD[@]}"
  "${ADD_NAMES_CMD[@]}"

  cat_updated="$(
    python3 "${HELPER}" merge-cat \
      --contig-to-genus "${contig_to_genus}" \
      --cat-translated "${cat_translated}"
  )"
  cat_total=$(($(wc -l < "${contig_to_genus}") - 1))
  echo "${sample_id} Number of contigs assigned a genus: ${cat_total}"
  echo "${sample_id} Number of contigs assigned with CAT (excluding ones propagated from bins): ${cat_updated}"

  MASTER_SHEET_CMD+=(--contig-to-genus "${contig_to_genus}")
  "${MASTER_SHEET_CMD[@]}"
  python3 "${HELPER}" update-manifest \
    --manifest "${MANIFEST}" \
    --sample "${sample_id}" \
    --cds-gff "${CDS_GFF}" \
    --ncrna-bed "${NCRNA_BED}" \
    --kofam "${kofam_tsv}" \
    --master-sheet "${master_sheet_tsv}"
}

function main() {
  local script_conda_env
  local config_shell

  script_conda_env="${METARIBO_ANNOT_CONDA_ENV:-metaribo_pipeline}"
  metaribo_reexec_in_conda "${script_conda_env}" "$@"
  metaribo_parse_master_config "${SCRIPT_ROOT}" "$@"
  load_master_values

  HELPER="${ROOT}/src/helpers/annotation_helper.py"
  ANNOTATION_CONFIG="${ANNOTATION_CONFIG:-${MASTER_CONTIG_ANNOTATION_CONFIG}}"
  MANIFEST="${METAGENOME_MANIFEST:-${MASTER_MANIFEST}}"
  SAMPLESHEET="${METAGENOME_SAMPLESHEET:-${MASTER_SAMPLESHEET}}"
  LOG="${METAGENOME_LOG:-${MASTER_LOG}}"

  metaribo_require_no_args "${PASSTHRU_ARGS[@]}"
  metaribo_start_log "${LOG}"

  cd "${ROOT}"

  echo "[$(date -Is)] - Starting metagenome annotation"
  config_shell="$(
    python3 "${HELPER}" config-shell --annotation-config "${ANNOTATION_CONFIG}"
  )"
  eval "${config_shell}"
  python3 "${HELPER}" check-paths \
    --manifest "${MANIFEST}" \
    --samplesheet "${SAMPLESHEET}" \
    --annotation-config "${ANNOTATION_CONFIG}"

  mkdir -p \
    "${CAT_DIR}" \
    "${BAKTA_DIR}" \
    "${KOFAM_DIR}" \
    "${ANNOT_WORKDIR}" \
    "${CONTIG_TO_GENUS_DIR}" \
    "${BAKTA_CDS_GFF_DIR}" \
    "${BAKTA_NCRNA_BED_DIR}" \
    "${BAKTA_CDS_BED_DIR}" \
    "${BAKTA_CDS_FAA_DIR}" \
    "${BAKTA_COVERAGE_DIR}" \
    "${MASTER_SHEET_DIR}"
  read -r -a CAT_PARAM_ARGS <<< "${CAT_PARAMS}"
  read -r -a BAKTA_PARAM_ARGS <<< "${BAKTA_PARAMS}"

  while IFS=$'\t' read -r SAMPLE_ID CONTIGS CONTIG_TO_GENUS BAM; do
    [[ "${SAMPLE_ID}" == "sample_id" ]] && continue
    run_annotation_for_sample "${SAMPLE_ID}" "${CONTIGS}" "${CONTIG_TO_GENUS}" "${BAM}"
  done < <(
    python3 "${HELPER}" list-inputs \
      --manifest "${MANIFEST}" \
      --samplesheet "${SAMPLESHEET}" \
      --annotation-config "${ANNOTATION_CONFIG}"
  )

  echo "[$(date -Is)] - Metagenome annotation complete"
}

main "$@"
