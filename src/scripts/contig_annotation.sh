#!/bin/bash

# Annotate metagenome contigs.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly SCRIPT_ROOT
source "${SCRIPT_ROOT}/src/scripts/common.sh"

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi

function run_annotation_for_sample() {
  local sample_id="$1"
  local contigs="$2"
  local bam="$3"
  local manifest_contig_to_genus="$4"
  local -a add_names_cmd
  local -a bakta_cmd
  local -a cat_cmd
  local -a master_sheet_cmd
  local -a merge_cat_cmd
  local -a update_manifest_cmd
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
  local cat_line_count
  local bakta_shell
  local has_genus_map

  has_genus_map=1
  contig_to_genus="${manifest_contig_to_genus}"
  if [[ "${manifest_contig_to_genus}" == "NA" ]]; then
    has_genus_map=0
    contig_to_genus="${CONTIG_TO_GENUS_DIR}/contig_to_genus_${sample_id}.tsv"
  fi
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

  bakta_cmd=(
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
  "${bakta_cmd[@]}"

  bakta_shell="$(
    python3 "${HELPER}" process-bakta \
      --sample "${sample_id}" \
      --bakta-run-dir "${bakta_run_dir}" \
      --bakta-outdir "${BAKTA_DIR}" \
      --local-outdir "${LOCAL_OUTDIR}" \
      --shell
  )"
  eval "${bakta_shell}"
  echo "${sample_id} Number of CDSs annotated by Bakta: ${NUM_CDS}"
  echo "${sample_id} Number of ncRNA regions annotated by Bakta: ${NUM_NCRNA}"

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

  master_sheet_cmd=(
    python3 "${HELPER}" make-master-sheet
    --cds-gff "${CDS_GFF}"
    --coverage "${coverage_tsv}"
    --kofam "${kofam_tsv}"
    --output "${master_sheet_tsv}"
  )
  if [[ "${RUN_CAT}" -eq 1 ]]; then
    cat_cmd=(
      "${CAT_PATH}" contigs
      "${CAT_PARAM_ARGS[@]}"
      -c "${contigs}"
      -d "${CAT_DB}"
      -t "${CAT_TAX}"
      -o "${cat_prefix}"
      -n "${THREADS}"
      --force
    )
    add_names_cmd=(
      "${CAT_PATH}" add_names
      -i "${cat_contig_calls}"
      -o "${cat_translated}"
      -t "${CAT_TAX}"
      --force
    )
    "${cat_cmd[@]}"
    "${add_names_cmd[@]}"

    merge_cat_cmd=(
      python3 "${HELPER}" merge-cat
      --contig-to-genus "${contig_to_genus}"
      --cat-translated "${cat_translated}"
    )
    if [[ "${has_genus_map}" -eq 0 ]]; then
      merge_cat_cmd+=(--reset-existing)
    fi
    cat_updated="$("${merge_cat_cmd[@]}")"
    has_genus_map=1
    cat_line_count="$(wc -l < "${contig_to_genus}")"
    cat_total=$((cat_line_count - 1))
    echo "${sample_id} Number of contigs assigned a genus: ${cat_total}"
    echo "${sample_id} Number of contigs assigned with CAT (excluding ones propagated from bins): ${cat_updated}"
  else
    echo "${sample_id} skipping CAT classification because cat_db/cat_tax is missing"
  fi

  if [[ "${has_genus_map}" -eq 1 && -s "${contig_to_genus}" ]]; then
    master_sheet_cmd+=(--contig-to-genus "${contig_to_genus}")
  fi
  "${master_sheet_cmd[@]}"
  update_manifest_cmd=(
    python3 "${HELPER}" update-manifest
    --manifest "${MANIFEST}"
    --sample "${sample_id}"
    --cds-gff "${CDS_GFF}"
    --ncrna-bed "${NCRNA_BED}"
    --kofam "${kofam_tsv}"
    --master-sheet "${master_sheet_tsv}"
  )
  if [[ "${has_genus_map}" -eq 1 ]]; then
    update_manifest_cmd+=(--contig-to-genus "${contig_to_genus}")
  fi
  "${update_manifest_cmd[@]}"
}

function main() {
  local -a master_values
  local config_shell
  local contigs
  local annotation_inputs
  local master_output
  local manifest_contig_to_genus
  local sample_id
  local bam

  metaribo_parse_master_config "${SCRIPT_ROOT}" "$@"
  master_output="$(
    python3 "${SCRIPT_ROOT}/src/helpers/master_config.py" annotation \
      --master-config "${MASTER_CONFIG}"
  )"
  mapfile -t master_values <<< "${master_output}"

  ROOT="${master_values[0]}"
  MASTER_MANIFEST="${master_values[1]}"
  MASTER_CONTIG_ANNOTATION_CONFIG="${master_values[2]}"
  MASTER_SAMPLESHEET="${master_values[3]}"
  MASTER_LOG="${master_values[4]}"

  HELPER="${ROOT}/src/helpers/annotation_helper.py"
  ANNOTATION_CONFIG="${ANNOTATION_CONFIG:-${MASTER_CONTIG_ANNOTATION_CONFIG}}"
  MANIFEST="${METAGENOME_MANIFEST:-${MASTER_MANIFEST}}"
  SAMPLESHEET="${METAGENOME_SAMPLESHEET:-${MASTER_SAMPLESHEET}}"
  LOG="${METAGENOME_LOG:-${MASTER_LOG}}"

  metaribo_start_log "${LOG}"

  cd "${ROOT}"

  echo "[$(date -Is)]: Starting metagenome annotation"
  python3 "${HELPER}" check-paths \
    --manifest "${MANIFEST}" \
    --samplesheet "${SAMPLESHEET}" \
    --annotation-config "${ANNOTATION_CONFIG}"
  config_shell="$(
    python3 "${HELPER}" config-shell --annotation-config "${ANNOTATION_CONFIG}"
  )"
  eval "${config_shell}"

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

  annotation_inputs="$(
    python3 "${HELPER}" list-inputs \
      --manifest "${MANIFEST}" \
      --samplesheet "${SAMPLESHEET}" \
      --annotation-config "${ANNOTATION_CONFIG}"
  )"
  while IFS=$'\t' read -r sample_id contigs manifest_contig_to_genus bam; do
    [[ "${sample_id}" == "sample_id" ]] && continue
    run_annotation_for_sample \
      "${sample_id}" \
      "${contigs}" \
      "${bam}" \
      "${manifest_contig_to_genus}"
  done <<< "${annotation_inputs}"

  echo "[$(date -Is)]: Metagenome annotation complete"
}

main "$@"
