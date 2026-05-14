#!/bin/bash

# Trim, align and filter Ribo-ITP reads.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly SCRIPT_ROOT

source "${SCRIPT_ROOT}/src/scripts/common.sh"

function count_fastq() {
  python3 "${HELPER}" count-fastq --fastq "$1"
}

function count_bam() {
  "${SAMTOOLS}" view -c "$1"
}

function append_stats() {
  python3 "${HELPER}" append-stats \
    --stats "$1" \
    --sample "$2" \
    --stage "$3" \
    --before "$4" \
    --after "$5"
}

function run_ribo_sample() {
  local ribo_id="$1"
  local meta_id="$2"
  local ribo_fq="$3"
  local contigs="$4"
  local ncrna_bed="$5"
  local -a cutadapt_cmd
  local -a umi_dedup_cmd
  local -a umi_extract_cmd
  local sample_dir
  local index_dir
  local index_prefix
  local stats
  local trimmed_fq
  local cutadapt_log
  local umi_fq
  local reads_for_align
  local sorted_bam
  local unique_bam
  local filtered_bam
  local dedup_bam
  local ref_fasta
  local input_reads
  local trimmed_reads
  local align_input_reads
  local umi_reads
  local sorted_reads
  local alignment_bam
  local alignment_reads
  local unique_reads
  local filtered_reads
  local dedup_reads

  sample_dir="${RIBO_DIR}/${ribo_id}"
  index_dir="${BT2_DIR}/${ribo_id}"
  index_prefix="${index_dir}/${ribo_id}"
  stats="${sample_dir}/${ribo_id}.alignment_stats.txt"
  trimmed_fq="${sample_dir}/${ribo_id}.trimmed.fq.gz"
  cutadapt_log="${sample_dir}/${ribo_id}.cutadapt.log"
  umi_fq="${sample_dir}/${ribo_id}.trimmed.umi.fq.gz"
  reads_for_align="${trimmed_fq}"
  sorted_bam="${sample_dir}/${ribo_id}.sorted.bam"
  unique_bam="${sample_dir}/${ribo_id}.no_multi.sorted.bam"
  filtered_bam="${sample_dir}/${ribo_id}.nomulti.filtered.sorted.bam"
  dedup_bam="${sample_dir}/${ribo_id}.nomulti.filtered.sorted.deduped.bam"

  if [[ "${CONTAINS_UMI}" -eq 1 ]]; then
    sorted_bam="${sample_dir}/${ribo_id}.sorted.umi.bam"
    unique_bam="${sample_dir}/${ribo_id}.no_multi.sorted.umi.bam"
    filtered_bam="${sample_dir}/${ribo_id}.nomulti.filtered.sorted.umi.bam"
    reads_for_align="${umi_fq}"
  fi

  mkdir -p "${sample_dir}" "${index_dir}"
  echo "${ribo_id} preprocessing against metagenome ${meta_id}"

  cutadapt_cmd=(
    "${CUTADAPT}"
    "${CUTADAPT_ARG_ARRAY[@]}"
    -j "${THREADS}"
    -o "${trimmed_fq}"
    "${ribo_fq}"
  )
  umi_extract_cmd=(
    "${UMI_TOOLS}" extract
    "${UMITOOLS_EXTRACT_ARG_ARRAY[@]}"
    --stdin="${trimmed_fq}"
    --stdout="${umi_fq}"
  )

  ref_fasta="$(
    python3 "${HELPER}" prepare-reference \
      --contigs "${contigs}" \
      --sample "${ribo_id}" \
      --reference-dir "${REFERENCE_DIR}"
  )"
  "${cutadapt_cmd[@]}" 2>&1 | tee "${cutadapt_log}"
  input_reads="$(count_fastq "${ribo_fq}")"
  trimmed_reads="$(count_fastq "${trimmed_fq}")"
  python3 "${HELPER}" init-stats --stats "${stats}" --sample "${ribo_id}"
  append_stats "${stats}" "${ribo_id}" "cutadapt" \
    "${input_reads}" "${trimmed_reads}"

  align_input_reads="${trimmed_reads}"
  if [[ "${CONTAINS_UMI}" -eq 1 ]]; then
    "${umi_extract_cmd[@]}"
    umi_reads="$(count_fastq "${umi_fq}")"
    append_stats "${stats}" "${ribo_id}" "umi_extract" \
      "${trimmed_reads}" "${umi_reads}"
    align_input_reads="${umi_reads}"
  fi

  "${BOWTIE2_BUILD}" --threads "${THREADS}" -f "${ref_fasta}" "${index_prefix}"
  "${BOWTIE2}" \
    --seed "${SEED}" \
    -x "${index_prefix}" \
    "${BT2_ARG_ARRAY[@]}" \
    -U "${reads_for_align}" \
    | "${SAMTOOLS}" sort -@ "${THREADS}" -o "${sorted_bam}" -
  "${SAMTOOLS}" index "${sorted_bam}"
  sorted_reads="$(count_bam "${sorted_bam}")"
  append_stats "${stats}" "${ribo_id}" "alignment" \
    "${align_input_reads}" "${sorted_reads}"

  alignment_bam="${sorted_bam}"
  alignment_reads="${sorted_reads}"
  if [[ "${REMOVE_MULTIMAPPERS}" -eq 1 ]]; then
    "${SAMTOOLS}" view -h -F 256 "${sorted_bam}" \
      | awk 'BEGIN{FS=OFS="\t"} /^@/ || $0 !~ /\tXS:i:/' \
      | "${SAMTOOLS}" view -b -o "${unique_bam}" -
    "${SAMTOOLS}" index "${unique_bam}"
    unique_reads="$(count_bam "${unique_bam}")"
    append_stats "${stats}" "${ribo_id}" "remove_multimappers" \
      "${alignment_reads}" "${unique_reads}"
    alignment_bam="${unique_bam}"
    alignment_reads="${unique_reads}"
  fi

  if [[ "${REMOVE_NCRNA}" -eq 1 ]]; then
    "${BEDTOOLS}" intersect -abam "${alignment_bam}" -b "${ncrna_bed}" -v \
      | "${SAMTOOLS}" sort -@ "${THREADS}" -o "${filtered_bam}" -
    "${SAMTOOLS}" index "${filtered_bam}"
    filtered_reads="$(count_bam "${filtered_bam}")"
    append_stats "${stats}" "${ribo_id}" "remove_ncrna" \
      "${alignment_reads}" "${filtered_reads}"
    alignment_bam="${filtered_bam}"
    alignment_reads="${filtered_reads}"
  fi

  if [[ "${CONTAINS_UMI}" -eq 1 ]]; then
    umi_dedup_cmd=(
      "${UMI_TOOLS}" dedup
      "${UMITOOLS_DEDUP_ARG_ARRAY[@]}"
      --stdin="${alignment_bam}"
      --stdout="${dedup_bam}"
    )
    "${umi_dedup_cmd[@]}"
    "${SAMTOOLS}" index "${dedup_bam}"
    dedup_reads="$(count_bam "${dedup_bam}")"
    append_stats "${stats}" "${ribo_id}" "umi_dedup" \
      "${alignment_reads}" "${dedup_reads}"
    alignment_bam="${dedup_bam}"
  fi

  python3 "${HELPER}" update-manifest \
    --manifest "${MANIFEST}" \
    --sample "${ribo_id}" \
    --alignment-bam "${alignment_bam}" \
    --stats "${stats}" \
    --local-outdir "${LOCAL_OUTDIR}"
  echo "${ribo_id} finished ribo preprocessing."
}

function main() {
  local -a master_values
  local config_shell
  local contigs
  local master_output
  local meta_id
  local ncrna_bed
  local ribo_fq
  local ribo_id

  metaribo_parse_master_config "${SCRIPT_ROOT}" "$@"
  master_output="$(
    python3 "${SCRIPT_ROOT}/src/helpers/master_config.py" ribo \
      --master-config "${MASTER_CONFIG}"
  )"
  mapfile -t master_values <<< "${master_output}"

  ROOT="${master_values[0]}"
  MASTER_MANIFEST="${master_values[1]}"
  MASTER_RIBO_PREPROCESS_CONFIG="${master_values[2]}"
  MASTER_METAGENOME_PREPROCESS_CONFIG="${master_values[3]}"
  MASTER_CONTIG_ANNOTATION_CONFIG="${master_values[4]}"
  MASTER_LOG="${master_values[5]}"

  HELPER="${ROOT}/src/helpers/ribo_helper.py"
  RIBO_CONFIG="${RIBO_CONFIG:-${MASTER_RIBO_PREPROCESS_CONFIG}}"
  MANIFEST="${RIBO_MANIFEST:-${MASTER_MANIFEST}}"
  LOG="${RIBO_LOG:-${MASTER_LOG}}"

  metaribo_require_no_args "${PASSTHRU_ARGS[@]}"
  metaribo_start_log "${LOG}"

  cd "${ROOT}"

  echo "[$(date -Is)]: Starting ribo preprocessing"
  python3 "${HELPER}" check-paths \
    --manifest "${MANIFEST}" \
    --ribo-config "${RIBO_CONFIG}"

  config_shell="$(
    python3 "${HELPER}" config-shell --ribo-config "${RIBO_CONFIG}"
  )"
  eval "${config_shell}"
  mkdir -p \
    "${RIBO_DIR}" \
    "${BT2_DIR}" \
    "${REFERENCE_DIR}" \
    "${LOCAL_ALIGNMENT_DIR}" \
    "${LOCAL_STATS_DIR}"

  mapfile -t CUTADAPT_ARG_ARRAY < <(
    python3 "${HELPER}" split-args --args="${CUTADAPT_ARGS}"
  )
  mapfile -t UMITOOLS_EXTRACT_ARG_ARRAY < <(
    python3 "${HELPER}" split-args --args="${UMITOOLS_EXTRACT_ARGS}"
  )
  mapfile -t BT2_ARG_ARRAY < <(
    python3 "${HELPER}" split-args --args="${BT2_ARGS}"
  )
  mapfile -t UMITOOLS_DEDUP_ARG_ARRAY < <(
    python3 "${HELPER}" split-args --args="${UMITOOLS_DEDUP_ARGS}"
  )

  while IFS=$'\t' read -r ribo_id meta_id ribo_fq contigs ncrna_bed; do
    [[ "${ribo_id}" == "ribo_id" ]] && continue
    run_ribo_sample \
      "${ribo_id}" \
      "${meta_id}" \
      "${ribo_fq}" \
      "${contigs}" \
      "${ncrna_bed}"
  done < <(
    python3 "${HELPER}" list-inputs \
      --manifest "${MANIFEST}" \
      --ribo-config "${RIBO_CONFIG}"
  )

  echo "[$(date -Is)]: Ribo preprocessing complete"
}

main "$@"
