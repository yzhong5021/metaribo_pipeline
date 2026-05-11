#!/bin/bash
#
# Shared shell helpers for MetaRibo pipeline entrypoints.

function metaribo_die() {
  echo "ERROR: $*" >&2
  exit 2
}

function metaribo_parse_master_config() {
  local script_root="$1"
  shift

  MASTER_CONFIG="${MASTER_CONFIG:-${script_root}/config/master.yaml}"
  PASSTHRU_ARGS=()

  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --master-config)
        if [[ "$#" -lt 2 ]]; then
          metaribo_die "--master-config requires a path"
        fi
        MASTER_CONFIG="$2"
        shift 2
        ;;
      --master-config=*)
        MASTER_CONFIG="${1#*=}"
        shift
        ;;
      *)
        PASSTHRU_ARGS+=("$1")
        shift
        ;;
    esac
  done

  export MASTER_CONFIG
}

function metaribo_require_no_args() {
  if [[ "$#" -gt 0 ]]; then
    echo "Unknown argument: $1" >&2
    exit 2
  fi
}

function metaribo_start_log() {
  local log_path="$1"

  mkdir -p "$(dirname "${log_path}")"
  exec > >(tee -a "${log_path}") 2>&1
}
