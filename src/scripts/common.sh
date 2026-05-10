#!/bin/bash
#
# Shared shell helpers for MetaRibo pipeline entrypoints.

function metaribo_die() {
  echo "ERROR: $*" >&2
  exit 2
}

function metaribo_find_conda() {
  local conda_exe="${CONDA_EXE:-}"

  if [[ -z "${conda_exe}" ]]; then
    conda_exe="$(command -v conda || true)"
  fi
  if [[ -n "${conda_exe}" ]]; then
    export CONDA_EXE="${conda_exe}"
  fi
  printf '%s\n' "${conda_exe}"
}

function metaribo_reexec_in_conda() {
  local env_name="$1"
  local conda_exe
  shift

  conda_exe="$(metaribo_find_conda)"
  if [[ "${METARIBO_ACTIVE_SCRIPT_ENV:-}" == "${env_name}" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" ]]; then
      export PATH="${CONDA_PREFIX}/bin:${PATH}"
    fi
    return 0
  fi

  if [[ -z "${conda_exe}" ]]; then
    metaribo_die "conda is required to run $(basename "$0")"
  fi

  export METARIBO_ACTIVE_SCRIPT_ENV="${env_name}"
  exec "${conda_exe}" run --no-capture-output -n "${env_name}" "$0" "$@"
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
