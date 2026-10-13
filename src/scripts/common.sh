# Common shell functions.

# Parse the shared --master-config option.
function metaribo_parse_master_config() {
  local script_root="$1"
  shift

  MASTER_CONFIG="${MASTER_CONFIG:-${script_root}/config/master.yaml}"

  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --master-config)
        if [[ "$#" -lt 2 ]]; then
          echo "ERROR: specify a path for --master-config." >&2
          exit 2
        fi
        MASTER_CONFIG="$2"
        shift 2
        ;;
      --master-config=*)
        MASTER_CONFIG="${1#*=}"
        shift
        ;;
      *)
        echo "Unknown argument: $1" >&2
        exit 2
        ;;
    esac
  done

  MASTER_CONFIG="$(realpath -- "${MASTER_CONFIG}")"
  export MASTER_CONFIG
}

function metaribo_start_log() {
  local log_path="$1"

  mkdir -p "$(dirname "${log_path}")"
  exec > >(tee -a "${log_path}") 2>&1
}
