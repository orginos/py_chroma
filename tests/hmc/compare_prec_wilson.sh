#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

CHROMA_OUT="${ROOT_DIR}/src/chroma/tests/hmc/hmc.prec_wilson.out.xml"
PY_OUT="${SCRIPT_DIR}/py_hmc.prec_wilson.out.xml"

python3 "${SCRIPT_DIR}/compare_hmc_outputs.py" "${CHROMA_OUT}" "${PY_OUT}"
