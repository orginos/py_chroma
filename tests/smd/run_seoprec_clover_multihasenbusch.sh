#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

export OMPI_MCA_btl=self,sm
export PYCHROMA_POOL_GB=0.25
export PYCHROMA_LIB="${ROOT_DIR}/src/py_chroma/build/libpychroma.dylib"
export PYTHONPATH="${ROOT_DIR}/src/py_chroma"

INPUT_XML="${ROOT_DIR}/src/chroma/tests/smd/smd.seoprec_clover_multihasenbusch.ini.xml"
OUTPUT_XML="${SCRIPT_DIR}/py_smd.seoprec_clover_multihasenbusch.out.xml"

python3 "${SCRIPT_DIR}/run_smd_pychroma.py" "${INPUT_XML}" "${OUTPUT_XML}"
