#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../../../.." && pwd)"

CHROMA_BIN="${CHROMA_BIN:-${ROOT_DIR}/install/chroma-qdpxx-double-nd4/bin/chroma}"
PYCHROMA_LIB="${PYCHROMA_LIB:-${ROOT_DIR}/src/py_chroma/build/libpychroma.dylib}"
PYCHROMA_POOL_GB="${PYCHROMA_POOL_GB:-0.05}"

cd "${SCRIPT_DIR}"
rm -f fred.xml fred_py.xml chroma.out.xml

OMPI_MCA_btl=self,sm mpirun -np 1 "${CHROMA_BIN}" -i plaq_density.ini.xml -o chroma.out.xml

OMPI_MCA_btl=self,sm PYCHROMA_POOL_GB="${PYCHROMA_POOL_GB}" \
  PYCHROMA_LIB="${PYCHROMA_LIB}" PYTHONPATH="${ROOT_DIR}/src/py_chroma" \
  python3 run_pychroma.py plaq_density.ini.xml

diff -u fred.xml fred_py.xml
echo "plaq_density match: fred.xml == fred_py.xml"
