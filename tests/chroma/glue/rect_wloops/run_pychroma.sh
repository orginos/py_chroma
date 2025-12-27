#!/usr/bin/env bash
set -euo pipefail

OMPI_MCA_btl=self,sm \
  PYCHROMA_POOL_GB="${PYCHROMA_POOL_GB:-0.05}" \
  PYCHROMA_LIB="${PYCHROMA_LIB:-/Users/kostas/Work/qcd_codes/chromaform/src/py_chroma/build/libpychroma.dylib}" \
  PYTHONPATH="${PYTHONPATH:-/Users/kostas/Work/qcd_codes/chromaform/src/py_chroma}" \
  python3 run_pychroma.py "${1:-rect_wloops.ini.xml}"
