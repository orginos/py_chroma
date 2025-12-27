#!/bin/bash
# need to set these external libs for the code to run
export PYCHROMA_LIB="${PYCHROMA_LIB:-/Users/kostas/Work/qcd_codes/chromaform/install/py_chroma-chroma-qdpxx-double-nd4-py_qdp/lib/libpychroma.dylib}" \
export PYTHONPATH="${PYTHONPATH:-/Users/kostas/Work/qcd_codes/chromaform/install/py_chroma-chroma-qdpxx-double-nd4-py_qdp/python}"

OMPI_MCA_btl=self,sm \
PYCHROMA_POOL_GB=0.05 \
python3 gauge_transform_test.py

