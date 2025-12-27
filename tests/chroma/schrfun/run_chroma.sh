#!/usr/bin/env bash
OMPI_MCA_btl=self,sm mpirun -np 1 "${CHROMA_BIN:-/Users/kostas/Work/qcd_codes/chromaform/install/chroma-qdpxx-double-nd4/bin/chroma}" -i "${1:?input.ini.xml}" -o "${2:-chroma.candidate.xml}"
