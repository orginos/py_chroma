# py_chroma

Minimal Python entry point for running Chroma inline measurements without
writing XML files to disk. This first cut still accepts XML strings so we can
drive Chroma from Python and then iterate toward typed parameter APIs.

## Build

```
mkdir -p build
cmake -S . -B build \
  -DCHROMA_PREFIX=/Users/kostas/Work/qcd_codes/chromaform/install/chroma-qdpxx-double-nd4-superbblas-cpu \
  -DQDP_PREFIX=/Users/kostas/Work/qcd_codes/chromaform/install/qdpxx-double-nd4
cmake --build build
```

The shared library will be written to `build/libpychroma.dylib` on macOS.

## Python Usage (ctypes)

```
export PYCHROMA_LIB=/Users/kostas/Work/qcd_codes/chromaform/src/py_chroma/build/libpychroma.dylib
python - <<'PY'
from py_chroma import (
    initialize, finalize, register_inline, set_lattice, set_rng_seed,
    set_gauge, run_inline_xml, run_plaquette,
)

initialize()
register_inline()
set_lattice([4, 4, 4, 8])
set_rng_seed(11)

# cfg_xml and inline_xml are XML strings that Chroma already accepts.
# Replace these with your XML strings or generated snippets.
cfg_xml = "<Cfg><cfg_type>...etc...</cfg_type></Cfg>"
set_gauge("cfg_type_here", cfg_xml, "/Cfg")

inline_xml = "<InlineMeasurements>...</InlineMeasurements>"
run_inline_xml(inline_xml)
run_plaquette(update_no=0, frequency=1)
finalize()
PY
```

Next step is to replace XML strings with Python-side parameter objects that
map directly onto Chroma C++ param structs.
