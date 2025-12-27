import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from py_chroma import (
    finalize,
    initialize,
    register_inline,
    run_inline_xml,
    set_gauge,
    set_lattice,
    set_rng_seed,
)


def load_input(path: Path):
    tree = ET.parse(path)
    root = tree.getroot()

    param = root.find("Param")
    if param is None:
        raise ValueError("Missing Param node in input XML.")

    nrow_text = (param.findtext("nrow") or "").strip()
    if not nrow_text:
        raise ValueError("Missing nrow in input XML.")
    nrow = [int(x) for x in nrow_text.split()]

    inline = param.find("InlineMeasurements")
    if inline is None:
        raise ValueError("Missing InlineMeasurements in input XML.")

    out_file = inline.find(".//NamedObject/out_file")
    if out_file is None:
        raise ValueError("Missing InlineMeasurements/elem/NamedObject/out_file.")
    out_file.text = "fred_py.xml"

    inline_xml = ET.tostring(inline, encoding="unicode")

    cfg = root.find("Cfg")
    if cfg is None:
        raise ValueError("Missing Cfg node in input XML.")
    cfg_type = (cfg.findtext("cfg_type") or "").strip()
    if not cfg_type:
        raise ValueError("Missing cfg_type in input XML.")

    cfg_xml = ET.tostring(cfg, encoding="unicode")
    return nrow, cfg_type, cfg_xml, inline_xml


def main() -> int:
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("plaq_density.ini.xml")
    nrow, cfg_type, cfg_xml, inline_xml = load_input(input_path)

    initialize()
    register_inline()
    set_lattice(nrow)
    set_rng_seed(11)
    set_gauge(cfg_type, cfg_xml, "/Cfg")
    run_inline_xml(inline_xml)
    finalize()
    print("py_chroma plaq_density ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
