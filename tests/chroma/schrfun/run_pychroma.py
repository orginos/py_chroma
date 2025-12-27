import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from py_chroma import (
    finalize,
    initialize,
    register_inline,
    run_inline_xml,
    set_gauge,
    set_input_xml,
    set_lattice,
    set_rng_seed,
    set_rng_seed_xml,
)


def load_input(path: Path):
    xml_text = path.read_text()
    tree = ET.fromstring(xml_text)
    root = tree

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
    inline_xml = ET.tostring(inline, encoding="unicode")

    rng = root.find("RNG")
    rng_xml = None
    if rng is not None:
        rng_xml = ET.tostring(rng, encoding="unicode")

    cfg = root.find("Cfg")
    if cfg is None:
        raise ValueError("Missing Cfg node in input XML.")
    cfg_type = (cfg.findtext("cfg_type") or "").strip()
    if not cfg_type:
        raise ValueError("Missing cfg_type in input XML.")
    cfg_xml = ET.tostring(cfg, encoding="unicode")

    return xml_text, nrow, rng_xml, cfg_type, cfg_xml, inline_xml


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit("Usage: run_pychroma.py <input.ini.xml> <output.xml>")

    input_path = Path(sys.argv[1])
    output_path = sys.argv[2]

    os.environ["PYCHROMA_OUT"] = output_path

    input_xml, nrow, rng_xml, cfg_type, cfg_xml, inline_xml = load_input(input_path)

    initialize()
    register_inline()
    set_input_xml(input_xml)
    set_lattice(nrow)
    if rng_xml:
        set_rng_seed_xml(rng_xml, "/RNG")
    else:
        set_rng_seed(11)
    set_gauge(cfg_type, cfg_xml, "/Cfg")
    run_inline_xml(inline_xml)
    finalize()
    print(f"py_chroma ok: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
