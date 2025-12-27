import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from py_chroma import (
    finalize,
    initialize,
    register_inline,
    rect_wilson_loops,
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

    rect = param.find("RectWLoops")
    if rect is None:
        raise ValueError("Missing RectWLoops in input XML.")

    def get_int(tag):
        text = (rect.findtext(tag) or "").strip()
        if not text:
            raise ValueError(f"Missing RectWLoops/{tag} in input XML.")
        return int(text)

    t_dir = get_int("t_dir")
    z_dir = get_int("z_dir")
    l_max = get_int("L_max")
    r_max = get_int("R_max")
    out_file = (rect.findtext("out_file") or "").strip()
    if not out_file:
        raise ValueError("Missing RectWLoops/out_file in input XML.")

    cfg = root.find("Cfg")
    if cfg is None:
        raise ValueError("Missing Cfg node in input XML.")
    cfg_type = (cfg.findtext("cfg_type") or "").strip()
    if not cfg_type:
        raise ValueError("Missing cfg_type in input XML.")
    cfg_xml = ET.tostring(cfg, encoding="unicode")

    return nrow, cfg_type, cfg_xml, t_dir, z_dir, l_max, r_max, out_file


def main() -> int:
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("rect_wloops.ini.xml")
    nrow, cfg_type, cfg_xml, t_dir, z_dir, l_max, r_max, out_file = load_input(input_path)

    initialize()
    register_inline()
    set_lattice(nrow)
    set_rng_seed(11)
    set_gauge(cfg_type, cfg_xml, "/Cfg")
    rect_wilson_loops(t_dir, z_dir, l_max, r_max, out_file)
    finalize()
    print("py_chroma rect_wloops ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
