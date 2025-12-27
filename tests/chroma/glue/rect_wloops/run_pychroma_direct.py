import argparse

from py_chroma import (
    finalize,
    initialize,
    register_inline,
    rect_wilson_loops,
    set_gauge,
    set_lattice,
    set_rng_seed,
)


def build_cfg_xml(cfg_type: str, cfg_file: str) -> str:
    return (
        "<Cfg>"
        f"<cfg_type>{cfg_type}</cfg_type>"
        f"<cfg_file>{cfg_file}</cfg_file>"
        "</Cfg>"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run rectangular Wilson loops via py_chroma.")
    parser.add_argument("--nrow", nargs=4, type=int, default=[4, 4, 4, 8])
    parser.add_argument("--t-dir", type=int, default=3)
    parser.add_argument("--z-dir", type=int, default=2)
    parser.add_argument("--l-max", type=int, default=3)
    parser.add_argument("--r-max", type=int, default=2)
    parser.add_argument("--out", default="rect_wloops_py.h5")
    parser.add_argument("--cfg-type", default="WEAK_FIELD")
    parser.add_argument("--cfg-file", default="dummy")
    parser.add_argument("--seed", type=int, default=11)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg_xml = build_cfg_xml(args.cfg_type, args.cfg_file)

    initialize()
    register_inline()
    set_lattice(args.nrow)
    set_rng_seed(args.seed)
    set_gauge(args.cfg_type, cfg_xml, "/Cfg")
    rect_wilson_loops(args.t_dir, args.z_dir, args.l_max, args.r_max, args.out)
    finalize()
    print("py_chroma rect_wloops ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
