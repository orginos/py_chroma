import argparse
from pathlib import Path

import h5py
import numpy as np

from py_chroma import (
    finalize,
    initialize,
    random_gauge_transform,
    rect_wilson_loops,
    register_inline,
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
    parser = argparse.ArgumentParser(
        description="Random gauge transform test for rect_wloops."
    )
    parser.add_argument("--nrow", nargs=4, type=int, default=[4, 4, 4, 8])
    parser.add_argument("--t-dir", type=int, default=3)
    parser.add_argument("--z-dir", type=int, default=2)
    parser.add_argument("--l-max", type=int, default=3)
    parser.add_argument("--r-max", type=int, default=2)
    parser.add_argument("--out-1", default="rect_wloops_before.h5")
    parser.add_argument("--out-2", default="rect_wloops_after.h5")
    parser.add_argument("--cfg-type", default="WEAK_FIELD")
    parser.add_argument("--cfg-file", default="dummy")
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--tol-abs", type=float, default=1e-12)
    parser.add_argument("--tol-rel", type=float, default=1e-10)
    return parser.parse_args()


def load_dataset(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h5f:
        data = h5f["wilson_loops"][...]
    if data.dtype.fields and "r" in data.dtype.fields and "i" in data.dtype.fields:
        return data["r"] + 1j * data["i"]
    return data.astype(np.complex128)


def compare(a: np.ndarray, b: np.ndarray, tol_abs: float, tol_rel: float) -> bool:
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} != {b.shape}")
    diff = np.abs(a - b)
    max_abs = float(np.max(diff))
    denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), tol_abs)
    max_rel = float(np.max(diff / denom))
    ok = (max_abs <= tol_abs) or (max_rel <= tol_rel)
    print(f"max_abs={max_abs:.3e} max_rel={max_rel:.3e} -> {'OK' if ok else 'FAIL'}")
    return ok


def main() -> int:
    args = parse_args()
    cfg_xml = build_cfg_xml(args.cfg_type, args.cfg_file)

    initialize()
    register_inline()
    set_lattice(args.nrow)
    set_rng_seed(args.seed)
    set_gauge(args.cfg_type, cfg_xml, "/Cfg")

    rect_wilson_loops(args.t_dir, args.z_dir, args.l_max, args.r_max, args.out_1)
    random_gauge_transform()
    rect_wilson_loops(args.t_dir, args.z_dir, args.l_max, args.r_max, args.out_2)
    finalize()

    data_before = load_dataset(Path(args.out_1))
    data_after = load_dataset(Path(args.out_2))
    ok = compare(data_before, data_after, args.tol_abs, args.tol_rel)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
