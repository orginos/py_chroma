#!/usr/bin/env python3
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple


def _norm_text(text: Optional[str]) -> Optional[str]:
    if text is None:
        return None
    cleaned = " ".join(text.split())
    return cleaned or None


def _parse_float_tokens(text: str) -> Optional[Tuple[float, ...]]:
    tokens = text.split()
    if not tokens:
        return None
    values = []
    for tok in tokens:
        try:
            values.append(float(tok))
        except ValueError:
            return None
    return tuple(values)


def _compare_text(a: Optional[str], b: Optional[str], tol: float) -> bool:
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    a_norm = _norm_text(a) or ""
    b_norm = _norm_text(b) or ""
    if a_norm == b_norm:
        return True
    a_vals = _parse_float_tokens(a_norm)
    b_vals = _parse_float_tokens(b_norm)
    if a_vals is None or b_vals is None or len(a_vals) != len(b_vals):
        return False
    return all(abs(x - y) <= tol for x, y in zip(a_vals, b_vals))


SKIP_TAGS = {"run_date", "Source_file_info", "Prop_file_info"}


def _nodes_equal(a: ET.Element, b: ET.Element, tol: float) -> bool:
    if a.tag in SKIP_TAGS:
        return True
    if a.tag != b.tag:
        return False
    if sorted(a.attrib.items()) != sorted(b.attrib.items()):
        return False
    if not _compare_text(a.text, b.text, tol):
        return False
    if len(a) != len(b):
        return False
    return all(_nodes_equal(ca, cb, tol) for ca, cb in zip(list(a), list(b)))


def _first_diff(a: ET.Element, b: ET.Element, tol: float, path: str = "") -> Optional[str]:
    if a.tag in SKIP_TAGS:
        return None
    if a.tag != b.tag:
        return f"{path} tag {a.tag} != {b.tag}"
    if sorted(a.attrib.items()) != sorted(b.attrib.items()):
        return f"{path} attrib {a.attrib} != {b.attrib}"
    if not _compare_text(a.text, b.text, tol):
        return f"{path} text {(_norm_text(a.text) or '')} != {(_norm_text(b.text) or '')}"
    if len(a) != len(b):
        return f"{path} children {len(a)} != {len(b)}"
    for i, (ca, cb) in enumerate(zip(list(a), list(b))):
        diff = _first_diff(ca, cb, tol, f"{path}/{a.tag}[{i}]")
        if diff:
            return diff
    return None


def _find_all(root: ET.Element, tag: str) -> Iterable[ET.Element]:
    return root.findall(f".//{tag}")


def _gauge_diag(root: ET.Element) -> dict:
    diag = {}
    cfg = root.find(".//Input//Cfg")
    if cfg is None:
        cfg = root.find(".//Cfg")
    if cfg is None:
        return diag
    cfg_type = cfg.findtext("cfg_type")
    if cfg_type:
        diag["cfg_type"] = cfg_type.strip()
    gauge_state = cfg.findtext(".//GaugeState/Name")
    if gauge_state:
        diag["gauge_state"] = gauge_state.strip()
    gauge_bc = cfg.findtext(".//GaugeBC/Name")
    if gauge_bc:
        diag["gauge_bc"] = gauge_bc.strip()
    return diag


def _compare_nodes(label: str, a_nodes: Iterable[ET.Element], b_nodes: Iterable[ET.Element], tol: float) -> bool:
    a_list = list(a_nodes)
    b_list = list(b_nodes)
    ok = True
    if len(a_list) != len(b_list):
        print(f"{label}: count mismatch {len(a_list)} != {len(b_list)}")
        ok = False
    for i, (a_node, b_node) in enumerate(zip(a_list, b_list)):
        if not _nodes_equal(a_node, b_node, tol):
            diff = _first_diff(a_node, b_node, tol, f"{label}[{i}]")
            print(f"{label}: entry {i} differs")
            if diff:
                print(f"  {diff}")
            ok = False
            break
    if ok:
        print(f"{label}: match ({len(a_list)} sections)")
    return ok


def _inline_root(root: ET.Element) -> Optional[ET.Element]:
    if root.tag == "InlineObservables":
        return root
    return root.find(".//InlineObservables")


def _top_level_observables(root: ET.Element) -> Iterable[ET.Element]:
    return [child for child in list(root) if child.tag == "Observables"]


def main() -> int:
    if len(sys.argv) not in (3, 4):
        print("Usage: compare_outputs.py <chroma_output.xml> <pychroma_output.xml> [tol]")
        return 2

    chroma_path, py_path = sys.argv[1], sys.argv[2]
    chroma_root = ET.parse(chroma_path).getroot()
    py_root = ET.parse(py_path).getroot()
    tol = float(sys.argv[3]) if len(sys.argv) == 4 else 1e-10

    print("Gauge diagnostics:")
    print(f"  chroma:  {_gauge_diag(chroma_root) or 'n/a'}")
    print(f"  pychroma: {_gauge_diag(py_root) or 'n/a'}")

    ok = True
    chroma_inline = _inline_root(chroma_root)
    py_inline = _inline_root(py_root)

    ok &= _compare_nodes(
        "InlineObservables",
        [chroma_inline] if chroma_inline is not None else [],
        [py_inline] if py_inline is not None else [],
        tol,
    )

    ok &= _compare_nodes(
        "InlineObservables/Observables",
        _find_all(chroma_inline, "Observables") if chroma_inline is not None else [],
        _find_all(py_inline, "Observables") if py_inline is not None else [],
        tol,
    )

    ok &= _compare_nodes(
        "TopLevel/Observables",
        _top_level_observables(chroma_root),
        _top_level_observables(py_root),
        tol,
    )

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
