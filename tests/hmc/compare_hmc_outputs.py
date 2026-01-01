import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def parse_bool(text):
    if text is None:
        return None
    return text.strip().lower() in ("true", "1", "yes")


def parse_float(text):
    if text is None:
        return None
    return float(text.strip())


def extract_updates(root):
    updates = []
    for update in root.findall(".//MCUpdates/elem/Update"):
        update_no = update.findtext("update_no")
        warmup = parse_bool(update.findtext("WarmUpP"))
        traj = update.find("HMCTrajectory")
        delta_h = parse_float(traj.findtext("deltaH")) if traj is not None else None
        acc_prob = parse_float(traj.findtext("AccProb")) if traj is not None else None
        accept_p = parse_bool(traj.findtext("AcceptP")) if traj is not None else None
        inline_tags = []
        inline_block = update.find("InlineObservables")
        if inline_block is not None:
            for elem in inline_block.findall("elem"):
                if list(elem):
                    inline_tags.append(list(elem)[0].tag)
        updates.append(
            {
                "update_no": int(update_no) if update_no is not None else None,
                "warmup": warmup,
                "delta_h": delta_h,
                "acc_prob": acc_prob,
                "accept_p": accept_p,
                "inline_tags": inline_tags,
            }
        )
    return updates


def float_close(a, b, atol, rtol):
    if a is None or b is None:
        return a == b
    return math.isclose(a, b, abs_tol=atol, rel_tol=rtol)


def compare_updates(left, right, atol, rtol):
    errors = []
    if len(left) != len(right):
        errors.append(f"Update count mismatch {len(left)} != {len(right)}")
        return errors

    for idx, (l_upd, r_upd) in enumerate(zip(left, right)):
        if l_upd["update_no"] != r_upd["update_no"]:
            errors.append(
                f"Update[{idx}] update_no mismatch {l_upd['update_no']} != {r_upd['update_no']}"
            )
        if l_upd["warmup"] != r_upd["warmup"]:
            errors.append(
                f"Update[{idx}] WarmUpP mismatch {l_upd['warmup']} != {r_upd['warmup']}"
            )
        if not float_close(l_upd["delta_h"], r_upd["delta_h"], atol, rtol):
            errors.append(
                f"Update[{idx}] deltaH mismatch {l_upd['delta_h']} != {r_upd['delta_h']}"
            )
        if not float_close(l_upd["acc_prob"], r_upd["acc_prob"], atol, rtol):
            errors.append(
                f"Update[{idx}] AccProb mismatch {l_upd['acc_prob']} != {r_upd['acc_prob']}"
            )
        if l_upd["accept_p"] != r_upd["accept_p"]:
            errors.append(
                f"Update[{idx}] AcceptP mismatch {l_upd['accept_p']} != {r_upd['accept_p']}"
            )
        if l_upd["inline_tags"] != r_upd["inline_tags"]:
            errors.append(
                f"Update[{idx}] Inline tags mismatch {l_upd['inline_tags']} != {r_upd['inline_tags']}"
            )
    return errors


def main():
    parser = argparse.ArgumentParser(description="Compare Chroma HMC XML outputs.")
    parser.add_argument("left", help="Left XML output")
    parser.add_argument("right", help="Right XML output")
    parser.add_argument("--atol", type=float, default=1e-2, help="Absolute tolerance")
    parser.add_argument("--rtol", type=float, default=1e-2, help="Relative tolerance")
    args = parser.parse_args()

    left_path = Path(args.left)
    right_path = Path(args.right)
    if not left_path.is_file():
        raise SystemExit(f"Left XML not found: {left_path}")
    if not right_path.is_file():
        raise SystemExit(f"Right XML not found: {right_path}")

    left_root = ET.parse(left_path).getroot()
    right_root = ET.parse(right_path).getroot()

    left_updates = extract_updates(left_root)
    right_updates = extract_updates(right_root)
    errors = compare_updates(left_updates, right_updates, args.atol, args.rtol)

    if errors:
        print("Mismatch:")
        for err in errors:
            print(f"  - {err}")
        raise SystemExit(1)

    print("OK: outputs match within tolerances.")


if __name__ == "__main__":
    main()
