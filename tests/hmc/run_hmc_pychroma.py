import argparse
import os
from pathlib import Path

from py_chroma import finalize, initialize
from py_chroma import run_hmc_xml


def main():
    parser = argparse.ArgumentParser(description="Run HMC via py_chroma.")
    parser.add_argument("input_xml", help="Path to HMC XML input")
    parser.add_argument("output_xml", help="Path to output XML file")
    args = parser.parse_args()

    input_path = Path(args.input_xml)
    if not input_path.is_file():
        raise SystemExit(f"Input XML not found: {input_path}")

    output_path = Path(args.output_xml)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    os.environ["PYCHROMA_OUT"] = str(output_path)

    initialize()
    try:
        run_hmc_xml(input_path.read_text())
    finally:
        finalize()


if __name__ == "__main__":
    main()
