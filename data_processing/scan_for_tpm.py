import os
import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan a directory tree for TPM/RNA/supplement files")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Root directory to scan (default: current working directory)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()

    matches = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            low = fn.lower()
            if (
                "supplement" in low
                or "tpm" in low
                or "rna" in low
                or "supplementary data 2" in low
                or "supplementary_data_2" in low
            ):
                matches.append(os.path.join(dirpath, fn))

    for m in matches:
        print(m)
    print("Done. Found", len(matches), "matches")


if __name__ == "__main__":
    main()
