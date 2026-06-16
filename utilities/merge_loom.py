# Merge multiple velocyto loom files.


import argparse
from pathlib import Path

import loompy


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge velocyto loom files.")
    parser.add_argument("loom_files", nargs="+", help="Input loom files")
    parser.add_argument("--output", default="D45_D50_D55.loom", help="Output loom file")
    parser.add_argument("--key", default="Accession", help="Merge key used by loompy.combine")
    args = parser.parse_args()

    files = [str(Path(f)) for f in args.loom_files]
    loompy.combine(files, args.output, key=args.key)
    print(f"Merged {len(files)} loom files into {args.output}")


if __name__ == "__main__":
    main()
