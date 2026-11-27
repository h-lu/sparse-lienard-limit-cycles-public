"""Build the article and supplement twice in fresh trees and compare with the delivered pair."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import shutil
import tempfile
import tex_checks as checker

SOURCE_DIR = Path(__file__).resolve().parents[1]
PACKAGE_DIR = SOURCE_DIR.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="keep the rebuilt PDFs in a new directory outside the package")
    parser.add_argument("--compare-dir", type=Path,
                        help="compare with the two PDFs in this directory (default: the delivered pair, if present)")
    args = parser.parse_args()
    output = args.output_dir.resolve() if args.output_dir else None
    if output and (output.is_relative_to(PACKAGE_DIR) or output.exists()):
        parser.error("--output-dir must be a new directory outside the package")
    comparison = args.compare_dir.resolve() if args.compare_dir else PACKAGE_DIR
    with tempfile.TemporaryDirectory(prefix="lienard-pdf-build-") as tmp:
        directory = Path(tmp)
        pdfs = [job + ".pdf" for job in checker.JOBS]
        compare = args.compare_dir is not None or any((comparison / name).exists() for name in pdfs)
        if compare and not all((comparison / name).is_file() for name in pdfs):
            parser.error("the comparison directory must contain both PDFs")
        with ThreadPoolExecutor(max_workers=2) as pool:
            first, second = pool.map(checker.build, (directory / "first", directory / "second"))
        checker.inspected_links(first)
        different = False
        for name in pdfs:
            result = digest(first / name)
            if result != digest(second / name):
                raise RuntimeError(f"The two clean builds differ for {name}")
            print(f"PASS two clean builds, {name}, SHA-256 {result}")
            if compare:
                delivered = digest(comparison / name)
                if result != delivered:
                    different = True
                    print(f"DIFFER {name}: rebuilt {result}, delivered {delivered}")
                else:
                    print(f"MATCH delivered {name}, byte for byte")
        if output:
            output.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=output.parent) as pending:
                pair = Path(pending) / "pdfs"
                pair.mkdir()
                for name in pdfs:
                    shutil.copyfile(first / name, pair / name)
                pair.rename(output)
        if different:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
