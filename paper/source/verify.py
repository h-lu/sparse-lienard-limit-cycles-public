"""Run the article and supplement certificate invocations and compare their complete outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def certificate(item):
    args = item["argv"]
    run = subprocess.run([sys.executable, "-B", *args], cwd=ROOT,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    header = ("$ python3 -B " + " ".join(args) + "\n").encode()
    actual = header + run.stdout + b"\nExit status 0\n"
    reference = (ROOT / "checks" / (item["name"] + ".txt")).read_bytes()
    if run.stderr or actual != reference:
        raise RuntimeError(f'Raw output differs for {item["name"]}.')
    return f'PASS {item["name"]} (full output matches)'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificates-only", action="store_true",
                        help="run the certificate suite (the default)")
    parser.parse_args()
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    checks = json.loads((ROOT / "inputs/certificates.json").read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(certificate, checks):
            print(result)
    print(f"PASS {len(checks)} certificate transcripts.")


if __name__ == "__main__":
    main()
