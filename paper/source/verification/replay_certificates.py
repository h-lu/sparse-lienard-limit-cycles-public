"""Replay every archived invocation and retain the complete reproduced output."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def replay(directory, record):
    argv = record["argv"]
    run = subprocess.run([sys.executable, "-B", *argv], cwd=directory,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    header = ("$ python3 -B " + " ".join(argv) + "\n").encode()
    transcript = header + run.stdout + b"\nExit status 0\n"
    expected = (directory / "checks" / (record["name"] + ".txt")).read_bytes()
    if run.stderr or transcript != expected:
        raise ValueError(f"Raw output differs for {record['name']}.")
    return record["name"], transcript.decode("utf-8")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    directory = ROOT
    records = json.loads((directory / "inputs" / "certificates.json").read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        for name, output in pool.map(lambda r: replay(directory, r), records):
            print(f"===== {name} =====")
            print(output, flush=True)
    print(f"PASS {len(records)} fresh runs, all complete transcripts identical.")
