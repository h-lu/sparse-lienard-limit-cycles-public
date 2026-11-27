"""Recreate both scientific archives from the editable extraction."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

SOURCE = Path(__file__).resolve().parents[1]


def manifest(files):
    return "".join(f"{hashlib.sha256(body).hexdigest()}  {name}\n"
                   for name, body in sorted(files.items())).encode()


def write_archive(path, files):
    with ZipFile(path, "x") as archive:
        for name, body in sorted({**files, "SHA256SUMS": manifest(files)}.items()):
            entry = ZipInfo(name, (2026, 11, 23, 0, 0, 0))
            entry.create_system = 3
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, body, compress_type=ZIP_DEFLATED, compresslevel=9)


def pack(output):
    # Use the distributed inventory. Local build files and unrelated documents stay outside it.
    files = {}
    for line in (SOURCE.parent / "SHA256SUMS").read_text().splitlines():
        _, name = line.split("  ", 1)
        target = (SOURCE.parent / name).resolve()
        if not name.startswith("source/") or not target.is_relative_to(SOURCE) or name in files:
            raise ValueError(f"Unexpected path in source inventory: {name}")
        files[name] = target.read_bytes()
    aliases = {"README.txt": "COMPUTATION-README.txt", "requirements.txt": "requirements-numerical.txt"}
    replay = set(["verify.py", "inputs/certificates.json", "inputs/supports.json",
                  "inputs/degree-ten-seed.json", "verification/certificate-summary.txt",
                  "verification/crosscheck.py", "verification/crosscheck.txt",
                  "verification/positive-valleys-direct.txt", "verification/mechanism-replay.txt",
                  "verification/refined-quadrature.txt", "verification/check_stride_two.py",
                  "verification/stride-two-fixed.txt", "inputs/palomar-samples.json",
                  "verification/check_full_support.py", "verification/full-support.txt", *aliases])
    for entry in json.loads(files["source/inputs/certificates.json"]):
        replay.update([entry["argv"][0], "checks/" + entry["name"] + ".txt"])
    numerical = {name: files["source/" + aliases.get(name, name)] for name in replay}
    output.mkdir(parents=True, exist_ok=False)
    write_archive(output / "sparse-lienard-editable-source.zip", files)
    write_archive(output / "sparse-lienard-replay.zip", numerical)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="new output directory for the two archives")
    pack(parser.parse_args().directory)
