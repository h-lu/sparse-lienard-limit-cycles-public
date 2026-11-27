"""Check the rational data printed in the Palomar preprint, Appendix A."""
import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def compare_source(path, data):
    body = path.read_bytes()
    assert hashlib.sha256(body).hexdigest() == data["source_sha256"]
    appendix = body.decode().split(r"\section{Data}\label{app:data}", 1)[1]
    coefficients = re.findall(r"^a_\{(\d+)\} &= (-?\d+/\d+)", appendix, re.M)
    assert coefficients == [(str(i), a) for i, a in enumerate(data["coefficients"], 1)]
    heights = re.findall(r"^h_\{?(\d+)\}? &= (\d+/\d+)", appendix, re.M)
    samples = data["samples"]
    assert heights == [(str(i), s["height"]) for i, s in enumerate(samples)]
    rows = re.findall(r"^\$A_(\d+)\^([+-])\$ & \\texttt\{([.\d]+)\} & \\texttt\{([.\d]+)\}",
                      appendix, re.M)
    expected = [(str(i), sign, *sample[side]) for i, sample in enumerate(samples)
                for side, sign in (("plus", "+"), ("minus", "-"))]
    assert rows == expected
    print("Pinned TeX digest and Appendix A agree: 6 coefficients, 7 heights, 28 endpoints.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-tex", type=Path,
                        help="also compare the frozen public TeX file with the archived transcription")
    args = parser.parse_args()
    data = json.loads((ROOT / "inputs/palomar-samples.json").read_text())
    if args.source_tex:
        compare_source(args.source_tex, data)
    coefficients = list(map(Fraction, data["coefficients"]))
    radius = Fraction(1, 10**40)
    assert len(coefficients) == 6 and all(abs(a) > radius for a in coefficients)
    print("All six coefficients stay nonzero throughout the printed sextic box.")

    previous_height = Fraction(0)
    signs = []
    margins = []
    for sample in data["samples"]:
        height = Fraction(sample["height"])
        assert height > previous_height
        previous_height = height
        plus_low, plus_high = map(Fraction, sample["plus"])
        minus_low, minus_high = map(Fraction, sample["minus"])
        assert 0 < plus_low <= plus_high and 0 < minus_low <= minus_high
        lower = plus_low - minus_high
        upper = plus_high - minus_low
        assert lower > 0 or upper < 0
        signs.append("+" if lower > 0 else "-")
        margins.append(min(abs(lower), abs(upper)))
    assert signs == ["+", "-", "+", "-", "+", "-", "+"]
    assert min(margins) == Fraction(2300260502, 10**42)
    print("Increasing heights, depth difference A+ - A-:", " ".join(signs))
    print("Least certified separation of the printed intervals: 2300260502 / 10^42.")

    brackets = len(signs) - 1
    counts = [brackets]
    for new_degree in range(7, 11):
        if new_degree % 2:
            brackets += 1
        counts.append(brackets)
    assert counts == [6, 7, 7, 8, 8]
    print("Brackets from the degree-raising rule at degrees 6 through 10:", *counts)

    support = list(range(1, len(coefficients) + 1))
    cycles = len(signs) - 1
    for power in (7, 9, 11, 13):
        assert power > support[-1] and power % 2 == 1
        support.append(power)
        cycles += 1
    assert len(support) == cycles == 10 and max(support) == 13
    print("Article Proposition A.1: four admissible odd extensions; 10 terms, degree 13, count 10.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()
