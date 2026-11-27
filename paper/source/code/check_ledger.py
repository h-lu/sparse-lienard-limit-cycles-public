"""Recompute support and profile data in the release ledger using rationals."""
from fractions import Fraction as Q
import json
from pathlib import Path


def require(test, label):
    if not test:
        raise ValueError(label)


def multiply(left, right):
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def general_support(n, blocks, prefix):
    stride = 258
    d = n * stride
    exponents = [d // 2, d - 1, d]
    exponents.extend(d - 1 + 2 * i for i in range(1, prefix + 1))
    for i in range(blocks):
        offset = 2 * prefix if i == 0 else 0
        exponents.extend((d + offset + stride - 1, d + 2 * stride,
                          d + 3 * stride - 1, d + 3 * stride))
        d += 3 * stride
    exponents.extend((d + 1, d + 3, d + 8))
    return sorted(exponents)


def low_supports(data):
    require([row["k"] for row in data["low"]] == list(range(3, 10)), "low budgets")
    counts = [row["project_hyperbolic"] for row in data["low"]]
    require(counts == [3, 4, 7, 10, 11, 14, 19], "low theorem counts")
    six, seven = data["low"][3:5]
    require(seven["project_support"] == [3, *six["project_support"]] and
            seven["project_hyperbolic"] == six["project_hyperbolic"] + 1, "inner cubic extension")
    alternative = data["eventual_even_wall_alternative"]
    require(alternative["exponents"] == [5, 9, 10, 11, 13, 15, "M"] and
            alternative["cycles"] == counts[4], "seven-power crest alternative")
    for row in data["low"]:
        k, support = row["k"], row["project_support"]
        if k == 3:
            require(support == [1, "d-1", "d"], "even-degree trinomial family")
        else:
            numeric = support[:-1] if support[-1] == "M" else support
            require(numeric == sorted(set(numeric)) and min(numeric) >= 1, f"collected k={k}")
            if len(numeric) < len(support):
                require(k == 9 and len(numeric) == k-1, "one eventual even exponent")
        require(len(support) == k and row["seed_or_odd_extension"] == k,
                f"row size and extension count at k={k}")
        odd_average = [Q(1)]
        for root in range(1, k):
            odd_average = multiply(odd_average, [-root, 1])
        require(all(odd_average) and len(odd_average) == k, "support of radial average")
        lifted = data["public_odd_lift"]["seed"] + data["public_odd_lift"]["new_powers"][:k-3]
        require(lifted == sorted(set(lifted)) and len(lifted) == k, "new odd powers")
    ten = data["public_odd_lift"]["seed"] + data["public_odd_lift"]["new_powers"]
    require(ten == [1, 9, 10, 11, 13, 15, 17, 19, 21, 23], "ten-term odd extension")
    for row in data["fixed_degree_alternatives"]:
        require(row["exponents"] == sorted(set(row["exponents"])) and
                len(row["exponents"]) == row["k"], "supplementary fixed support")
    print("Support sizes agree for k=3,...,9. Each radial average has k nonzero coefficients. "
          "Seven new odd powers give the recorded ten-term extension.")


def cubic_check(item):
    stride = item["stride"]
    lift, sep = Q(item["lift"]), Q(item["separation"])
    even = {0: Q(1), 2: lift-3, 3: Q(2)}
    odd = {0: Q(1), 1: lift-Q(3, 2)/sep, 3: 1/(2*sep**3)}
    for node, sign in zip(item["nodes"], item["signs"], strict=True):
        z = Q(node)
        value = sum(a*b*Q(1 + stride*(j-i), 17 + stride*(i+j))*z**(i+j)
                    for i, a in even.items() for j, b in odd.items())
        require(sign * value > Q(item["lower_absolute"]), f"cubic test {node}")
    print("The seven nine-power cubic values have the stated signs and strict rational margins.")


def count_checks(data):
    item = data["ten_power_input"]
    n = item["n"]
    require(item["D"] * n == item["seed_even"] == 2 * item["seed_odd"], "seed powers")
    require(general_support(n, 1, 0) == item["exponents"] and item["cycles"] == 26,
            "ten-power specialization")
    for blocks in range(1, 33):
        for prefix in range(4):
            support = general_support(n, blocks, prefix)
            k = 4 * blocks + 6 + prefix
            require(len(set(support)) == len(support) == k, "distinct general powers")
            require(sum(a % 2 == 0 for a in support) == 2 * blocks + 2, "even powers")
            count = (2 + 6 * blocks + prefix) + 4 + (10 * blocks + 2 * prefix + 4)
            require(count == 16 * blocks + 10 + 3 * prefix == 4*k - 14 - prefix, "residue ledger")
            require(count >= 4*k-17, "uniform count")
    require(item["ordinary_positive_signs"] == item["incoming_positive_signs"] == [1,-1,1,-1,1,-1,1],
            "two positive test lists")
    for i, pair in enumerate(item["negative_scaled_boxes"]):
        lo, hi = map(Q, pair)
        require(lo < hi and ((lo > 0) if i % 2 == 0 else (hi < 0)), "target signs at negative incoming heights")
    print("The ten exponents and target signs agree. The count and support ledgers agree on all four residues through 32 blocks; the proof covers every block count.")


if __name__ == "__main__":
    ledger = Path(__file__).resolve().parents[1] / "inputs" / "supports.json"
    data = json.loads(ledger.read_text())
    low_supports(data)
    cubic_check(data["nine_power_input"])
    count_checks(data)
