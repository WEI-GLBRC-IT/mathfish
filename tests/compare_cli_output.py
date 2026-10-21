"""Compare run_cli_tests.sh output against tests/expected_cli.out.

Text must match exactly; numbers must match to 1e-9 relative. The tolerance
absorbs last-digit drift between numpy builds: numpy 1.24 vs the 1.19 that
produced expected_cli.out moved 342 of 972 numbers, by at most 4.2e-14.

Usage: python3 compare_cli_output.py expected_cli.out tests.out
"""
import re
import sys

REL_TOL = 1e-9
NUMBER = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+")
CASE = re.compile(r"^### (\w+)\n(.*?)(?=^### |\Z)", re.S | re.M)


def cases(path):
    with open(path) as f:
        return {m.group(1): m.group(2) for m in CASE.finditer(f.read())}


def main(expected_path, actual_path):
    expected, actual = cases(expected_path), cases(actual_path)
    failures = []
    if list(expected) != list(actual):
        failures.append(f"cases differ: expected {list(expected)}, got {list(actual)}")
    compared = drifted = 0
    worst = 0.0
    for name in expected:
        if name not in actual:
            continue
        exp_text, act_text = expected[name], actual[name]
        if NUMBER.sub("#", exp_text) != NUMBER.sub("#", act_text):
            failures.append(f"[{name}] non-numeric text differs")
            continue
        for e, a in zip(map(float, NUMBER.findall(exp_text)), map(float, NUMBER.findall(act_text))):
            compared += 1
            if e == a:
                continue
            drifted += 1
            rel = abs(e - a) / max(abs(e), abs(a))
            worst = max(worst, rel)
            if rel > REL_TOL:
                failures.append(f"[{name}] {e!r} -> {a!r} (relative {rel:.1e})")
    print(f"cases: {len(actual)}/{len(expected)}, numbers compared: {compared}, "
          f"drifted within {REL_TOL:g}: {drifted - sum('relative' in f for f in failures)}, "
          f"worst relative difference: {worst:.1e}")
    for f in failures[:20]:
        print("FAIL", f)
    return 1 if failures or compared == 0 else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
