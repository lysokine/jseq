#!/usr/bin/env python3
"""Check a finished run against numbers produced earlier by other code.

  python tests/regression.py delivered <out dir> # vs the numbers delivered to the lab
  python tests/regression.py jph3      <out dir> # vs the 2023 CRIS.py output behind the JPH3 paper

Three checks against the delivered numbers:
  1. 3' junction: read counts and net-zero counts equal the delivered ones exactly.
  2. 5' junction, delivered settings (the 8 fixed tiles, landmark-2 check off): counts in
     every class equal scripts/final_v2.py exactly, apart from reads where two offsets tie for
     the best landmark-1 score (jseq takes the one nearest zero, final_v2 took the first).
  3. 5' junction, package defaults (plateau tile rule, landmark-2 veto): prints how far
     each library moves from the delivered value. This is a report, not a pass/fail.
For JPH3: Total, indel-0 and designed-sequence read counts equal CRIS.py's 2023 CSV exactly.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import jseq  # noqa: E402

WITHHELD = """The configuration, counts and expected values of the study project are not
distributed with the code: that study is not published yet. This script needs them. They will be
added when the paper is out; until then it can only be run inside the lab's copy."""


def need(path):
    if not os.path.exists(path):
        sys.exit(WITHHELD)
    return path




def counts(out, lib):
    with open(f"{out}/counts/{lib}.json") as f:
        return json.load(f)


def total(d, fields, tiles=None):
    return {k: sum(t[k] for key, t in d["tiles"].items()
                   if tiles is None or jseq.tile_number(key) in tiles) for k in fields}


def delivered(out):
    exp = json.load(open(need(f"{HERE}/expected/delivered.json")))
    legacy = set(exp["legacy_5p_tiles"])
    ok = True
    print("1. 3' junction, all tiles, vs delivered counts")
    for lib, e in exp["libraries"].items():
        if not lib.startswith("3J_"):
            continue
        d = counts(out, lib)
        got = total(d, ("counted",))
        zero = sum(t["hist"].get("0", 0) for t in d["tiles"].values())
        same = got["counted"] == e["counted"] and zero == e["zero"]
        ok &= same
        print(f"   {lib:<12s} counted {got['counted']:>9,} (delivered {e['counted']:>9,})  "
              f"net-zero {zero:>9,} ({e['zero']:>9,})  {'OK' if same else 'DIFFERENT'}")
    print("\n2. 5' junction, delivered settings, vs scripts/final_v2.py")
    for lib, e in exp["libraries"].items():
        if not lib.startswith("5J_"):
            continue
        d = counts(out, lib)
        g = total(d, jseq.LandmarkCounter.fields, legacy)
        mine = {"anchored": g["anchored"],
                "zero": g["zero"] + g["zero_l2_lowq"] + g["l2_displaced"] + g["l2_contradicted"],
                "small_indel": g["small"], "large_del": g["large"], "unresolved": g["unresolved"]}
        diff = {k: mine[k] - e[k] for k in mine if mine[k] != e[k]}
        same = mine["anchored"] == e["anchored"] and sum(abs(v) for v in diff.values()) <= 2 * g["tie_changed"]
        ok &= same
        print(f"   {lib:<12s} anchored {mine['anchored']:>9,}  net-zero {mine['zero']:>9,} "
              f"(final_v2 {e['zero']:>9,})  tie reads {g['tie_changed']:>3d}  "
              f"{'OK' if same else 'DIFFERENT ' + json.dumps(diff)}")
    print("\n3. 5' junction, package defaults, vs delivered % (report only)")
    rows = {r["library"]: r for r in csv.DictReader(open(f"{out}/results/per_library.csv"))}
    print(f"   {'library':<12s} {'delivered':>9s} {'now':>8s} {'change':>7s}   "
          f"{'L2 displaced':>12s} {'L2 contradicted':>15s} {'L2 low quality':>14s}  tiles")
    for lib, e in exp["libraries"].items():
        if not lib.startswith("5J_"):
            continue
        r = rows[lib]
        now = float(r["pct_net_zero"])
        print(f"   {lib:<12s} {e['pct_net_zero_rounded']:>8.2f}% {now:>7.2f}% {now - e['pct_net_zero_rounded']:>+7.2f}"
              f"   {float(r['pct_landmark2_displaced']):>11.3f}% {float(r['pct_landmark2_contradicted']):>14.3f}% "
              f"{float(r['pct_landmark2_low_quality']):>13.3f}%  "
              f"{r['tiles_kept']}/{r['tiles_total']}")
    return ok


def jph3(out):
    ok = True
    print("JPH3, all tiles, vs CRIS.py 2023 CSV")
    for e in csv.DictReader(open(f"{HERE}/expected/jph3_cris_py_2023.csv")):
        d = counts(out, e["library"])
        g = total(d, ("counted", "designed_exact"))
        zero = sum(t["hist"].get("0", 0) for t in d["tiles"].values())
        mine = (g["counted"], zero, g["designed_exact"])
        theirs = (int(e["Total"]), int(e["indel0_reads"]), int(e["designed_reads"]))
        same = mine == theirs
        ok &= same
        print(f"   {e['library']:<8s} Total {mine[0]:>9,} ({theirs[0]:>9,})  indel 0 {mine[1]:>9,} "
              f"({theirs[1]:>9,}) = {100 * mine[1] / mine[0]:.1f}%  designed {mine[2]:>9,} ({theirs[2]:>9,})  "
              f"{'OK' if same else 'DIFFERENT'}")
    return ok


if __name__ == "__main__":
    which, out = sys.argv[1], sys.argv[2]
    ok = {"delivered": delivered, "jph3": jph3}[which](out)
    print("\nREGRESSION", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)
