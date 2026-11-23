#!/bin/bash
# Check the package works outside its own project directory, on data it has never seen.
# Copies the package to a scratch directory, builds a small synthetic project there (one amplicon,
# paired-end so fastp is exercised), runs run_local.sh end to end, and checks the answer.
#
#   bash tests/portability.sh            (activate the junctionseq environment first)
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
T=$(mktemp -d "${JSEQ_TMP:-/tmp}/jseq_port_XXXX")
trap 'rm -rf "$T"' EXIT
echo "scratch: $T"
cp -r "$HERE" "$T/pkg"
rm -rf "$T/pkg/example_output" "$T/pkg/__pycache__" "$T/pkg/tests/__pycache__"

python - "$T" <<'PY'
import json, os, random, sys
T = sys.argv[1]
rng = random.Random(11)
tail = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
# a synthetic amplicon: random flanks around a designed junction
up, junc, down = tail(60), "GGACCTTAGCTAGGACCTTA", tail(60)
ref = up + junc + down
anchor, end_anchor = up[20:40], down[20:40]
cfg = {"project": "Portability", "samples": "s.tsv",
       "amplicons": {"amp": {"description": "synthetic", "reference": ref, "mode": "spacing",
                             "strands": "forward", "anchor": anchor, "end_anchor": end_anchor,
                             "designed_junction": junc, "tile_filter": "none"}},
       "heatmap": {"label": "Synthetic", "rows": [{"amplicon": "amp", "label": "5'"}]}}
d = f"{T}/proj"
os.makedirs(d)
json.dump(cfg, open(f"{d}/cfg.json", "w"))
# 600 reads: 400 designed length, 150 with +2, 50 with -6. Paired, overlapping by ~80 nt.
# edits land inside the designed junction (ref[60:80]), so they change its sequence as well
want = [(ref, 400), (ref[:70] + "GG" + ref[70:], 150), (ref[:70] + ref[76:], 50)]
r1, r2 = [], []
comp = str.maketrans("ACGT", "TGCA")
i = 0
for seq, n in want:
    for _ in range(n):
        s = tail(5) + seq + tail(5)
        r1.append(s[:120]); r2.append(s[-120:].translate(comp)[::-1]); i += 1
def write(path, reads):
    with open(path, "w") as f:
        for j, s in enumerate(reads):
            f.write(f"@M:1:FC:1:110{j%4+1}:{j}:1 1:N:0:A\n{s}\n+\n{'I'*len(s)}\n")
write(f"{d}/a_R1.fastq", r1); write(f"{d}/a_R2.fastq", r2)
with open(f"{d}/s.tsv", "w") as f:
    f.write("library\tsample\tculture\tamplicon\tr1\tr2\n")
    f.write(f"L1\ta\t1\tamp\t{d}/a_R1.fastq\t{d}/a_R2.fastq\n")
PY

cd "$T"                       # deliberately not inside the package or the project
bash "$T/pkg/run_local.sh" "$T/proj/cfg.json" "$T/out"

python - "$T" <<'PY'
import csv, sys
T = sys.argv[1]
r = list(csv.DictReader(open(f"{T}/out/results/per_library.csv")))[0]
assert r["denominator"] == "600", r["denominator"]
for k, want in (("pct_net_zero", 400/6), ("pct_designed_junction_exact", 400/6),  # 400 of 600
                ("pct_frame_preserving", 450/6), ("pct_large_indel_gt10bp", 0.0)):
    assert abs(float(r[k]) - want) < 1e-3, (k, r[k], want)   # CSV holds 4 decimals
spec = {int(x["net_length_change_bp"]): int(x["reads"]) for x in
        csv.DictReader(open(f"{T}/out/results/spectra.csv"))}
assert spec == {-6: 50, 0: 400, 2: 150}, spec
print("portability OK: merged, counted and summarized correctly outside the project")
PY
