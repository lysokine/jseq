#!/usr/bin/env python3
"""Generate config/demo.json: a synthetic project with the same shape as a real one.

The test suite runs against this rather than against any real experiment's config, so the tests
carry no unpublished sequence or design and can be run by anyone. The dimensions below are
arbitrary; they are only chosen to be realistic enough that a test written against this config
exercises the same code paths as a two-junction amplicon experiment would.

Deterministic: run it again and you get the same file.
"""
import json, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(20261120)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def unique(ref, s):
    """A landmark has to occur exactly once, or the config check rejects it."""
    return ref.count(s) == 1


# --- landmark amplicon: 283 bp, anchor at 81, landmark 1 at 190 (50 nt), landmark 2 at 241 (42 nt)
while True:
    ref5 = seq(290)
    anchor5, l1, l2 = ref5[70:90], ref5[195:245], ref5[248:]
    if all(unique(ref5, x) for x in (anchor5, l1, l2)):
        break
# the cut this config describes sits at index 150
# --- spacing amplicon: 197 bp, anchors at 103 and 177, designed junction 40 nt
while True:
    ref3 = seq(225)
    a3, e3 = ref3[95:115], ref3[203:223]
    if unique(ref3, a3) and unique(ref3, e3) and ref3.find(e3) > ref3.find(a3):
        break
junction3 = ref3[150:190]

cfg = {
    "project": "Demo",
    "_about": ("A synthetic project used by the test suite. The sequences are random; the layout "
               "matches a real two-junction experiment. Regenerate with `python config/make_demo.py`."),
    "samples": "demo_samples.tsv",
    "qc": {"min_reads_per_tile_pooled": 1000, "min_reads_per_tile_library": 5000,
           "max_tile_spread_pts": 2.0, "plateau_top_n": 8, "plateau_tolerance_pts": 1.0},
    "amplicons": {
        "Demo_5p": {
            "description": "Synthetic landmark-mode amplicon, 290 bp, cut between index 149 and 150.",
            "reference": ref5, "mode": "landmark", "strands": "forward",
            "anchor": anchor5,
            "landmark1": {"seq": l1, "offsets": [-60, 60], "min_match": 40},
            "landmark2": {"seq": l2, "deletion_offsets": [-190, -56], "min_match": 33,
                          "confirm": "veto", "min_base_quality": 30},
            "tile_filter": "plateau"},
        "Demo_3p": {
            "description": "Synthetic spacing-mode amplicon, 225 bp.",
            "reference": ref3, "mode": "spacing", "strands": "forward",
            "anchor": a3, "end_anchor": e3, "designed_junction": junction3,
            "tile_filter": "none"}},
    "heatmap": {"label": "Demo",
                "rows": [{"amplicon": "Demo_5p", "label": "5′"}, {"amplicon": "Demo_3p", "label": "3′"}]},
}
with open(f"{HERE}/demo.json", "w") as f:
    json.dump(cfg, f, indent=1, ensure_ascii=False)
    f.write("\n")
# The sample sheet is written by the tests, which make their own FASTQ files; this one only
# documents the columns.
with open(f"{HERE}/demo_samples.tsv", "w") as f:
    f.write("# The tests write their own sample sheet. This file documents the columns.\n")
    f.write("library\tsample\tculture\tamplicon\tr1\tr2\n")
print(f"wrote {HERE}/demo.json  (cut index 150, landmark 1 at 195, landmark 2 at 248)")
