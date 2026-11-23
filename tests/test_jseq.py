#!/usr/bin/env python3
"""Tests with synthetic reads whose correct answer is known by construction.

Run:  python tests/test_jseq.py        (no pytest needed; exits non-zero on the first failure)
Takes a few seconds and reads no real data.
"""
import copy, csv, json, os, random, sys, tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import jseq  # noqa: E402

# The tests run against config/demo.json: a synthetic project with the same layout as a real
# experiment (same amplicon lengths, anchor and landmark offsets, thresholds), so they exercise the
# same code paths without carrying anyone's unpublished sequence. Regenerate with config/make_demo.py.
CFG = json.load(open(os.path.join(HERE, "..", "config", "demo.json")))
jseq.validate(CFG, [])
A5 = CFG["amplicons"]["Demo_5p"]
A3 = CFG["amplicons"]["Demo_3p"]
REF5 = A5["reference"]
CUT = 150            # the cut lies between reference index 149 and 150
L1_END = REF5.find(A5["landmark1"]["seq"]) + len(A5["landmark1"]["seq"])
DUP = L1_END - CUT   # size of a tandem duplication that copies the cut through landmark 1
rng = random.Random(7)


def tail(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def se300(seq):
    """Pad or truncate to a 300 nt single-end read."""
    return (seq + tail(300))[:300]


def mutate(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
    return "".join(s)


def classify(seq, amp=A5, qual=None):
    c = jseq.LandmarkCounter(amp)
    t = jseq.new_tile(c)
    c.add(t, seq, qual if qual is not None else "I" * len(seq))
    cls = [k for k in c.fields if k not in ("reads", "anchored", "tie_changed") and t[k]]
    return (cls[0] if cls else None), t


TESTS = []


def test(f):
    TESTS.append(f)
    return f


# ---------------------------------------------------------------- landmark mode

@test
def wild_type_is_net_zero():
    cls, t = classify(se300(REF5))
    assert cls == "zero", cls
    assert t["hist"] == Counter({0: 1})


@test
def errors_in_the_corrupted_stretch_do_not_matter():
    # 14 substitutions between the cut and landmark 1: the read is still the designed length
    cls, _ = classify(se300(mutate(REF5, range(CUT, 189, 3))))
    assert cls == "zero", cls


@test
def small_insertion_and_deletion_are_sized():
    cls, t = classify(se300(REF5[:CUT] + "GAT" + REF5[CUT:]))
    assert cls == "small" and t["hist"] == Counter({3: 1}), (cls, t["hist"])
    cls, t = classify(se300(REF5[:CUT] + REF5[CUT + 5:]))
    assert cls == "small" and t["hist"] == Counter({-5: 1}), (cls, t["hist"])


@test
def large_deletion_removing_landmark1_is_found_by_landmark2():
    # a deletion long enough to remove most of landmark 1
    cls, t = classify(se300(REF5[:CUT] + REF5[CUT + 84:]))
    assert cls == "large" and t["large_hist"] == Counter({-84: 1}), (cls, t["large_hist"])


@test
def tandem_duplication_is_not_called_net_zero():
    # Astra's counterexample: a tandem duplication of the cut through landmark 1. Landmark 1 then
    # sits at its reference distance with every position matching; only landmark 2 shows the read
    # is too long.
    dup = REF5[:L1_END] + REF5[CUT:L1_END] + REF5[L1_END:]
    cls, _ = classify(se300(dup))                 # 300 nt read: landmark 2 pushed off the end, but
    assert cls == "l2_contradicted", cls          # the well-read bases at its place are wrong
    cls, t = classify(dup + tail(60))             # long read: landmark 2 found DUP bp downstream
    assert cls == "l2_displaced" and t["displaced_hist"] == Counter({DUP: 1}), (cls, t)


@test
def truncated_read_cannot_hide_a_duplication():
    """Astra, 2026-11-20: a read that ends inside landmark 2 used to skip the check entirely, so a
    duplication with well-called mismatches in the part that WAS reached counted as net-zero."""
    dup = REF5[:L1_END] + REF5[CUT:L1_END] + REF5[L1_END:]
    cls, _ = classify(dup[:L1_END + 25])
    assert cls == "l2_contradicted", cls
    # a genuine net-zero read truncated at the same place is still only "unconfirmed", not wrong
    cls, _ = classify(se300(REF5)[:L1_END + 25])
    assert cls == "zero_l2_lowq", cls


@test
def poorly_read_landmark2_is_no_evidence_but_well_read_mismatch_is():
    # 12 mismatches in landmark 2 (needs 33 of 42). At quality 12 they say nothing; at quality 40
    # they contradict the net-zero call.
    pos = list(range(242, 282, 3))[:12]
    seq = se300(mutate(REF5, pos))
    lowq = "".join("-" if i in pos else "I" for i in range(len(seq)))
    cls, _ = classify(seq, qual=lowq)
    assert cls == "zero_l2_lowq", cls
    cls, _ = classify(seq)
    assert cls == "l2_contradicted", cls


@test
def confirm_policies_decide_what_counts():
    t = {"reads": 100, "anchored": 100, "zero": 60, "zero_l2_lowq": 5, "l2_displaced": 2,
         "l2_contradicted": 3, "small": 20, "large": 4, "unresolved": 6, "tie_changed": 0,
         "hist": {"0": 70, "3": 20}, "large_hist": {"-84": 4}, "displaced_hist": {"93": 2}}
    d = {"library": "x", "sample": "x", "culture": "1", "amplicon": "Demo_5p", "tiles": {"F:1:1101": t}}
    keep = {"F:1:1101": True}
    want = {"veto": (65, 6), "strict": (60, 11), "off": (70, 6)}     # (net-zero %, unresolved %)
    for pol, (nz, unres) in want.items():
        amp = dict(A5, landmark2=dict(A5["landmark2"], confirm=pol))
        row, spec, den = jseq.summarize_library(d, amp, keep, jseq.DEFAULT_QC)
        assert abs(row["pct_net_zero"] - nz) < 1e-9, (pol, row["pct_net_zero"])
        assert abs(row["pct_unresolved"] - unres) < 1e-9, (pol, row["pct_unresolved"])
        assert spec[0] == nz and spec[-84] == 4, (pol, spec)
        assert (spec.get(93, 0) == 2) == (pol != "off"), (pol, spec)
    row, _, _ = jseq.summarize_library(d, A5, keep, jseq.DEFAULT_QC)
    assert abs(row["pct_net_zero_ceiling"] - (100 - 20 - 4 - 5)) < 1e-9, row["pct_net_zero_ceiling"]


@test
def landmark2_tolerates_late_cycle_errors():
    # 8 mismatches in the 42 nt landmark 2 (threshold 33 matches) still confirm net-zero
    cls, _ = classify(se300(mutate(REF5, range(242, 282, 5))))
    assert cls == "zero", cls


@test
def anchor_rules_follow_cris_py():
    _, t = classify(se300(REF5[81:]))             # anchor at read position 0: CRIS.py ignores it
    assert t["anchored"] == 0
    _, t = classify(se300(REF5[100:]))            # no anchor
    assert t["anchored"] == 0
    _, t = classify(jseq.rc(se300(REF5)))         # reverse strand, strands = forward
    assert t["anchored"] == 0
    both = dict(A5, strands="both")
    cls, _ = classify(jseq.rc(se300(REF5)), both)
    assert cls == "zero", cls


@test
def unreadable_read_is_unresolved_not_dropped():
    junk = REF5[:101] + tail(199)                 # anchor present, nothing downstream readable
    cls, t = classify(junk)
    assert cls == "unresolved" and t["anchored"] == 1, cls


@test
def tie_breaking_prefers_the_smaller_length_change():
    # build landmark-1 hits of equal score at -4 and +4 is hard to do naturally; instead check the
    # scorer directly on a periodic sequence where every shift scores the same
    import numpy as np
    arr = np.frombuffer(("AC" * 100).encode(), np.uint8)
    lm = np.frombuffer(("AC" * 10).encode(), np.uint8)
    o, m, legacy = jseq.best_offset(arr, 50, np.arange(-10, 11), lm)
    assert (o, m, legacy) == (0, 20, -10), (o, m, legacy)


# ---------------------------------------------------------------- spacing mode (CRIS.py)

def cris_py_reference(reads, amp):
    """Independent re-statement of CRIS.py v2's loop (CRISpy_v2.py lines 205-221)."""
    ref, ss, se = amp["reference"], amp["anchor"], amp["end_anchor"]
    wt = ref.find(se) + len(se) - ref.find(ss)
    total, hist, test = 0, Counter(), 0
    for line in reads:
        if line.find(ss) > 0 and line.find(se) > 0:
            total += 1
            hist[line.find(se) + len(se) - line.find(ss) - wt] += 1
            if amp["designed_junction"] in line:
                test += 1
    return total, hist, test


@test
def spacing_mode_equals_cris_py_on_random_edits():
    ref = A3["reference"]
    reads = []
    for i in range(3000):
        s = ref
        k = rng.random()
        p = rng.randrange(len(ref))
        if k < 0.3:
            s = s[:p] + tail(rng.randint(1, 15)) + s[p:]
        elif k < 0.6:
            s = s[:p] + s[p + rng.randint(1, 30):]
        elif k < 0.7:
            s = mutate(s, [p])
        lead = tail(rng.choice([0, 0, 3, 10]))
        reads.append(lead + s + tail(rng.randint(0, 20)))
    c = jseq.SpacingCounter(A3)
    t = jseq.new_tile(c)
    for r in reads:
        c.add(t, r)
    total, hist, test = cris_py_reference(reads, A3)
    assert t["counted"] == total and t["hist"] == hist and t["designed_exact"] == test, \
        (t["counted"], total, t["designed_exact"], test)


# ---------------------------------------------------------------- tile rule

@test
def tile_selection_never_pools_lanes():
    """Astra, 2026-11-20: two libraries on different lanes of one flowcell, each measured
    perfectly but genuinely different, must not have one lane selected away by the other."""
    def lib(lane, pct):
        return {"tiles": {f"FC:{lane}:11{n:02d}": {"reads": 10000, "anchored": 10000,
                                                   "zero": pct * 100, "zero_l2_lowq": 0,
                                                   "l2_displaced": 0, "l2_contradicted": 0}
                          for n in range(1, 9)}}
    amp = {"mode": "landmark", "tile_filter": "plateau", "landmark2": {"confirm": "veto"}}
    libs = [lib(1, 80), lib(2, 20)]
    for lane in ("1", "2"):
        sub = [{**d, "tiles": {k: t for k, t in d["tiles"].items() if jseq.lane_of(k) == f"FC:{lane}"}}
               for d in libs]
        keep, _, _ = jseq.choose_tiles(amp, sub, jseq.DEFAULT_QC)
        kept = [k for k, v in keep.items() if v]
        assert len(kept) == 8, (lane, len(kept))     # every tile of that lane is uniform: keep all
    assert sorted({jseq.lane_of(k) for d in libs for k in d["tiles"]}) == ["FC:1", "FC:2"]


@test
def plateau_rule_keeps_tiles_near_the_best_ones():
    rates = {"1101": 80.4, "1102": 80.2, "1103": 80.1, "1104": 80.0, "1201": 79.9, "1202": 79.8,
             "1203": 79.7, "1204": 79.6, "1301": 79.3, "1302": 78.5, "1303": 60.0, "1304": 3.0}
    libs = [{"tiles": {f"FC:1:{k}": {"reads": 10000, "anchored": 10000, "zero": int(round(100 * v)),
                                     "zero_l2_lowq": 0, "l2_displaced": 0, "l2_contradicted": 0}
                       for k, v in rates.items()}}]
    amp = {"mode": "landmark", "tile_filter": "plateau", "landmark2": {"confirm": "veto"}}
    keep, _, plateau = jseq.choose_tiles(amp, libs, jseq.DEFAULT_QC)
    kept = sorted(k.split(":")[-1] for k, v in keep.items() if v)
    # best 8: 80.4 ... 79.6, median 79.95; keep >= 78.95
    assert abs(plateau - 79.95) < 1e-9, plateau
    assert kept == ["1101", "1102", "1103", "1104", "1201", "1202", "1203", "1204", "1301"], kept


# ---------------------------------------------------------------- end to end

def write_fastq(path, reads, tiles):
    with open(path, "w") as f:
        for i, (r, tl) in enumerate(zip(reads, tiles)):
            f.write(f"@INST:1:FCX:1:{tl}:{i}:1 1:N:0:A\n{r}\n+\n{'F' * len(r)}\n")


@test
def file_digest_covers_the_whole_of_a_small_file():
    """Astra, 2026-11-20 r3: the digest used to hash only the first megabyte of a file between
    1 and 2 MiB, so an edit at the end of such a FASTQ left the fingerprint unchanged."""
    d = tempfile.mkdtemp(prefix="jseq_digest_")
    for size in (1 << 19, (3 << 19), (1 << 20) + 10, 3 << 20):      # 0.5, 1.5, ~1, 3 MiB
        a = os.path.join(d, f"a{size}")
        with open(a, "wb") as f:
            f.write(b"A" * size)
        before = jseq.file_digest(a)
        with open(a, "r+b") as f:                                   # change the last byte only
            f.seek(-1, os.SEEK_END)
            f.write(b"C")
        assert jseq.file_digest(a) != before, f"end edit invisible at {size} bytes"
        assert os.path.getsize(a) == size


@test
def end_to_end_two_libraries_one_culture():
    d = tempfile.mkdtemp(prefix="jseq_test_")
    wt = se300(REF5)
    ins = se300(REF5[:CUT] + "GA" + REF5[CUT:])
    dele = se300(REF5[:CUT] + REF5[CUT + 84:])
    bad = se300(mutate(REF5, range(CUT, 240)))             # corrupted beyond landmark 1
    # good tiles 1101-1108 at 70 % net-zero; tile 1199 is a "bad" tile at 0 %
    reads, tiles = [], []
    for tl in [f"11{n:02d}" for n in range(1, 9)]:
        reads += [wt] * 70 + [ins] * 20 + [dele] * 10
        tiles += [tl] * 100
    reads += [bad] * 100
    tiles += ["1199"] * 100
    write_fastq(f"{d}/a.fastq", reads, tiles)
    write_fastq(f"{d}/b.fastq", reads, tiles)
    cfg = copy.deepcopy(CFG)
    cfg["amplicons"] = {"Demo_5p": cfg["amplicons"]["Demo_5p"]}
    cfg["qc"] = dict(cfg["qc"], min_reads_per_tile_pooled=50, min_reads_per_tile_library=50)
    cfg["heatmap"]["rows"] = [{"amplicon": "Demo_5p", "label": "5′"}]
    cfg["samples"] = "samples.tsv"
    json.dump(cfg, open(f"{d}/cfg.json", "w"))
    with open(f"{d}/samples.tsv", "w") as f:
        f.write("library\tsample\tculture\tamplicon\tr1\tr2\n")
        f.write(f"L1\ta\t1\tDemo_5p\t{d}/a.fastq\t\nL2\tb\t1\tDemo_5p\t{d}/b.fastq\t\n")
    for lib in ("L1", "L2"):
        jseq.main(["run-one", f"{d}/cfg.json", f"{d}/out", "--library", lib])
    jseq.main(["summarize", f"{d}/cfg.json", f"{d}/out"])
    rows = list(csv.DictReader(open(f"{d}/out/results/per_library.csv")))
    assert len(rows) == 2
    r = rows[0]
    assert r["tiles_kept"] == "8" and r["tiles_total"] == "9", (r["tiles_kept"], r["tiles_total"])
    assert abs(float(r["pct_net_zero"]) - 70.0) < 1e-9, r["pct_net_zero"]
    assert abs(float(r["pct_small_indel"]) - 20.0) < 1e-9
    assert abs(float(r["pct_large_deletion"]) - 10.0) < 1e-9
    assert r["tile_qc"] == "PASS"
    hm = list(csv.DictReader(open(f"{d}/out/results/heatmap_by_culture.csv")))
    assert abs(float(hm[0]["0"]) - 70.0) < 1e-6 and abs(float(hm[0]["<-10"]) - 10.0) < 1e-6
    assert abs(float(hm[0]["2"]) - 20.0) < 1e-6
    assert os.path.getsize(f"{d}/out/results/heatmap_by_culture.png") > 10000


if __name__ == "__main__":
    for f in TESTS:
        f()
        print(f"ok   {f.__name__}")
    print(f"\nall {len(TESTS)} tests passed")
