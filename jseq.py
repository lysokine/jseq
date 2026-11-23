#!/usr/bin/env python3
"""jseq: measure correct integration at knock-in junctions from amplicon sequencing.

Each junction (5' or 3') is one PCR amplicon. For every library the program counts how many reads
carry the junction at its designed length ("net-zero"), and the full spectrum of length changes.
Two counting modes are available, chosen per amplicon in the project config:

  spacing   CRIS.py v2 counting (Connelly and Pruett-Miller 2019). A read counts when both the start
            anchor and the end anchor occur in it exactly. Its length change is the anchor spacing in
            the read minus the anchor spacing in the reference. Denominator: reads with both anchors.

  landmark  For junctions whose own sequence corrupts the base calls downstream of it, such as a
            designed junction that is a perfect inverted repeat. One exact upstream anchor,
            then two error-tolerant landmarks downstream of the corrupted stretch:
              landmark 1  searched at offsets near its reference distance (small indels);
              landmark 2  at the end of the amplicon, searched further back (large deletions that
                          remove landmark 1), and used to confirm every net-zero call.
            Denominator: every read that carries the upstream anchor, so no read is silently dropped.

Subcommands
  run-one       merge (if paired) and count one library            -> OUT/counts/<library>.json
  summarize     tile QC, per-library and per-culture tables, heatmap -> OUT/results/
  list-samples  print the libraries in a project, one per line
  check         validate a project config and its sample sheet without reading any FASTQ

Run `python jseq.py <subcommand> -h` for arguments. See README.md for the full description.
"""
import argparse, csv, gzip, hashlib, json, math, os, shutil, statistics, subprocess, sys, tempfile
from collections import Counter, defaultdict

try:
    import numpy as np
except ImportError:      # list-samples and check work without numpy, e.g. on a login node
    np = None

VERSION = "1.1.0"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


# --------------------------------------------------------------------------- configuration

DEFAULT_QC = {"min_reads_per_tile_pooled": 1000, "min_reads_per_tile_library": 5000,
              "max_tile_spread_pts": 2.0, "plateau_top_n": 8, "plateau_tolerance_pts": 1.0}
DEFAULT_MERGE = {"overlap_len_require": 10, "overlap_diff_limit": 8,
                 "overlap_diff_percent_limit": 25, "threads": 4,
                 "trim_poly_g": True, "poly_g_min_len": 10}


def read_samples(path):
    """Relative FASTQ paths are resolved against the sample sheet, so a sheet can be moved with
    its data and run from any directory."""
    base = os.path.dirname(os.path.abspath(path))
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    need = {"library", "sample", "culture", "amplicon", "r1"}
    if not rows:
        sys.exit(f"{path}: no samples")
    miss = need - set(rows[0])
    if miss:
        sys.exit(f"{path}: missing column(s) {sorted(miss)}")
    bad = []
    for r in rows:
        r["r2"] = (r.get("r2") or "").strip()
        for k in need:
            if r.get(k) is None:
                sys.exit(f"{path}: a row is missing the {k} column")
            r[k] = r[k].strip()
        for k in ("r1", "r2"):
            if r[k] and not os.path.isabs(r[k]):
                r[k] = os.path.normpath(os.path.join(base, r[k]))
        if not r["library"] or any(c.isspace() for c in r["library"]) or "/" in r["library"]:
            bad.append(r["library"])
    if bad:
        sys.exit(f"{path}: library names must have no spaces or slashes, they become file names: {bad}")
    libs = [r["library"] for r in rows]
    dup = [l for l, n in Counter(libs).items() if n > 1]
    if dup:
        sys.exit(f"{path}: duplicate library names {dup}")
    return rows


def load_project(cfg_path):
    with open(cfg_path) as f:
        cfg = json.load(f)
    base = os.path.dirname(os.path.abspath(cfg_path))
    sp = cfg["samples"]
    samples = read_samples(sp if os.path.isabs(sp) else os.path.join(base, sp))
    cfg["qc"] = {**DEFAULT_QC, **cfg.get("qc", {})}
    cfg["merge"] = {**DEFAULT_MERGE, **cfg.get("merge", {})}
    problems = validate(cfg, samples)
    if problems:
        sys.exit("config problems:\n  " + "\n  ".join(problems))
    return cfg, samples


def validate(cfg, samples):
    out = []
    for name, a in cfg["amplicons"].items():
        ref = a["reference"].upper()
        a["reference"] = ref
        for k in ("anchor", "end_anchor", "designed_junction"):
            if a.get(k):
                a[k] = a[k].upper()
        if a["mode"] not in ("spacing", "landmark"):
            out.append(f"{name}: mode must be 'spacing' or 'landmark'")
            continue
        if a.get("strands", "forward") not in ("forward", "reverse", "both"):
            out.append(f"{name}: strands must be forward, reverse or both")
        if ref.count(a["anchor"]) != 1:
            out.append(f"{name}: anchor must occur exactly once in the reference "
                       f"(found {ref.count(a['anchor'])})")
        if a["mode"] == "spacing":
            if ref.count(a["end_anchor"]) != 1:
                out.append(f"{name}: end_anchor must occur exactly once in the reference")
            elif ref.find(a["end_anchor"]) <= ref.find(a["anchor"]):
                out.append(f"{name}: end_anchor must lie downstream of anchor")
            if a.get("designed_junction") and a["designed_junction"] not in ref:
                out.append(f"{name}: designed_junction is not in the reference")
        else:
            for k in ("landmark1", "landmark2"):
                s = a[k]["seq"].upper()
                a[k]["seq"] = s
                if ref.count(s) != 1:
                    out.append(f"{name}: {k} must occur exactly once in the reference")
                elif ref.find(s) <= ref.find(a["anchor"]):
                    out.append(f"{name}: {k} must lie downstream of anchor")
                if not 0 < a[k]["min_match"] <= len(s):
                    out.append(f"{name}: {k}.min_match must be between 1 and {len(s)}")
            if a["landmark2"].get("confirm", "veto") not in ("veto", "strict", "off"):
                out.append(f"{name}: landmark2.confirm must be veto, strict or off")
        tf = a.get("tile_filter", "none")
        if not (tf in ("none", "plateau") or isinstance(tf, list)):
            out.append(f"{name}: tile_filter must be 'none', 'plateau' or a list of tile numbers")
    for r in samples:
        if r["amplicon"] not in cfg["amplicons"]:
            out.append(f"{r['library']}: amplicon '{r['amplicon']}' is not defined in the config")
    return out


# --------------------------------------------------------------------------- reading FASTQ

def fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        n = 0
        while True:
            h = f.readline()
            if not h:
                return
            s, plus, q = f.readline(), f.readline(), f.readline()
            n += 1
            if not h.startswith("@"):
                sys.exit(f"{path}: record {n} does not start with '@'. This is not FASTQ, or the "
                         f"file is misaligned.")
            if not q or not plus.startswith("+"):
                sys.exit(f"{path}: record {n} is truncated or malformed. The file is incomplete; "
                         f"counting it would silently undercount this library.")
            s, q = s.rstrip(), q.rstrip()
            if len(s) != len(q):
                sys.exit(f"{path}: record {n} has {len(s)} bases and {len(q)} quality values")
            yield h, s, q


def tile_of(header):
    """Illumina read names are instrument:run:flowcell:lane:tile:x:y. The tile key keeps the
    flowcell and lane so libraries from different runs are never pooled together."""
    f = header[1:].split(None, 1)[0].split(":")
    return f"{f[2]}:{f[3]}:{f[4]}" if len(f) >= 7 else "NA"


def lane_of(tile_key):
    """The group a tile belongs to for tile selection: one flowcell lane. Lanes are kept apart
    because libraries are not always spread evenly over them, and a selection rule that pooled
    lanes could pick the lane carrying the better-performing libraries."""
    return tile_key.rsplit(":", 1)[0] if tile_key != "NA" else "NA"


def strand_views(s, strands, q=""):
    """(sequence, quality) pairs to search, in order; quality is reversed with the sequence."""
    if strands == "forward":
        return ((s, q),)
    if strands == "reverse":
        return ((rc(s), q[::-1]),)
    return ((s, q), (rc(s), q[::-1]))


# --------------------------------------------------------------------------- counting modes

class SpacingCounter:
    """CRIS.py v2 logic. `find(...) > 0` is deliberate: CRIS.py ignores an anchor that starts at
    read position 0, and keeping that makes the counts identical to CRIS.py output."""
    fields = ("reads", "counted", "designed_exact")

    def __init__(self, amp):
        ref = amp["reference"]
        self.ss, self.se = amp["anchor"], amp["end_anchor"]
        self.wt = ref.find(self.se) + len(self.se) - ref.find(self.ss)
        self.test = amp.get("designed_junction") or None
        self.strands = amp.get("strands", "forward")

    def add(self, t, s, q=""):
        t["reads"] += 1
        for v, _ in strand_views(s, self.strands):
            a = v.find(self.ss)
            b = v.find(self.se)
            if a > 0 and b > 0:
                t["counted"] += 1
                t["hist"][b + len(self.se) - a - self.wt] += 1
                if self.test and self.test in v:
                    t["designed_exact"] += 1
                return


def best_offset(arr, base, offs, lm):
    """Score landmark `lm` (uint8 array) at read position base+o for every o in offs that fits.
    Returns (offset, matches, legacy_offset). Ties are broken towards the offset nearest zero
    (the smallest length change that explains the read); legacy_offset is the first maximum,
    which is what the delivered November 2026 numbers used."""
    n = len(lm)
    o = offs[(offs + base >= 0) & (offs + base + n <= len(arr))]
    if not len(o):
        return None, 0, None
    idx = (base + o)[:, None] + np.arange(n)[None, :]
    m = (arr[idx] == lm[None, :]).sum(1)
    mx = int(m.max())
    cand = o[m == mx]
    return int(cand[np.argmin(np.abs(cand))]), mx, int(o[int(m.argmax())])


def n_match(a, b):
    return sum(x == y for x, y in zip(a, b))


class LandmarkCounter:
    """Classes of anchored reads (see README, 'landmark'):
      zero              landmark 1 at net zero, landmark 2 confirms it
      zero_l2_lowq      landmark 1 at net zero, landmark 2 too poorly read to check either way
      l2_displaced      landmark 1 at net zero, landmark 2 found at another distance
      l2_contradicted   landmark 1 at net zero, high-quality bases at landmark 2's place disagree
      small             landmark 1 found at a non-zero offset (small indel)
      large             landmark 1 not found, landmark 2 found further upstream (large deletion)
      unresolved        neither landmark found
    Which of the zero* / l2_* classes count as net-zero is decided in `summarize`."""
    fields = ("reads", "anchored", "zero", "zero_l2_lowq", "l2_displaced", "l2_contradicted",
              "small", "large", "unresolved", "tie_changed")
    hists = ("hist", "large_hist", "displaced_hist")

    def __init__(self, amp):
        ref = amp["reference"]
        a0 = ref.find(amp["anchor"])
        self.anchor = amp["anchor"]
        L1, L2 = amp["landmark1"], amp["landmark2"]
        self.l1, self.l2 = L1["seq"], L2["seq"]
        self.n1, self.n2 = len(self.l1), len(self.l2)
        self.e1, self.e2 = ref.find(self.l1) - a0, ref.find(self.l2) - a0
        self.t1, self.t2 = L1["min_match"], L2["min_match"]
        self.o1 = np.arange(L1["offsets"][0], L1["offsets"][1] + 1)
        self.o2 = np.arange(L2["deletion_offsets"][0], L2["deletion_offsets"][1] + 1)
        self.ow = np.arange(L2["deletion_offsets"][0], 2001)   # where landmark 2 went, if not home
        self.b1 = np.frombuffer(self.l1.encode(), np.uint8)
        self.b2 = np.frombuffer(self.l2.encode(), np.uint8)
        self.hq = chr(33 + int(L2.get("min_base_quality", 30)))   # Phred+33 character
        self.strands = amp.get("strands", "forward")

    def add(self, t, s, q=""):
        t["reads"] += 1
        for v, vq in strand_views(s, self.strands, q):
            a = v.find(self.anchor)
            if a > 0:          # same convention as CRIS.py: an anchor at position 0 is not used
                break
        else:
            return
        t["anchored"] += 1
        p1 = a + self.e1
        arr = None
        if v[p1:p1 + self.n1] == self.l1:               # fast path: exact, at the reference distance
            o1, m1 = 0, self.n1
            lo = max(0, p1 - (-int(self.o1[0])))
            if v.find(self.l1, lo, p1 + self.n1 - 1) not in (-1, p1):
                t["tie_changed"] += 1                    # a second exact copy upstream (legacy pick)
        else:
            arr = np.frombuffer(v.encode(), np.uint8)
            o1, m1, legacy = best_offset(arr, p1, self.o1, self.b1)
            if o1 is not None and m1 >= self.t1 and legacy != o1:
                t["tie_changed"] += 1
        if o1 is not None and m1 >= self.t1:
            t["hist"][o1] += 1
            if o1 != 0:
                t["small"] += 1
                return
            p2 = a + self.e2
            home = v[max(p2, 0):p2 + self.n2]            # may be short: the read can end inside it
            if len(home) == self.n2 and n_match(home, self.l2) >= self.t2:
                t["zero"] += 1
                return
            if arr is None:
                arr = np.frombuffer(v.encode(), np.uint8)
            o2, m2, _ = best_offset(arr, p2, self.ow, self.b2)
            if o2 is not None and m2 >= self.t2:
                t["l2_displaced"] += 1                   # landmark 2 is elsewhere: rearranged
                t["displaced_hist"][o2] += 1
                return
            # A partial window still carries evidence: more well-read disagreements than the whole
            # landmark is allowed means this cannot be landmark 2, however much of it was reached.
            hq_mismatch = sum(1 for x, y, c in zip(home, self.l2[max(-p2, 0):], vq[max(p2, 0):p2 + self.n2])
                              if x != y and c >= self.hq)
            if hq_mismatch > self.n2 - self.t2:
                t["l2_contradicted"] += 1                # well-read bases say it is not landmark 2
            else:
                t["zero_l2_lowq"] += 1                   # not enough evidence either way
            return
        if arr is None:
            arr = np.frombuffer(v.encode(), np.uint8)
        o2, m2, _ = best_offset(arr, a + self.e2, self.o2, self.b2)
        if o2 is not None and m2 >= self.t2:
            t["large"] += 1
            t["large_hist"][o2] += 1
        else:
            t["unresolved"] += 1


COUNTING_KEYS = ("reference", "mode", "strands", "anchor", "end_anchor", "designed_junction",
                 "landmark1", "landmark2")


def counting_fingerprint(amp, merge, merged):
    """Identifies everything that changes the counts. Summary-only settings (tile_filter, qc,
    heatmap) are deliberately left out, so those can be changed and re-summarized without
    recounting; landmark2.confirm is a summary-time choice too. `merged` says whether this library
    goes through fastp, which depends on the sample sheet rather than the amplicon's mode."""
    a = {k: amp[k] for k in COUNTING_KEYS if k in amp}
    if "landmark2" in a:
        a["landmark2"] = {k: v for k, v in a["landmark2"].items() if k != "confirm"}
    a["merge"] = merge if merged else None
    return hashlib.sha256(json.dumps(a, sort_keys=True).encode()).hexdigest()[:16]


def file_digest(path, window=1 << 20):
    """Cheap identity check for a FASTQ: a digest of the whole file up to 2 * window, and of its
    first and last `window` bytes beyond that. Catches a file replaced in place, which size and
    mtime alone can miss. For a large file it is not a whole-file digest, so a replacement that
    kept the size and both ends would still slip past; on inputs of this size it costs 0.02 s a
    file, and hashing hundreds of gigabytes on every summarize would not be worth that gap."""
    h = hashlib.sha256()
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        if size <= 2 * window:
            h.update(f.read())
        else:
            h.update(f.read(window))
            f.seek(-window, os.SEEK_END)
            h.update(f.read(window))
    return h.hexdigest()[:16]


def input_fingerprint(row):
    """Identifies the FASTQ files a count was made from, without reading them end to end."""
    out = {}
    for k in ("r1", "r2"):
        if row.get(k):
            st = os.stat(row[k])
            out[k] = {"path": os.path.abspath(row[k]), "bytes": st.st_size,
                      "mtime_ns": st.st_mtime_ns, "digest": file_digest(row[k])}
    return out


def new_tile(counter):
    t = {k: 0 for k in counter.fields}
    for h in getattr(counter, "hists", ("hist",)):
        t[h] = Counter()
    return t


def count_file(path, amp):
    if np is None:
        sys.exit("numpy is missing: activate the junctionseq conda environment")
    counter = SpacingCounter(amp) if amp["mode"] == "spacing" else LandmarkCounter(amp)
    tiles = {}
    for h, s, q in fastq(path):
        k = tile_of(h)
        t = tiles.get(k)
        if t is None:
            t = tiles[k] = new_tile(counter)
        counter.add(t, s, q)
    for t in tiles.values():
        for h in list(t):
            if isinstance(t[h], Counter):
                t[h] = {str(x): n for x, n in sorted(t[h].items())}
    return tiles


# --------------------------------------------------------------------------- run-one

def merge_pairs(row, cfg, out_dir):
    """Merge a read pair with fastp. R1 is given as read 1, so merged reads keep R1's orientation
    and its read name (with the tile). Quality, length and adapter filtering are switched off.
    Poly-G trimming is NOT off: fastp enables it by default on two-colour instruments, where an
    unlit cluster reads as G, and the delivered 3' numbers were produced with it on. It is passed
    explicitly here (`trim_poly_g` in the config) so the behaviour does not depend on fastp
    guessing the platform from the read names."""
    if not shutil.which("fastp"):
        sys.exit("fastp not on PATH (activate the junctionseq conda environment)")
    m = cfg["merge"]
    os.makedirs(f"{out_dir}/merged", exist_ok=True)
    os.makedirs(f"{out_dir}/logs", exist_ok=True)
    final = f"{out_dir}/merged/{row['library']}.merged.fq.gz"
    scratch = os.environ.get("JSEQ_TMP") or (f"/scratch/{os.environ.get('USER', '')}"
                                             if os.path.isdir(f"/scratch/{os.environ.get('USER', '')}")
                                             else None)
    tmpd = tempfile.mkdtemp(prefix=f"jseq_{row['library']}_", dir=scratch)
    try:
        tmp = f"{tmpd}/m.fq.gz"
        cmd = ["fastp", "-i", row["r1"], "-I", row["r2"], "--merge", "--merged_out", tmp,
               "--overlap_len_require", str(m["overlap_len_require"]),
               "--overlap_diff_limit", str(m["overlap_diff_limit"]),
               "--overlap_diff_percent_limit", str(m["overlap_diff_percent_limit"]),
               "--disable_quality_filtering", "--disable_length_filtering",
               "--disable_adapter_trimming", "--thread", str(m["threads"]),
               "--trim_poly_g" if m["trim_poly_g"] else "--disable_trim_poly_g",
               "--poly_g_min_len", str(m["poly_g_min_len"]),
               "--json", f"{out_dir}/logs/{row['library']}.fastp.json", "--html", os.devnull]
        with open(f"{out_dir}/logs/{row['library']}.fastp.err", "w") as err:
            subprocess.run(cmd, check=True, stderr=err)
        shutil.move(tmp, final)
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    return final


def cmd_run_one(args):
    cfg, samples = load_project(args.config)
    if args.index:
        row = samples[args.index - 1]
    else:
        hits = [r for r in samples if r["library"] == args.library]
        if not hits:
            sys.exit(f"library {args.library} not in sample sheet")
        row = hits[0]
    amp = cfg["amplicons"][row["amplicon"]]
    out = args.out
    os.makedirs(f"{out}/counts", exist_ok=True)
    path = row["r1"]
    if row["r2"]:
        path = merge_pairs(row, cfg, out)
    tiles = count_file(path, amp)
    res = {"jseq_version": VERSION, "library": row["library"], "sample": row["sample"],
           "culture": row["culture"], "amplicon": row["amplicon"], "mode": amp["mode"],
           "input": path, "r1": row["r1"], "r2": row["r2"],
           "counting_fingerprint": counting_fingerprint(amp, cfg["merge"], bool(row["r2"])),
           "input_fingerprint": input_fingerprint(row), "tiles": tiles}
    tmp = f"{out}/counts/{row['library']}.json.tmp"
    with open(tmp, "w") as f:
        json.dump(res, f)
    os.replace(tmp, f"{out}/counts/{row['library']}.json")
    if not tiles:
        print(f"{row['library']}: WARNING, no reads in {path}")
        return
    first = next(iter(tiles.values()))
    tot = {k: sum(t[k] for t in tiles.values()) for k in first if not isinstance(first[k], dict)}
    print(f"{row['library']}: {json.dumps(tot)}")


# --------------------------------------------------------------------------- summarize

def denominator(t, mode):
    return t["anchored"] if mode == "landmark" else t["counted"]


def policy(amp):
    """How landmark 2 is used on net-zero calls (landmark mode only):
      veto    (default) drop a net-zero call when landmark 2 shows a rearrangement: found at another
              distance, or contradicted by well-read bases. Reads too poor to check stay net-zero.
      strict  also move net-zero calls that landmark 2 could not check into 'unresolved'.
      off     landmark 1 alone decides (the delivered November 2026 numbers)."""
    return amp.get("landmark2", {}).get("confirm", "veto") if amp["mode"] == "landmark" else None


def net_zero(t, mode, pol=None):
    if mode == "spacing":
        return int(t["hist"].get("0", 0))
    z = t["zero"]
    if pol in ("veto", "off"):
        z += t["zero_l2_lowq"]
    if pol == "off":
        z += t["l2_displaced"] + t["l2_contradicted"]
    return z


def tile_number(key):
    return key.rsplit(":", 1)[-1]


def choose_tiles(amp, libs, qc):
    """Return ({tile_key: kept?}, pooled_rates, plateau) for one amplicon on one flowcell.
    Plateau rule: pool all libraries per tile, take the median net-zero rate of the best
    `plateau_top_n` tiles, keep every tile within `plateau_tolerance_pts` of it."""
    mode = amp["mode"]
    pol = policy(amp)
    pool = defaultdict(lambda: [0, 0])
    for d in libs:
        for k, t in d["tiles"].items():
            pool[k][0] += net_zero(t, mode, pol)
            pool[k][1] += denominator(t, mode)
    rates = {k: 100.0 * z / n for k, (z, n) in pool.items() if n >= qc["min_reads_per_tile_pooled"]}
    tf = amp.get("tile_filter", "none")
    plateau = None
    if tf == "none" or set(pool) == {"NA"}:
        keep = {k: True for k in pool}
    elif tf == "plateau":
        top = sorted(rates.values(), reverse=True)[:qc["plateau_top_n"]]
        plateau = statistics.median(top) if top else float("nan")
        keep = {k: (k in rates and rates[k] >= plateau - qc["plateau_tolerance_pts"]) for k in pool}
    else:
        wanted = {str(x) for x in tf}
        keep = {k: tile_number(k) in wanted for k in pool}
    return keep, rates, plateau


def lib_tile_spread(d, keep, mode, pol, min_reads):
    z = [100.0 * net_zero(t, mode, pol) / denominator(t, mode)
         for k, t in d["tiles"].items() if keep.get(k) and denominator(t, mode) >= min_reads]
    return (max(z) - min(z)) if len(z) > 1 else float("nan"), len(z)


def summarize_library(d, amp, keep, qc):
    mode = amp["mode"]
    pol = policy(amp)
    S = Counter()
    spec = Counter()
    for k, t in d["tiles"].items():
        if not keep.get(k):
            continue
        for f, v in t.items():
            if not isinstance(v, dict):
                S[f] += v
        spec.update({int(x): n for x, n in t["hist"].items()})
        if mode == "landmark":
            spec.update({int(x): n for x, n in t["large_hist"].items()})
            if pol != "off":       # hist[0] holds every landmark-1 net-zero read; take out the vetoed
                spec[0] -= t["l2_displaced"] + t["l2_contradicted"]
                spec.update({int(x): n for x, n in t["displaced_hist"].items()})
            if pol == "strict":
                spec[0] -= t["zero_l2_lowq"]
    if spec.get(0) == 0:
        del spec[0]
    den = denominator(S, mode)
    zero = net_zero(S, mode, pol) if mode == "landmark" else spec.get(0, 0)
    pc = lambda x: 100.0 * x / den if den else float("nan")
    spread_kept, _ = lib_tile_spread(d, keep, mode, pol, qc["min_reads_per_tile_library"])
    spread_all, _ = lib_tile_spread(d, {k: True for k in d["tiles"]}, mode, pol,
                                    qc["min_reads_per_tile_library"])
    row = {"library": d["library"], "sample": d["sample"], "culture": d["culture"],
           "amplicon": d["amplicon"], "mode": mode,
           "reads_total": sum(t["reads"] for t in d["tiles"].values()),
           "reads_on_kept_tiles": S["reads"], "denominator": den,
           "pct_net_zero": pc(zero)}
    if mode == "landmark":
        rearranged = S["l2_displaced"] + S["l2_contradicted"] if pol != "off" else 0
        row.update({
            "confirm_policy": pol,
            "pct_small_indel": pc(S["small"]), "pct_large_deletion": pc(S["large"]),
            "pct_landmark2_displaced": pc(S["l2_displaced"]),
            "pct_landmark2_contradicted": pc(S["l2_contradicted"]),
            "pct_landmark2_low_quality": pc(S["zero_l2_lowq"]),
            "pct_unresolved": pc(S["unresolved"] + (S["zero_l2_lowq"] if pol == "strict" else 0)),
            "pct_net_zero_ceiling": pc(den - S["small"] - S["large"] - rearranged),
            "w1_tie_changed_reads": S["tie_changed"]})
    else:
        row.update({
            "pct_designed_junction_exact": pc(S["designed_exact"]) if amp.get("designed_junction") else "",
            "pct_frame_preserving": pc(sum(v for k, v in spec.items() if k % 3 == 0)),
            "pct_large_indel_gt10bp": pc(sum(v for k, v in spec.items() if abs(k) > 10))})
    row.update({"tiles_kept": sum(1 for k in d["tiles"] if keep.get(k)), "tiles_total": len(d["tiles"]),
                "tile_spread_kept_pts": spread_kept, "tile_spread_all_pts": spread_all,
                "tile_qc": ("n/a" if math.isnan(spread_kept) else
                            "PASS" if spread_kept <= qc["max_tile_spread_pts"] else "FAIL")})
    return row, spec, den


def fmt(v):
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.4f}"
    return v


def write_csv(path, rows, cols=None):
    cols = cols or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "")) for k in cols})


def stale_reasons(d, row, cfg):
    """Why a saved count no longer matches the sample sheet and config it is about to be
    summarized under. Empty list means it is current."""
    why = []
    if d.get("jseq_version") != VERSION:
        why.append(f"counted by jseq {d.get('jseq_version')}, this is {VERSION}")
    for k in ("sample", "culture", "amplicon"):
        if d.get(k) != row[k]:
            why.append(f"{k} was {d.get(k)!r}, sample sheet now says {row[k]!r}")
    if d.get("amplicon") == row["amplicon"]:
        want = counting_fingerprint(cfg["amplicons"][row["amplicon"]], cfg["merge"], bool(row["r2"]))
        if d.get("counting_fingerprint") != want:
            why.append("amplicon or merge settings changed since counting")
    try:
        if d.get("input_fingerprint") != input_fingerprint(row):
            why.append("FASTQ files changed since counting")
    except FileNotFoundError as e:
        why.append(f"FASTQ missing: {e.filename}")
    return why


def cmd_summarize(args):
    cfg, samples = load_project(args.config)
    out = args.out
    res_dir = f"{out}/results"
    os.makedirs(res_dir, exist_ok=True)
    data, missing, stale = {}, [], []
    for r in samples:
        p = f"{out}/counts/{r['library']}.json"
        if not os.path.exists(p):
            missing.append(r["library"])
            continue
        with open(p) as f:
            d = json.load(f)
        why = stale_reasons(d, r, cfg)
        if why and not args.ignore_stale:
            stale.append(f"{r['library']}: {'; '.join(why)}")
            continue
        if why:
            print(f"warning: using stale counts for {r['library']} ({'; '.join(why)})")
        data[r["library"]] = d
    if stale:
        sys.exit("these counts were made from different inputs or settings, so summarizing them "
                 "would report the old run under the new labels. Re-run `run-one` for them, or "
                 "pass --ignore-stale if you are sure.\n  " + "\n  ".join(stale))
    if missing and not args.allow_missing:
        sys.exit(f"no counts for {missing}; run run-one for them first (or pass --allow-missing)")
    qc = cfg["qc"]
    lib_rows, spec_rows, tile_rows, specs = [], [], [], {}
    report = [f"# {cfg['project']}: junction sequencing summary", "",
              f"jseq {VERSION}. Config `{os.path.abspath(args.config)}`.", ""]
    for amp_name, amp in cfg["amplicons"].items():
        libs = [data[r["library"]] for r in samples if r["amplicon"] == amp_name and r["library"] in data]
        if not libs:
            continue
        lanes = sorted({lane_of(k) for d in libs for k in d["tiles"]})
        report += [f"## {amp_name} ({amp['mode']} mode, tile filter: {amp.get('tile_filter', 'none')})", ""]
        keep = {}
        for fc in lanes:
            sub = [{**d, "tiles": {k: t for k, t in d["tiles"].items()
                                   if lane_of(k) == fc}} for d in libs]
            kfc, rates, plateau = choose_tiles(amp, sub, qc)
            keep.update(kfc)
            for k in sorted(kfc):
                tile_rows.append({"amplicon": amp_name, "tile": k, "pooled_pct_net_zero": rates.get(k, float("nan")),
                                  "kept": "yes" if kfc[k] else "no",
                                  "plateau_pct": plateau if plateau is not None else float("nan")})
            nk = sum(kfc.values())
            line = f"- flowcell/lane `{fc}`: kept {nk} of {len(kfc)} tiles"
            if plateau is not None:
                line += (f" (plateau {plateau:.2f}% = median of the best {qc['plateau_top_n']} tiles; "
                         f"kept tiles within {qc['plateau_tolerance_pts']} pt of it)")
            report.append(line)
        report.append("")
        for d in libs:
            row, spec, den = summarize_library(d, amp, keep, qc)
            lib_rows.append(row)
            specs[d["library"]] = (spec, den)
            for k in sorted(spec):
                spec_rows.append({"library": d["library"], "amplicon": amp_name, "net_length_change_bp": k,
                                  "reads": spec[k], "pct_of_denominator": 100.0 * spec[k] / den if den else float("nan")})
        here = [r for r in lib_rows if r["amplicon"] == amp_name]
        na = [r["library"] for r in here if r["tile_qc"] == "n/a"]
        if na:
            report.append(f"Tile QC could not be computed for {', '.join(na)}: fewer than two tiles "
                          f"with at least {qc['min_reads_per_tile_library']} reads. This is not a pass.")
        bad = [r["library"] for r in here if r["tile_qc"] == "FAIL"]
        report.append("Tile QC: " + (("all other libraries PASS" if na else "all libraries PASS")
                      if not bad else
                      f"FAIL for {', '.join(bad)}: net-zero differs by more than "
                      f"{qc['max_tile_spread_pts']} pt between tiles of the same library. "
                      "Do not use these numbers; see README, 'When tile QC fails'."))
        report.append("")
    write_csv(f"{res_dir}/per_library.csv", lib_rows)
    write_csv(f"{res_dir}/spectra.csv", spec_rows)
    write_csv(f"{res_dir}/tiles.csv", tile_rows)
    # per culture: mean of per-library values (each library weighted equally)
    cult_rows = []
    for amp_name in cfg["amplicons"]:
        by = defaultdict(list)
        for r in lib_rows:
            if r["amplicon"] == amp_name:
                by[r["culture"]].append(r)
        cult_means = []
        for c in sorted(by, key=natural):
            v = [r["pct_net_zero"] for r in by[c]]
            m = statistics.mean(v)
            cult_means.append(m)
            cult_rows.append({"amplicon": amp_name, "culture": c, "n_libraries": len(v),
                              "libraries": ";".join(r["library"] for r in by[c]),
                              "mean_pct_net_zero": m,
                              "sd_between_libraries": statistics.stdev(v) if len(v) > 1 else float("nan")})
        if cult_means:
            n = len(cult_means)
            sd = statistics.stdev(cult_means) if n > 1 else float("nan")
            cult_rows.append({"amplicon": amp_name, "culture": "ALL", "n_cultures": n,
                              "mean_pct_net_zero": statistics.mean(cult_means),
                              "sd_between_cultures": sd,
                              "sem_between_cultures": sd / math.sqrt(n) if n > 1 else float("nan")})
    write_csv(f"{res_dir}/per_culture.csv", cult_rows,
              ["amplicon", "culture", "n_libraries", "libraries", "mean_pct_net_zero",
               "sd_between_libraries", "n_cultures", "sd_between_cultures", "sem_between_cultures"])
    report += ["## Net-zero by culture",
               "Culture rows: mean of that culture's libraries, SD between them. ALL row: mean, SD and SEM "
               "of the culture means, so each culture counts once however many libraries it has.", "",
               "| amplicon | culture | libraries | mean % | SD | SEM |", "|---|---|---|---|---|---|"]
    for r in cult_rows:
        all_row = r["culture"] == "ALL"
        report.append(f"| {r['amplicon']} | {r['culture']} | "
                      f"{r['n_cultures'] if all_row else r['n_libraries']}{' cultures' if all_row else ''} | "
                      f"{r['mean_pct_net_zero']:.2f} | "
                      f"{fmt(r['sd_between_cultures'] if all_row else r['sd_between_libraries'])} | "
                      f"{fmt(r['sem_between_cultures']) if all_row else ''} |")
    report.append("")
    if cfg.get("heatmap"):
        try:
            p = heatmap(cfg, samples, specs, res_dir)
            report.append(f"Heatmap: `{p}.png` / `.pdf`; values in `{p}.csv`.")
        except ImportError as e:
            report.append(f"Heatmap skipped: {e}")
    with open(f"{res_dir}/report.md", "w") as f:
        f.write("\n".join(report) + "\n")
    print("\n".join(report))


def natural(s):
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


# --------------------------------------------------------------------------- heatmap

BIN_LABELS = ["<-10"] + [str(b) for b in range(-10, 11)] + [">10"]


def binned(spec, den):
    out = [sum(v for k, v in spec.items() if k < -10)]
    out += [spec.get(b, 0) for b in range(-10, 11)]
    out += [sum(v for k, v in spec.items() if k > 10)]
    return [100.0 * x / den for x in out]


def heatmap(cfg, samples, specs, res_dir):
    """One row per culture and junction, each the mean of that culture's libraries, in the layout
    of the published JPH3 figure (Benedetti et al., Cell 2025, Fig. S1F). In landmark mode a row sums to
    less than 100%: reads that carry the anchor but cannot be sized stay in the denominator."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fd = os.environ.get("JSEQ_FONT_DIR")
    if fd:
        for f in os.listdir(fd):
            if f.lower().endswith(".ttf"):
                font_manager.fontManager.addfont(os.path.join(fd, f))
    names = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in names else "DejaVu Sans", "font.size": 12})
    hm = cfg["heatmap"]
    cultures = sorted({r["culture"] for r in samples if r["library"] in specs}, key=natural)
    M, labels, rows_out = [], [], []
    for rep, c in enumerate(cultures, 1):
        for jr in hm["rows"]:
            libs = [r["library"] for r in samples
                    if r["culture"] == c and r["amplicon"] == jr["amplicon"] and r["library"] in specs]
            if not libs:
                continue
            m = np.mean([binned(*specs[l]) for l in libs], axis=0)
            M.append(m)
            labels.append((rep, jr["label"], m[11]))
            rows_out.append([f"Replicate {rep}", c, jr["label"], ";".join(libs)] + [f"{x:.3f}" for x in m])
    M = np.array(M)
    nrow = len(M)
    fig = plt.figure(figsize=(12.4, 1.15 + 0.32 * nrow))
    L, W = 0.215, 0.695
    H = 0.32 * nrow / (1.15 + 0.32 * nrow)
    B = 0.62 / (1.15 + 0.32 * nrow)
    ax = fig.add_axes([L, B, W, H])
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=100, interpolation="nearest")
    ax.set_yticks([])
    ax.set_xticks([0, 1, 6, 11, 16, 21, 22])
    ax.set_xticklabels(["<-10", "-10", "-5", "0", "5", "10", ">10"])
    ax.set_xticks(range(23), minor=True)
    ax.set_xlabel("Indel (bp)", labelpad=2)
    ax.set_title("Indel frequency heatmap", fontsize=12, pad=6)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = fig.add_axes([0.922, B, 0.011, H])
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([0, 20, 40, 60, 80, 100])
    cb.outline.set_visible(False)
    fig.canvas.draw()
    R = fig.canvas.get_renderer()

    def width(t):
        tt = fig.text(0, 0, t)
        w = tt.get_window_extent(R).width / fig.bbox.width
        tt.remove()
        return w
    green = "#3E9E4F"
    rowh = H / nrow
    xg = L - 0.008
    xb = xg - width("(00.0%)") - 0.004
    xbr = xb - width("5′") - 0.008
    for i, (rep, jn, pct) in enumerate(labels):
        y = B + H - (i + 0.5) * rowh
        fig.text(xg, y, f"({pct:.1f}%)", ha="right", va="center", color=green)
        fig.text(xb, y, jn, ha="right", va="center")
    for rep in range(1, len(cultures) + 1):
        idx = [i for i, l in enumerate(labels) if l[0] == rep]
        if not idx:
            continue
        y0 = B + H - idx[0] * rowh - 0.012
        y1 = B + H - (idx[-1] + 1) * rowh + 0.012
        fig.add_artist(plt.Line2D([xbr, xbr], [y1, y0], color="black", lw=1.0, transform=fig.transFigure))
        fig.text(xbr - 0.006, (y0 + y1) / 2, f"Replicate {rep}", ha="right", va="center")
    fig.text(xg, B + H + 0.04, "(% no indels)", ha="right", va="bottom", color=green)
    fig.text(xg - width("(% no indels)") - 0.004, B + H + 0.04, hm["label"], ha="right", va="bottom")
    stem = f"{res_dir}/heatmap_by_culture"
    fig.savefig(stem + ".png", dpi=300, facecolor="white")
    fig.savefig(stem + ".pdf", facecolor="white")
    plt.close(fig)
    with open(stem + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["replicate", "culture", "junction", "libraries_averaged"] + BIN_LABELS)
        w.writerows(rows_out)
    return stem


# --------------------------------------------------------------------------- small commands

def cmd_list(args):
    _, samples = load_project(args.config)
    for r in samples:
        print(r["library"])


def cmd_check(args):
    cfg, samples = load_project(args.config)
    bad = [p for r in samples for p in (r["r1"], r["r2"]) if p and not os.path.exists(p)]
    for name, a in cfg["amplicons"].items():
        ref = a["reference"]
        i = ref.find(a["anchor"])
        if a["mode"] == "spacing":
            j = ref.find(a["end_anchor"]) + len(a["end_anchor"])
            print(f"{name}: spacing mode, anchors at {i} and {j - len(a['end_anchor'])}, reference spacing {j - i}")
        else:
            print(f"{name}: landmark mode, anchor at {i}, landmark 1 at {ref.find(a['landmark1']['seq'])}, "
                  f"landmark 2 at {ref.find(a['landmark2']['seq'])} (reference length {len(ref)})")
    print(f"{len(samples)} libraries, {len(bad)} missing FASTQ files")
    for p in bad:
        print("  missing:", p)
    sys.exit(1 if bad else 0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=VERSION)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-one", help="merge (if paired) and count one library")
    p.add_argument("config")
    p.add_argument("out")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--index", type=int, help="1-based row in the sample sheet (for LSF array jobs)")
    g.add_argument("--library")
    p.set_defaults(func=cmd_run_one)
    p = sub.add_parser("summarize", help="tile QC, tables and heatmap from the counts")
    p.add_argument("config")
    p.add_argument("out")
    p.add_argument("--allow-missing", action="store_true",
                   help="summarize the libraries that have counts, ignoring the rest")
    p.add_argument("--ignore-stale", action="store_true",
                   help="summarize counts even if the inputs or counting settings have changed")
    p.set_defaults(func=cmd_summarize)
    p = sub.add_parser("list-samples")
    p.add_argument("config")
    p.set_defaults(func=cmd_list)
    p = sub.add_parser("check")
    p.add_argument("config")
    p.set_defaults(func=cmd_check)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
