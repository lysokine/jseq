"""Does choosing tiles on the pooled net-zero rate leak into the per-library numbers?

The plateau rule selects tiles by the quantity it then reports, so in principle it could inflate
that quantity by selecting on noise. Two checks:

  leave-one-out   choose the tiles from the other 11 libraries only, then measure this library.
                  If the selection carried information about this library, the held-out value
                  would be lower than the in-sample one.
  noise floor     binomial standard error of the pooled per-tile rate, to compare against the
                  1 point plateau tolerance.
"""
import csv, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import jseq

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the package root
CFG = sys.argv[1] if len(sys.argv) > 1 else f'{P}/config/study.json'
if not os.path.exists(CFG):
    sys.exit(f"usage: {sys.argv[0]} <project config.json>\n"
             "The study project's config and counts are not distributed with the code, because that "
             "study is not published yet. Point this at a project config with landmark-mode results.")
cfg, samples = jseq.load_project(CFG)
AMP = next(k for k, v in cfg['amplicons'].items() if v['mode'] == 'landmark')
amp = cfg['amplicons'][AMP]
libs = [r['library'] for r in samples if r['amplicon'] == AMP]
OUT = sys.argv[2] if len(sys.argv) > 2 else f'{P}/example_output/{cfg["project"].lower()}'
data = {l: json.load(open(f'{OUT}/counts/{l}.json')) for l in libs}

keep_all, rates, plateau = jseq.choose_tiles(amp, list(data.values()), cfg['qc'])
n = [sum(t['anchored'] for d in data.values() for k, t in d['tiles'].items() if jseq.tile_number(k) == jseq.tile_number(kk))
     for kk in [next(iter(keep_all))]]
print(f'in-sample plateau {plateau:.3f}%, {sum(keep_all.values())} of {len(keep_all)} tiles kept')
pooled_n = {k: sum(d['tiles'][k]['anchored'] for d in data.values() if k in d['tiles']) for k in keep_all}
se = {k: 100 * math.sqrt(max(rates.get(k, 0), 1e-9) / 100 * (1 - rates.get(k, 0) / 100) / pooled_n[k])
      for k in keep_all if pooled_n[k]}
print(f'pooled reads per tile {min(pooled_n.values()):,}-{max(pooled_n.values()):,}; '
      f'binomial SE of a tile rate {min(se.values()):.4f}-{max(se.values()):.4f} points '
      f'(plateau tolerance {cfg["qc"]["plateau_tolerance_pts"]} points)')

print(f'\n{"library":<12s} {"in-sample":>10s} {"held-out":>9s} {"diff":>6s}  tiles(in/held)')
worst = 0.0
rows = []
for l in libs:
    others = [d for k, d in data.items() if k != l]
    keep_lo, _, pl_lo = jseq.choose_tiles(amp, others, cfg['qc'])
    a, _, _ = jseq.summarize_library(data[l], amp, keep_all, cfg['qc'])
    b, _, _ = jseq.summarize_library(data[l], amp, keep_lo, cfg['qc'])
    d = b['pct_net_zero'] - a['pct_net_zero']
    worst = max(worst, abs(d))
    rows.append((l, a['pct_net_zero'], b['pct_net_zero'], d))
    print(f'{l:<12s} {a["pct_net_zero"]:>9.3f}% {b["pct_net_zero"]:>8.3f}% {d:>+6.3f}  '
          f'{sum(keep_all.values())}/{sum(keep_lo.values())}')
print(f'\nlargest leave-one-out change: {worst:.3f} points')
json.dump({'plateau_pct': plateau, 'tiles_kept': sum(keep_all.values()), 'tiles_total': len(keep_all),
           'pooled_reads_per_tile': [min(pooled_n.values()), max(pooled_n.values())],
           'binomial_se_pts': [min(se.values()), max(se.values())],
           'leave_one_out': [{'library': r[0], 'in_sample': r[1], 'held_out': r[2], 'diff': r[3]} for r in rows],
           'max_abs_diff_pts': worst},
          open(f'{P}/validation/tile_rule_leakage.json', 'w'), indent=1)
print('wrote results/tile_rule_leakage.json')
