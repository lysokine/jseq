"""Effect of the landmark-2 confirm policy on the 5' numbers, with the tile set held fixed.
Separates the policy's effect from the tile rule's, which an earlier comparison conflated."""
import json, os, sys
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
keep, _, _ = jseq.choose_tiles(amp, list(data.values()), cfg['qc'])

hdr = ('library', 'veto', 'strict', 'str-veto', 'off', 'off-veto', 'ceiling', 'ceil-veto')
print(f'{hdr[0]:<12s}' + ''.join(f'{h:>10s}' for h in hdr[1:]))
rows, ds, do, dc = [], [], [], []
for l in libs:
    r = {}
    for pol in ('veto', 'strict', 'off'):
        a = dict(amp, landmark2=dict(amp['landmark2'], confirm=pol))
        r[pol] = jseq.summarize_library(data[l], a, keep, cfg['qc'])[0]
    v, s, o = (r[k]['pct_net_zero'] for k in ('veto', 'strict', 'off'))
    c = r['veto']['pct_net_zero_ceiling']
    ds.append(s - v); do.append(o - v); dc.append(c - v)
    rows.append({'library': l, 'veto': v, 'strict': s, 'off': o, 'ceiling': c})
    print(f'{l:<12s}{v:>9.2f}%{s:>9.2f}%{s-v:>+10.2f}{o:>9.2f}%{o-v:>+10.2f}{c:>9.2f}%{c-v:>+10.2f}')
rng = lambda x: (min(x), max(x))
print(f'\nstrict - veto : {rng(ds)[0]:+.2f} to {rng(ds)[1]:+.2f} points')
print(f'off    - veto : {rng(do)[0]:+.2f} to {rng(do)[1]:+.2f} points')
print(f'ceiling- veto : {rng(dc)[0]:.2f} to {rng(dc)[1]:.2f} points')
# How much could `veto` be wrong by? Among net-zero calls whose landmark 2 IS readable, a known
# fraction turn out rearranged. If the unreadable-end reads carry rearrangements at that same rate,
# the error introduced by counting them as net-zero is that rate times the size of the class.
print(f'\n{"library":<12s}{"rearranged rate":>16s}{"unreadable %":>14s}{"implied error":>15s}')
err = []
for l in libs:
    a = dict(amp, landmark2=dict(amp['landmark2'], confirm='veto'))
    r = jseq.summarize_library(data[l], a, keep, cfg['qc'])[0]
    rear = r['pct_landmark2_displaced'] + r['pct_landmark2_contradicted']
    lowq = r['pct_landmark2_low_quality']
    readable_zero = r['pct_net_zero'] - lowq
    rate = rear / (readable_zero + rear)
    err.append(100 * rate * lowq / 100)
    print(f'{l:<12s}{100*rate:>15.3f}%{lowq:>13.3f}%{100*rate*lowq/100:>14.4f}')
print(f'implied error from counting unreadable-end reads as net-zero: '
      f'{min(err):.4f} to {max(err):.4f} points')

json.dump({'tiles_kept': sum(keep.values()), 'per_library': rows,
           'implied_error_pts': [min(err), max(err)],
           'strict_minus_veto_pts': rng(ds), 'off_minus_veto_pts': rng(do),
           'ceiling_minus_veto_pts': rng(dc)},
          open(f'{P}/validation/policy_effect.json', 'w'), indent=1)
print('wrote results/policy_effect.json')
