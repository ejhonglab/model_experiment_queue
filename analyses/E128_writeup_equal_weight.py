"""Write-up checks for E128 (no simulation): equal weight with groups below k pooled, and the amplitude failures in part B."""
import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E

flies = E.fly_means()
rows = []
for acq, X in sorted(flies.items()):
    M = X.mean(axis=0)
    lab = E.fly_groups(M)
    ids, S, n = E._group_sums(M, lab)
    g = ids != 0
    r = dict(acq=acq[:21], measured=E.cv(S.sum(axis=0)), equal=E.cv((S / n[:, None]).mean(axis=0)))
    for kmin in (2, 3, 5, 10, 20):
        big = g & (n >= kmin)
        parts = list(S[big] / n[big, None])
        if (~big).any():
            parts.append(S[~big].sum(axis=0) / n[~big].sum())
        r[f'equal_ge{kmin}'] = E.cv(np.mean(parts, axis=0))
    # equal weight by square root of size (between measured and equal)
    w = np.sqrt(n)
    r['sqrt_weight'] = E.cv(((S / n[:, None]) * w[:, None]).sum(axis=0) / w.sum())
    # equal weight among groups >=5, drop small groups and remainder entirely
    b5 = g & (n >= 5)
    r['equal_ge5_only'] = E.cv((S[b5] / n[b5, None]).mean(axis=0))
    rows.append(r)
T = pd.DataFrame(rows)
pd.set_option('display.width', 200)
print(T.round(3).to_string(index=False))
print('flies below 0.230 per reading:', {c: int((T[c] < 0.230).sum()) for c in T.columns if c != 'acq'})

D = pd.read_parquet('outputs/E128_bouton_groups/b_deals.parquet')
sp = E.flies_range()
lo, hi = sp['totals_cv']
D['tin'] = D.totals_cv.between(lo, hi)
D['resp_above'] = D.resp_cv > sp['resp_cv'][1]
D['top_above'] = D.top1 > sp['top1'][1]
for v in ('all', 'big'):
    sub = D[(D.variant == v) & D.tin]
    print(v, 'in-range deals', len(sub), 'resp above', round(sub.resp_above.mean(), 3), 'top1 above', round(sub.top_above.mean(), 3),
          'both', round((sub.resp_above & sub.top_above).mean(), 3), 'either', round((sub.resp_above | sub.top_above).mean(), 3))
    allv = D[D.variant == v]
    print('  overall resp above', round(allv.resp_above.mean(), 3), 'top1 above', round(allv.top_above.mean(), 3),
          ' totals >= lo', round((allv.totals_cv >= lo).mean(), 3), ' totals > hi', round((allv.totals_cv > hi).mean(), 3))
    # per wiring: fraction totals in range alone
    print('  totals-only fraction per wiring:', allv.groupby('wiring', sort=False).tin.mean().round(4).to_dict())
    print('  wirings >= 5% totals-only:', int((allv.groupby('wiring').tin.mean() >= 0.05).sum()))
    # per wiring passing if only resp or only top1 condition
    pr = allv.assign(p_resp=lambda d: d.tin & ~d.resp_above, p_top=lambda d: d.tin & ~d.top_above).groupby('wiring', sort=False)[['p_resp', 'p_top', 'passes']].mean()
    print(pr.round(4).to_string())
# median n_rois of deals, n_gloms
print(D.groupby('variant').n_gloms.median(), D.groupby('variant').n_rois.median())
# correlation in deals: top-weight share vs totals_cv
