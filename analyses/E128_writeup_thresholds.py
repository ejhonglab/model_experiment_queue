"""Write-up check for E128 part B (no simulation): the passing fraction under other threshold rules for the responders' CV.
Deals re-drawn with the script's own seeds, so rule 'grouped' reproduces b_wirings.csv."""
import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E
E119 = E._e119()
flies = E.fly_means()
sp = E.flies_range()
info = {}
for acq, X in sorted(flies.items()):
    M = X.mean(axis=0); lab = E.fly_groups(M)
    info[acq] = dict(sizes=np.unique(lab[lab > 0], return_counts=True)[1], rate_g=E.fly_rate(M, lab),
                     rate_all=float((M > E.CRITERION).mean()))
    print(acq[:21], 'rate grouped', round(info[acq]['rate_g'], 4), 'rate all ROIs', round(info[acq]['rate_all'], 4))
mean_all = float(np.mean([f['rate_all'] for f in info.values()]))
rows = []
for wi, w in enumerate(E119.WIRINGS):
    U, gl = E.model_boutons(w); gloms = np.unique(gl)
    thr119 = E119.matched_threshold(pd.DataFrame(U, columns=list(E119.ODORS)), mean_all)
    for fi, (acq, f) in enumerate(sorted(info.items())):
        rng = np.random.default_rng([E.SEED, wi, fi, 0])
        for d in range(E.N_DEALS):
            P, _ = E.deal_population(U, gl, f['sizes'], rng, gloms)
            r = dict(wiring=w, acq=acq)
            for name, thr in (('grouped', float(np.quantile(P.ravel(), 1 - f['rate_g']))),
                              ('all_roi_rate', float(np.quantile(P.ravel(), 1 - f['rate_all']))),
                              ('e119_fixed', thr119)):
                st = E.three_stats(P, thr)
                r[name] = E.model_pass(st, sp)
                r[f'rcv_{name}'] = st['resp_cv']
            rows.append(r)
D = pd.DataFrame(rows)
T = D.groupby('wiring', sort=False)[['grouped', 'all_roi_rate', 'e119_fixed', 'rcv_grouped', 'rcv_all_roi_rate', 'rcv_e119_fixed']].agg(
    {'grouped': 'mean', 'all_roi_rate': 'mean', 'e119_fixed': 'mean', 'rcv_grouped': 'median', 'rcv_all_roi_rate': 'median', 'rcv_e119_fixed': 'median'})
print(T.round(4).to_string())
print('wirings >= 5%:', {c: int((T[c] >= 0.05).sum()) for c in ('grouped', 'all_roi_rate', 'e119_fixed')})
print('mean all-ROI rate', round(mean_all, 4))
