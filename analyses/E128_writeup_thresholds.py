"""Write-up check for E128 part B (no simulation): the passing fraction under other threshold rules for the responders' CV.
Deals re-drawn with the script's own seeds, so rule 'grouped' reproduces b_wirings.csv.

Prints as before and, unless --no-save, also saves its tables to results/E128_bouton_groups/ (`E.save_writeup`):
writeup_thresholds_rates.csv (each fly's response rate among its grouped ROIs and among all its ROIs),
writeup_thresholds.csv (per wiring, the passing fraction and the deals' median responders' CV under each rule) and
writeup_thresholds_summary.csv (per rule, the wirings with at least 5 % of deals passing and the range of the passing
fraction over the wirings, in percent), with provenance_writeup_thresholds.json.

Usage (from the repository root): ./venv/bin/python analyses/E128_writeup_thresholds.py [--no-save] [--out DIR]
"""
import argparse
from pathlib import Path

import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E
from kcq import harness as H

RULES = ('grouped', 'all_roi_rate', 'e119_fixed')


def summary_table(T: pd.DataFrame) -> pd.DataFrame:
    """Per rule: the number of wirings (of T's) with at least 5 % of deals passing, the lowest and highest passing
    fraction over the wirings in percent with the wirings where they are, and the range of the deals' median
    responders' CV."""
    return pd.DataFrame([dict(rule=c, n_wirings_ge_5pct=int((T[c] >= 0.05).sum()), n_wirings=len(T),
                              pct_pass_min=100 * T[c].min(), wiring_min=T[c].idxmin(),
                              pct_pass_max=100 * T[c].max(), wiring_max=T[c].idxmax(),
                              resp_cv_median_lo=T[f'rcv_{c}'].min(), resp_cv_median_hi=T[f'rcv_{c}'].max())
                         for c in RULES])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--no-save', action='store_true', help='print only, as before; save no table')
    ap.add_argument('--out', type=Path, default=None, help='folder for the tables (default: the results folder)')
    args = ap.parse_args(argv)

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

    if args.no_save:
        return
    rates = pd.DataFrame([dict(acq=acq, rate_grouped=f['rate_g'], rate_all_rois=f['rate_all'], n_groups=len(f['sizes']))
                          for acq, f in sorted(info.items())])
    E.save_writeup(dict(thresholds_rates=rates, thresholds=T.reset_index(), thresholds_summary=summary_table(T)),
                   H.script_path(__file__), args.out, mean_all_roi_rate=mean_all, n_deals=E.N_DEALS, seed=E.SEED,
                   flies_range=sp, rules=dict(
                       grouped="the threshold matched to the fly's response rate among its grouped ROIs (part B's rule)",
                       all_roi_rate="the threshold matched to the fly's response rate among all its ROIs",
                       e119_fixed="one threshold per wiring, matched on every bouton once to the flies' mean rate "
                                  "among all ROIs (E119's rule)"))


if __name__ == '__main__':
    main()
