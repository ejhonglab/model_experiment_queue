"""Write-up check for E128: the model's per-glomerulus bouton profile (real wiring, APL on) of the lead glomeruli.

Prints as before and, unless --no-save, also saves its tables to results/E128_bouton_groups/ (`E.save_writeup`):
writeup_leads.csv (per glomerulus, unrounded: boutons, mean response, CV over odors, top odor and its share),
writeup_leads_profiles.csv (every glomerulus's mean bouton response per odor) and writeup_model_sample.csv (how the
wiring's boutons are spread over its glomeruli: the share on the three largest, in percent, and the effective number
of glomeruli by size), with provenance_writeup_leads.json.

Usage (from the repository root): ./venv/bin/python analyses/E128_writeup_leads.py [--no-save] [--out DIR]
"""
import argparse
from pathlib import Path

import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E
from kcq import harness as H

WIRING = 'real'


def sample_table(nb: pd.Series, wiring: str = WIRING) -> pd.DataFrame:
    """How a wiring's boutons (`nb`, the number per glomerulus) are spread over its glomeruli."""
    p = nb / nb.sum()
    return pd.DataFrame([dict(wiring=wiring, n_glomeruli=len(nb), n_boutons=int(nb.sum()),
                              n_boutons_top3=int(nb.nlargest(3).sum()), pct_boutons_top3=100 * nb.nlargest(3).sum() / nb.sum(),
                              n_effective=float(1 / (p ** 2).sum()))])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--no-save', action='store_true', help='print only, as before; save no table')
    ap.add_argument('--out', type=Path, default=None, help='folder for the tables (default: the results folder)')
    args = ap.parse_args(argv)

    E119 = E._e119()
    U, gl = E.model_boutons(WIRING)
    df = pd.DataFrame(U, columns=list(E119.ODORS)); df['g'] = gl
    m = df.groupby('g').mean()
    tot = m.mean(axis=1)
    cvg = m.std(axis=1, ddof=0) / m.mean(axis=1)
    nb = df.groupby('g').size()
    top_share = m.max(axis=1) / m.clip(lower=0).sum(axis=1)
    out = pd.DataFrame({'n_boutons': nb, 'mean': tot.round(2), 'cv_over_odors': cvg.round(2),
                        'top_odor': m.idxmax(axis=1), 'top_share': top_share.round(2)})
    print(out.sort_values('cv_over_odors', ascending=False).to_string())
    for g in ('DM5', 'DM3', 'DM4'):
        print(g, m.loc[g].round(1).sort_values(ascending=False).to_dict())
    print('odors', list(E119.ODORS))

    if args.no_save:
        return
    raw = pd.DataFrame({'n_boutons': nb, 'mean': tot, 'cv_over_odors': cvg, 'top_odor': m.idxmax(axis=1),
                        'top_share': top_share})
    raw = raw.loc[out.sort_values('cv_over_odors', ascending=False).index].rename_axis('glomerulus').reset_index()
    E.save_writeup(dict(leads=raw.assign(wiring=WIRING), leads_profiles=m.rename_axis('glomerulus').reset_index().assign(wiring=WIRING),
                        model_sample=sample_table(nb)),
                   H.script_path(__file__), args.out, wiring=WIRING, cond='on', reading='rate',
                   e90_archive=str(E119.archive_dir()))


if __name__ == '__main__':
    main()
