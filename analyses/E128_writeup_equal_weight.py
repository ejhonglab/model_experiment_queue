"""Write-up checks for E128 (no simulation): equal weight with groups below k pooled, and the amplitude failures in part B.

Prints as before and, unless --no-save, also saves its tables to results/E128_bouton_groups/ (`E.save_writeup`):
writeup_equal_weight.csv (per fly, the readings printed first), writeup_equal_weight_below.csv (flies below 0.230 per
reading), writeup_group_sizes.csv and writeup_group_sizes_range.csv (how a fly's grouped ROIs are spread over its
groups, per fly and as the range over the five flies), writeup_amplitude.csv and writeup_amplitude_summary.csv (part
B's deals inside the flies' range of totals: which amplitude condition fails, per wiring and per variant), with
provenance_writeup_equal_weight.json. Part B's deals are read from outputs/E128_bouton_groups/b_deals.parquet where the
experiment was run, and are otherwise drawn again with the script's seeds (`E.part_b`), checked against b_wirings.csv.

Usage (from the repository root): ./venv/bin/python analyses/E128_writeup_equal_weight.py [--no-save] [--out DIR]
"""
import argparse
import contextlib
import io
from pathlib import Path

import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E
from kcq import harness as H

KMINS = (2, 3, 5, 10, 20)
DEALS = Path(__file__).resolve().parents[1] / 'outputs' / E.EXPERIMENT / 'b_deals.parquet'


def equal_weight_table(flies) -> pd.DataFrame:
    """Per fly (full acquisition names): the totals' CV measured; at equal weight; at equal weight with the groups
    under k ROIs pooled with the remainder into one group (equal_ge<k>); weighted by the square root of size; and at
    equal weight among the groups of >= 5 ROIs only (the small groups and the remainder dropped)."""
    rows = []
    for acq, X in sorted(flies.items()):
        M = X.mean(axis=0)
        lab = E.fly_groups(M)
        ids, S, n = E._group_sums(M, lab)
        g = ids != 0
        r = dict(acq=acq, measured=E.cv(S.sum(axis=0)), equal=E.cv((S / n[:, None]).mean(axis=0)))
        for kmin in KMINS:
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
    return pd.DataFrame(rows)


def below_table(T: pd.DataFrame) -> pd.DataFrame:
    """Per reading of `equal_weight_table`: the number of flies below the half-way point (0.230)."""
    return pd.DataFrame([dict(reading=c, n_flies_below=int((T[c] < E.HALF_WAY).sum()), n_flies=len(T),
                              half_way=E.HALF_WAY) for c in T.columns if c != 'acq'])


def group_sizes_table(flies) -> pd.DataFrame:
    """Per fly, how its grouped ROIs (the remainder out) are spread over its tuning groups, in percent of the grouped
    ROIs: in the three groups of largest variance share (the ones left out), in the three largest groups, and in the
    groups of fewer than five ROIs; and those small groups' share of the weight at equal weight (their number over the
    number of groups plus the remainder). n_effective is the effective number of groups by size, 1 / sum(p^2)."""
    rows = []
    for acq, X in sorted(flies.items()):
        M = X.mean(axis=0)
        lab = E.fly_groups(M)
        ids, n = np.unique(lab[lab > 0], return_counts=True)
        size = pd.Series(n, index=ids)
        N = int(size.sum())
        top = E.variance_shares(M, lab).drop(0, errors='ignore').sort_values(ascending=False).index[:3]
        small = size[size < E.MIN_GROUP]
        has_rem = int((lab == 0).any())
        rows.append(dict(acq=acq, n_grouped=N, n_groups=len(size), n_groups_small=len(small), remainder=bool(has_rem),
                         pct_grouped_top3_by_share=100 * size[top].sum() / N,
                         pct_grouped_top3_by_size=100 * size.nlargest(3).sum() / N,
                         pct_grouped_small=100 * small.sum() / N,
                         pct_weight_small=100 * len(small) / (len(size) + has_rem),
                         n_effective=float(1 / ((size / N) ** 2).sum())))
    return pd.DataFrame(rows)


def ranges(T: pd.DataFrame, by: str) -> pd.DataFrame:
    """Per numeric column of T: its lowest and highest value, and the rows (`by`) where they are."""
    return pd.DataFrame([dict(statistic=c, lo=T[c].min(), hi=T[c].max(), **{f'{by}_lo': T[by][T[c].idxmin()],
                                                                           f'{by}_hi': T[by][T[c].idxmax()]})
                         for c in T.columns if c != by and pd.api.types.is_numeric_dtype(T[c])
                         and not pd.api.types.is_bool_dtype(T[c])])


def deals() -> tuple:
    """Part B's deals: (deals, source), from the local parquet if there, else drawn again with the script's seeds;
    a redraw must reproduce b_wirings.csv's passing and in-range fractions."""
    if DEALS.exists():
        return pd.read_parquet(DEALS), str(DEALS.relative_to(DEALS.parents[2]))
    with contextlib.redirect_stdout(io.StringIO()):
        D = E.part_b()['deals']
    W = E.wiring_table(D).set_index(['wiring', 'variant'])
    S = pd.read_csv(E.out_dir() / 'b_wirings.csv').set_index(['wiring', 'variant']).loc[W.index]
    for c in ('frac_pass', 'frac_totals_in'):
        if not np.allclose(W[c], S[c], rtol=0, atol=1e-12):
            raise RuntimeError(f'the redrawn deals do not reproduce b_wirings.csv ({c})')
    return D, "drawn again with E128_bouton_groups.part_b's seeds (b_wirings.csv reproduced)"


def amplitude_tables(D: pd.DataFrame, sp: dict):
    """Part B's deals with the totals' CV inside the flies' range: per variant and wiring, the fractions of deals in
    range and passing with each amplitude condition alone (writeup_amplitude.csv); per variant, among the in-range deals
    the fractions above the flies' on the responders' CV, the top-1 % share, both or either, the same over all deals,
    the number of wirings with at least 5 % of deals in range, and the deals' median sizes (writeup_amplitude_summary)."""
    lo, hi = sp['totals_cv']
    rows, summ = [], []
    for v in ('all', 'big'):
        sub = D[(D.variant == v) & D.tin]
        allv = D[D.variant == v]
        pr = allv.assign(p_resp=lambda d: d.tin & ~d.resp_above, p_top=lambda d: d.tin & ~d.top_above).groupby(
            'wiring', sort=False)[['tin', 'p_resp', 'p_top', 'passes']].mean()
        rows.append(pr.rename(columns={'tin': 'frac_totals_in'}).reset_index().assign(variant=v))
        summ.append(dict(variant=v, n_deals=len(allv), n_deals_totals_in=len(sub),
                         frac_resp_above=sub.resp_above.mean(), frac_top1_above=sub.top_above.mean(),
                         frac_both_above=(sub.resp_above & sub.top_above).mean(),
                         frac_either_above=(sub.resp_above | sub.top_above).mean(),
                         overall_resp_above=allv.resp_above.mean(), overall_top1_above=allv.top_above.mean(),
                         frac_totals_ge_lo=(allv.totals_cv >= lo).mean(), frac_totals_gt_hi=(allv.totals_cv > hi).mean(),
                         n_wirings_totals_in_ge_5pct=int((allv.groupby('wiring').tin.mean() >= 0.05).sum()),
                         n_wirings=int(allv.wiring.nunique()),
                         n_gloms_median=allv.n_gloms.median(), n_rois_median=allv.n_rois.median()))
    A = pd.concat(rows, ignore_index=True)
    return A[['variant'] + [c for c in A.columns if c != 'variant']], pd.DataFrame(summ)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--no-save', action='store_true', help='print only, as before; save no table')
    ap.add_argument('--out', type=Path, default=None, help='folder for the tables (default: the results folder)')
    args = ap.parse_args(argv)

    flies = E.fly_means()
    T = equal_weight_table(flies)
    pd.set_option('display.width', 200)
    print(T.assign(acq=T.acq.str[:21]).round(3).to_string(index=False))
    print('flies below 0.230 per reading:', {c: int((T[c] < 0.230).sum()) for c in T.columns if c != 'acq'})

    D, source = deals()
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

    if args.no_save:
        return
    G = group_sizes_table(flies)
    A, S = amplitude_tables(D, sp)
    E.save_writeup(dict(equal_weight=T, equal_weight_below=below_table(T), group_sizes=G,
                        group_sizes_range=ranges(G, 'acq'), amplitude=A, amplitude_summary=S),
                   H.script_path(__file__), args.out, deals=source, flies_range=sp, half_way=E.HALF_WAY,
                   min_group=E.MIN_GROUP)


if __name__ == '__main__':
    main()
