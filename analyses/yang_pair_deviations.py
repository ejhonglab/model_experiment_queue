#!/usr/bin/env python3
"""Which odor pairs deviate most, and most consistently, from the monotone (isotonic) fit of the control flies' KC
matrix on the matrix with APL output blocked (TNT); whether all nine fruity x small-ester (block) pairs do, and which
other pairs do. Yang's question on her version of the summary page (2026-11-16); this script is the source of the
numbers of that answer.

All values are cosine similarities (1 - the cosine distance of the main page), per recording from
`kcq.harness.yang_recording_rdms`, averaged over the recordings of a genotype (11 control, 8 TNT). The residual of a
pair is the control flies' similarity minus the isotonic fit of control on TNT at the pair's TNT similarity, so a
negative residual is a pair that APL decorrelates more than the trend (the block's pairs). The fit's residual is the
negative of `kcq.harness.isotonic_resid`'s on the distances (checked).

1. Per pair: the residual and the change of rank (rank among the 120 pairs, 1 = least similar; control minus TNT),
   each with a 90 % interval, a two-sided p and a Benjamini-Hochberg q over the 120 pairs from 4000 resamples of the
   recordings with replacement (each genotype separately, seed 1); the share of resamples with the residual's sign;
   and how many single recordings sit on the side of the residual (a control recording against the TNT mean, a TNT
   recording against the control mean).
2. Per odor, from all 120 pairs: the mean residual of the odor's 15 pairs, and the odor effects of an additive model
   (a pair's residual = the sum of its two odors' effects, least squares), with 90 % intervals from 3000 resamples
   (seed 2), and the share of the residuals' variance the odor effects explain.
3. The same odor effects fitted on the 111 pairs outside the block: what they predict for the nine block pairs, and
   the remainder, each block pair's residual beyond its odors' effects, with 90 % intervals from 2500 resamples
   (seed 3).

Writes <shared>/results/yang_pair_deviations/: pairs.csv, odors.csv, block_pairs.csv, summary.csv (statistic, value,
lo, hi) and provenance.json.
Usage: ./venv/bin/python analyses/yang_pair_deviations.py
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from typing import Callable, Dict, Iterator, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.isotonic import IsotonicRegression

from kcq import harness as H

SEEDS = dict(pairs=1, odor_effects=2, outside_block=3)
RESAMPLES = dict(pairs=4000, odor_effects=3000, outside_block=2500)
GENOTYPES = ('TNT3c', 'TNT')   # control, APL output blocked (TNT)


def iso_resid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """`y` minus its increasing isotonic fit on `x` (y: control similarities, x: TNT similarities)."""
    return y - IsotonicRegression(increasing=True, out_of_bounds='clip').fit(x, y).predict(x)


def rank_change(c: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Each pair's rank among the pairs in control minus its rank with APL output blocked (1 = least similar)."""
    return rankdata(c) - rankdata(t)


def resamples(C: np.ndarray, T: np.ndarray, n: int, seed: int) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
    """`n` resampled (control mean, TNT mean) pairs: the recordings of each genotype drawn with replacement, the
    control draw first, from `np.random.default_rng(seed)`.

    Args:
        C, T: recordings x pairs similarities of the control and the TNT flies.
    """
    rng = np.random.default_rng(seed)
    for _ in range(n):
        c = C[rng.integers(0, len(C), len(C))].mean(0)
        t = T[rng.integers(0, len(T), len(T))].mean(0)
        yield c, t


def two_sided_p(X: np.ndarray) -> np.ndarray:
    """Two-sided resampling p of zero per column of X (resamples x items): twice the smaller tail share, each tail
    counted with +1 in numerator and denominator, capped at 1."""
    B = len(X)
    lo = ((X <= 0).sum(0) + 1) / (B + 1)
    hi = ((X >= 0).sum(0) + 1) / (B + 1)
    return np.minimum(1.0, 2 * np.minimum(lo, hi))


def bh(p: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values (q), in the input's order."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def side_counts(below: np.ndarray, n: int, full: np.ndarray) -> np.ndarray:
    """Recordings on the side of the full-data residual: those below the fit if it is negative, else the others."""
    below = np.asarray(below)
    return np.where(np.asarray(full) < 0, below, n - below)


def single_recordings(C: np.ndarray, T: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(per pair, the number of control recordings below the fit when each is fitted alone against the TNT mean;
    the number of TNT recordings whose fit puts the control mean below it)."""
    sc, st = C.mean(0), T.mean(0)
    n_ctrl = sum((iso_resid(C[r], st) < 0).astype(int) for r in range(len(C)))
    n_tnt = sum((iso_resid(sc, T[r]) < 0).astype(int) for r in range(len(T)))
    return np.asarray(n_ctrl), np.asarray(n_tnt)


def odor_design(pairs: Sequence[Tuple[str, str]], odors: Sequence[str]) -> np.ndarray:
    """Pairs x odors indicator of each pair's two odors."""
    col = {o: k for k, o in enumerate(odors)}
    X = np.zeros((len(pairs), len(odors)))
    for i, (a, b) in enumerate(pairs):
        X[i, col[a]] = X[i, col[b]] = 1
    return X


def odor_fit(res: np.ndarray, X: np.ndarray, fit_rows: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """Additive odor effects (least squares, pinv) fitted on `fit_rows` (default all), and every pair's residual
    beyond them (res - X @ a)."""
    rows = np.ones(len(res), bool) if fit_rows is None else np.asarray(fit_rows)
    a = np.linalg.pinv(X[rows]) @ res[rows]
    return a, res - X @ a


def percentiles(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    lo, hi = np.percentile(X, [5, 95], axis=0)
    return lo, hi


def pair_table(idx: Sequence[Tuple[str, str]], C: np.ndarray, T: np.ndarray, block: np.ndarray, *,
        n: int = RESAMPLES['pairs'], seed: int = SEEDS['pairs']) -> Tuple[pd.DataFrame, np.ndarray]:
    """Part 1: pairs.csv's table, and each resample's block mean residual."""
    grp = lambda o: H.ODOR_GROUPS.get(o, o)   # noqa: E731
    sc, st = C.mean(0), T.mean(0)
    res, drank = iso_resid(sc, st), rank_change(sc, st)
    R, D = [], []
    for c, t in resamples(C, T, n, seed):
        R.append(iso_resid(c, t))
        D.append(rank_change(c, t))
    R, D = np.array(R), np.array(D)
    r_lo, r_hi = percentiles(R)
    d_lo, d_hi = percentiles(D)
    r_p, d_p = two_sided_p(R), two_sided_p(D)
    n_ctrl_below, n_tnt_below = single_recordings(C, T)
    tab = pd.DataFrame(dict(odor_a=[a for a, _ in idx], odor_b=[b for _, b in idx],
        group_a=[grp(a) for a, _ in idx], group_b=[grp(b) for _, b in idx], block=block,
        within_group=[grp(a) == grp(b) for a, b in idx], sim_ctrl=sc, sim_tnt=st, apl_effect=st - sc,
        iso_resid=res, resid_lo=r_lo, resid_hi=r_hi, resid_p=r_p, resid_q=bh(r_p),
        resid_same_sign=(np.sign(R) == np.sign(res)).mean(0), n_ctrl_below=n_ctrl_below, n_tnt_below=n_tnt_below,
        n_ctrl_side=side_counts(n_ctrl_below, len(C), res), n_tnt_side=side_counts(n_tnt_below, len(T), res),
        rank_ctrl=rankdata(sc), rank_tnt=rankdata(st), d_rank=drank, drank_lo=d_lo, drank_hi=d_hi, drank_p=d_p,
        drank_q=bh(d_p)))
    return tab, R[:, block].mean(1)


def odor_tables(idx: Sequence[Tuple[str, str]], C: np.ndarray, T: np.ndarray, block: np.ndarray,
        odors: Sequence[str] = H.APL16_ODORS) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, tuple]]:
    """Parts 2 and 3: odors.csv's and block_pairs.csv's tables, and their summary rows (statistic -> (value, lo, hi))."""
    X = odor_design(idx, odors)
    M = X / X.sum(0)
    res = iso_resid(C.mean(0), T.mean(0))
    a, e = odor_fit(res, X)
    means, effects = [], []
    for c, t in resamples(C, T, RESAMPLES['odor_effects'], SEEDS['odor_effects']):
        r = iso_resid(c, t)
        means.append(r @ M)
        effects.append(np.linalg.pinv(X) @ r)
    m_lo, m_hi = percentiles(np.array(means))
    a_lo, a_hi = percentiles(np.array(effects))
    a3, rem = odor_fit(res, X, ~block)
    rems, rem_means, effects3 = [], [], []
    for c, t in resamples(C, T, RESAMPLES['outside_block'], SEEDS['outside_block']):
        r = iso_resid(c, t)
        a_b, rem_b = odor_fit(r, X, ~block)
        rems.append(rem_b[block])
        rem_means.append(rem_b[block].mean())
        effects3.append(a_b)
    s_lo, s_hi = percentiles(np.array(rems))
    o_lo, o_hi = percentiles(np.array(effects3))
    mean_lo, mean_hi = np.percentile(rem_means, [5, 95])
    odors_tab = pd.DataFrame(dict(odor=list(odors), group=[H.ODOR_GROUPS.get(o, o) for o in odors], mean_resid=res @ M,
        mean_lo=m_lo, mean_hi=m_hi, effect_all_pairs=a, effect_all_lo=a_lo, effect_all_hi=a_hi,
        effect_outside_block=a3, effect_out_lo=o_lo, effect_out_hi=o_hi))
    bidx = [p for p, b in zip(idx, block) if b]
    block_tab = pd.DataFrame(dict(odor_a=[a_ for a_, _ in bidx], odor_b=[b_ for _, b_ in bidx], iso_resid=res[block],
        predicted_from_odors=(X @ a3)[block], pair_specific=rem[block], specific_lo=s_lo, specific_hi=s_hi))
    rows = dict(block_predicted_from_odors=((X @ a3)[block].mean(), None, None),
        block_pair_specific=(rem[block].mean(), mean_lo, mean_hi),
        odor_effect_variance_share=(1 - e.var() / res.var(), None, None))
    return odors_tab, block_tab, rows


def load() -> Tuple[List[Tuple[str, str]], np.ndarray, np.ndarray]:
    """(pairs, control recordings x pairs, TNT recordings x pairs) of Yang's cosine similarities."""
    rec = {g: 1 - H.yang_recording_rdms(g, 'cosine') for g in GENOTYPES}
    assert list(rec['TNT3c'].columns) == list(rec['TNT'].columns), 'the genotypes list their pairs differently'
    C, T = rec['TNT3c'].to_numpy(), rec['TNT'].to_numpy()
    assert not np.isnan(C).any() and not np.isnan(T).any()
    return list(rec['TNT3c'].columns), C, T


def main() -> None:
    idx, C, T = load()
    block = np.array([H.block_label(a, b) == H.FRUITY_SMALL_BLOCK for a, b in idx])
    sc, st = C.mean(0), T.mean(0)
    res = iso_resid(sc, st)
    # the same fit on the distances, as the rest of the program computes it
    pidx = pd.MultiIndex.from_tuples(idx)
    d_res, d_info = H.isotonic_resid(pd.Series(1 - sc, index=pidx), pd.Series(1 - st, index=pidx))
    assert np.allclose(-res, d_res.reindex(pidx).to_numpy(), atol=1e-12), 'the residual disagrees with H.isotonic_resid'
    assert abs(d_info['r2'] - (1 - res.var(ddof=1) / sc.var(ddof=1))) < 1e-12
    pairs, block_means = pair_table(idx, C, T, block)
    odors, block_pairs, rows = odor_tables(idx, C, T, block)
    b_lo, b_hi = np.percentile(block_means, [5, 95])
    excl = (pairs.resid_lo > 0) | (pairs.resid_hi < 0)
    stats = dict(iso_r2=(d_info['r2'], None, None), resid_sd=(d_info['resid_sd'], None, None),
        block_mean_resid=(res[block].mean(), b_lo, b_hi), **rows,
        n_pairs_interval_excludes_zero=(int(excl.sum()), None, None),
        n_pairs_below_interval=(int((pairs.resid_hi < 0).sum()), None, None),
        n_pairs_above_interval=(int((pairs.resid_lo > 0).sum()), None, None),
        **{'n_pairs_q_lt_0.10': (int((pairs.resid_q < 0.10).sum()), None, None),
            'n_pairs_q_lt_0.05': (int((pairs.resid_q < 0.05).sum()), None, None)},
        n_pairs_drank_interval_excludes_zero=(int(((pairs.drank_lo > 0) | (pairs.drank_hi < 0)).sum()), None, None),
        n_control_recordings=(len(C), None, None), n_tnt_recordings=(len(T), None, None))
    summary = pd.DataFrame([dict(statistic=k, value=v, lo=lo, hi=hi) for k, (v, lo, hi) in stats.items()])
    out = H.shared_dir() / 'results' / 'yang_pair_deviations'
    out.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out / 'pairs.csv', index=False)
    odors.to_csv(out / 'odors.csv', index=False)
    block_pairs.to_csv(out / 'block_pairs.csv', index=False)
    summary.to_csv(out / 'summary.csv', index=False)
    rdm_csv = H.yang_data_dir() / 'computed' / 'rdm_per_recording_long.csv'
    prov = dict(script='analyses/yang_pair_deviations.py', **H.code_commits(),
        inputs=dict(per_recording_rdms=str(rdm_csv), genotypes=dict(control='TNT3c', apl_output_blocked='TNT'),
            measure='1 - cosine distance (cosine similarity) over all KCs of a recording, mean over recordings'),
        seeds=SEEDS, resamples=RESAMPLES, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        host=subprocess.run(['hostname'], capture_output=True, text=True).stdout.strip(),
        question="Yang, 2026-11-16: which odor pairs deviate most, and most consistently, from the monotone fit of the "
            "control flies' KC matrix on the TNT matrix; whether all nine block pairs do, and which others")
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    pd.set_option('display.width', 220)
    print(summary.to_string(index=False))
    print(f'written to {out}')


if __name__ == '__main__':
    main()
