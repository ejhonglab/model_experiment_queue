#!/usr/bin/env python3
"""E128: the flies' bouton totals by tuning group (no simulations, no olfsysm change).

The entry is item E128 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-22 12:41; moved in under the idle-host rule 12:56). The question: whether a few groups of boutons make the PN
boutons' per-odor totals uneven across odors in Remy's five bouton flies (2022-10-21 fly02 excluded, E119's part 0),
and whether the model's boutons, sampled as unevenly as the flies' groups are, give the flies' totals. "The family" is
the breadth-2 APL-effect fit on E72 part B's eight wirings (E90's saved stages, APL on).

Statistics (E119's, `E119.roi_stats`): the per-odor totals' CV over the 12 megamat odors (the CV over odors of the mean
response per odor over the ROIs, signed), the responders' CV and the top-1 % share; the flies' trial means in
baseline-noise units, responders above 2 SDs (E119's VERDICT_CRITERION).

Groups (E119's A4, `E119.a4_groups`): a fly's responsive ROIs (any odor above 2 SDs) grouped by their 12-odor tuning
(average linkage on 1 - Pearson r, cut at 0.3); every group counts, whatever its size; the ROIs in no group (the
unresponsive ones) are one remainder, kept in every reading and never left out.

Parts (CLI `--part`):
- part0: each fly's totals' CV, responders' CV and top-1 % share are reproduced against E119's a2_roi_stats.csv and its
  tuning groups (group id and size of every group of >= 5 ROIs; the numbers of responsive ROIs, groups and big groups)
  against a4_groups.csv / a4_summary.csv. Any difference stops the run (part0.json says which).
- A, the flies: per fly the totals' CV as measured; with every group at equal weight (the total the mean of the groups'
  mean responses, the remainder one more group; without the remainder beside); and with the 1, 2, 3 and 5 groups of
  largest variance share left out. A group's share is cov(group sum, total) / var(total) across odors (they sum to 1
  with the remainder's, and may be negative). The null: the responsive ROIs dealt at random into groups of the same
  sizes (the remainder fixed), N_NULL deals per fly, the same readings (shares recomputed on each deal). The repeat
  check: groups (and the groups to leave out) made from the mean of two repeats, the readings taken on the third, each
  repeat held out in turn.
- B, the model: E90's saved boutons of the 28 shared glomeruli (E119's `model_units`, 'rate' reading, unpooled), eight
  wirings, APL on. Each of a fly's groups is dealt to one of the 28 glomeruli at random with replacement; a glomerulus's
  weight is the number of ROIs dealt to it; its ROIs are that many of its boutons drawn at random with replacement.
  Beside: the same with only the groups of >= 5 ROIs; the glomeruli at equal weight (EQUAL_PER_GLOM boutons drawn from
  each, EQUAL_DRAWS draws); E119's own reading (every bouton once). The threshold is matched to the fly's overall
  response rate among its grouped ROIs (`E119.matched_threshold`).
- score: the classes, the outcome letter (verdict.json, comparison_table.md, figure.png).

The classes (the entry's, fixed beforehand; half-way 0.230 between the family's 0.103 and the flies' lower end 0.357):
- flies, "a few groups carry it": in at least 4 of 5 flies, equal weights or leaving out up to three groups brings the
  totals' CV below 0.230 (the smallest of the four readings), and the null's deals do so (the smallest of the same four)
  in under 5 %;
- model, "the sample can account for it": on at least 6 of 8 wirings, at least 5 % of deals have the totals' CV inside
  the flies' range (E119's a2 range at 2 SDs) with neither the responders' CV nor the top-1 % share above the flies'
  range; a wiring's fraction is over all five flies' deals pooled (per fly beside).
- the repeat check holds if, in at least 4 of 5 flies, the held-out readings bring the totals' CV below 0.230 for at
  least 2 of the 3 held-out repeats.
Outcomes: flies fail -> (c); flies hold, the repeat check fails -> (d); both hold and the model holds -> (a); else (b).

Outputs (<shared>/results/E128_bouton_groups/): part0.json, a_readings.csv, a_shares.csv, a_null.csv,
a_null_summary.csv, a_repeat.csv, b_deals.csv (summary per fly, wiring, variant), b_wirings.csv, b_leads.csv,
verdict.json, comparison_table.md, figure.png, provenance.json; the deals themselves under outputs/E128_bouton_groups/.
The write-up checks (analyses/E128_writeup_*.py) add writeup_*.csv and provenance_writeup_*.json (`save_writeup`).

Usage (from the repository root; light, no simulation):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E128_bouton_groups.py --part all
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E128_bouton_groups'
CRITERION = 2.0
HALF_WAY = 0.230                         # between the family's 0.103 and the flies' lower end 0.357
LEAVE_OUT = (1, 2, 3, 5)
CLASS_READINGS = ('equal', 'leave1', 'leave2', 'leave3')   # "equal weights, or leaving out up to three groups"
N_NULL = 200
N_DEALS = 1000
NULL_FRAC = 0.05
MODEL_FRAC = 0.05
FLY_MIN = 4                              # of 5 flies
WIRING_MIN = 6                           # of 8 wirings
REPEAT_MIN = 2                           # of 3 held-out repeats
MIN_GROUP = 5
EQUAL_PER_GLOM = 30
EQUAL_DRAWS = 50
SEED = 0
TOL = 1e-9
THREE = ('totals_cv', 'resp_cv', 'top1')


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir() -> Path:
    return H.outputs_dir(EXPERIMENT)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def save_writeup(tables: Dict[str, pd.DataFrame], script: str, out: Optional[Path] = None, **extra) -> Path:
    """Saves the tables of a write-up check (analyses/E128_writeup_*.py; no simulation).

    Args:
        tables: {name: table}, each written as <out>/writeup_<name>.csv, without the index.
        script: the check's path from the repository root (`H.script_path(__file__)`).
        out: the folder; by default the results folder (`out_dir()`).
        **extra: more entries for the provenance file.

    Returns:
        The provenance file written beside the tables, provenance_writeup_<name>.json (<name> the script's name after
        "E128_writeup_"), with the code's commits and the tables' names.
    """
    out = Path(out) if out is not None else out_dir()
    out.mkdir(parents=True, exist_ok=True)
    for name, T in tables.items():
        T.to_csv(out / f'writeup_{name}.csv', index=False)
    p = out / f"provenance_writeup_{Path(script).stem.replace('E128_writeup_', '')}.json"
    _write_json(p, dict(script=script, **H.code_commits(), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
                        tables=[f'writeup_{n}.csv' for n in tables], **extra))
    return p


def _e119():
    from experiments import E119_stages_per_roi as E119
    return E119


# ---------------------------------------------------------------------------
# groups and readings

def fly_groups(M: np.ndarray, criterion: float = CRITERION, cut: float = 0.3) -> np.ndarray:
    """E119 A4's tuning groups of one fly's trial means M (ROIs x odors): per ROI its fcluster label (1, 2, ...) among
    the responsive ROIs (any odor above `criterion`), 0 for the remainder (the ROIs in no group)."""
    from scipy.cluster.hierarchy import fcluster, linkage
    M = np.asarray(M, dtype=float)
    keep = (M > criterion).any(axis=1)
    lab = np.zeros(M.shape[0], dtype=int)
    if keep.sum() >= 2:
        lab[keep] = fcluster(linkage(M[keep], 'average', metric='correlation'), cut, 'distance')
    elif keep.sum() == 1:
        lab[keep] = 1
    return lab


def _group_sums(M: np.ndarray, labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(group ids, groups x odors summed responses, sizes); the remainder (label 0) included when present."""
    ids, inv = np.unique(labels, return_inverse=True)
    S = np.zeros((len(ids), M.shape[1]))
    np.add.at(S, inv, M)
    return ids, S, np.bincount(inv, minlength=len(ids))


def variance_shares(M: np.ndarray, labels: np.ndarray) -> pd.Series:
    """Each group's share of the variance across odors of the summed total (remainder, label 0, included):
    cov(group sum, total) / var(total). The shares sum to 1."""
    ids, S, _ = _group_sums(np.asarray(M, dtype=float), labels)
    T = S.sum(axis=0)
    Tc = T - T.mean()
    v = float(Tc @ Tc)
    return pd.Series((S - S.mean(axis=1, keepdims=True)) @ Tc / v if v > 0 else np.full(len(ids), np.nan), index=ids)


def cv(t: np.ndarray) -> float:
    """E119's totals' CV: the SD over odors (ddof 0) over the mean."""
    t = np.asarray(t, dtype=float)
    return float(t.std() / t.mean()) if t.mean() != 0 else np.nan


def readings(M: np.ndarray, labels: np.ndarray, order: Optional[Sequence[int]] = None,
             leave_out: Sequence[int] = LEAVE_OUT) -> Dict[str, float]:
    """The totals' CV of ROIs x odors M: 'measured'; 'equal' (every group at equal weight, the remainder one more
    group); 'equal_no_rem' (the groups only); 'leave<k>' (the k groups first in `order`, by default the largest variance
    shares of M itself, left out; the remainder never left out)."""
    M = np.asarray(M, dtype=float)
    ids, S, n = _group_sums(M, labels)
    out = dict(measured=cv(S.sum(axis=0)), equal=cv((S / n[:, None]).mean(axis=0)))
    g = ids != 0
    out['equal_no_rem'] = cv((S[g] / n[g, None]).mean(axis=0)) if g.any() else np.nan
    # beside: only the groups of >= MIN_GROUP ROIs at equal weight, the smaller groups' ROIs pooled with the remainder
    # as one more group (no singleton carries a group's weight)
    big = g & (n >= MIN_GROUP)
    if big.any():
        parts = list(S[big] / n[big, None])
        if (~big).any():
            parts.append(S[~big].sum(axis=0) / n[~big].sum())
        out['equal_big'] = cv(np.mean(parts, axis=0))
    else:
        out['equal_big'] = np.nan
    if order is None:
        sh = variance_shares(M, labels)
        order = list(sh.drop(0, errors='ignore').sort_values(ascending=False).index)
    pos = {k: i for i, k in enumerate(ids)}
    order = [o for o in order if o != 0 and o in pos]
    for k in leave_out:
        drop = [pos[o] for o in order[:k]]
        keep = np.setdiff1d(np.arange(len(ids)), drop)
        out[f'leave{k}'] = cv(S[keep].sum(axis=0))
    return out


def class_min(r: Dict[str, float]) -> float:
    """The smallest of the class's readings (equal weight, leaving out 1 to 3 groups)."""
    return float(np.nanmin([r[k] for k in CLASS_READINGS]))


def three_stats(X: np.ndarray, thr: float) -> Dict[str, float]:
    """E119's totals' CV, responders' CV and top-1 % share of ROIs x odors X at threshold `thr` (as `E119.roi_stats`)."""
    X = np.asarray(X, dtype=float)
    resp = X[X > thr]
    k = max(1, int(np.ceil(0.01 * X.shape[0])))
    P = np.clip(X, 0, None)
    s = P.sum(axis=0)
    top = -np.sort(-P, axis=0)[:k].sum(axis=0)
    ok = s > 0
    return dict(totals_cv=cv(X.mean(axis=0)),
                resp_cv=float(resp.std() / resp.mean()) if resp.size > 1 and resp.mean() > 0 else np.nan,
                top1=float(np.mean(top[ok] / s[ok])) if ok.any() else np.nan)


# ---------------------------------------------------------------------------
# the flies

def fly_means() -> Dict[str, np.ndarray]:
    """Per bouton fly, the trials (repeats x ROIs x odors, baseline-noise units; E119's `fly_trials`)."""
    return {a: d['noise'] for a, d in _e119().fly_trials('bouton').items()}


def part0() -> dict:
    """Reproduces E119's a2 statistics (criterion 2, the flies' rows) and A4's tuning groups per fly; raises on any
    difference (written to part0.json first)."""
    E119 = _e119()
    r119 = H.shared_dir() / 'results' / 'E119_stages_per_roi'
    a2 = pd.read_csv(r119 / 'a2_roi_stats.csv')
    a2 = a2[(a2.reading == 'flies') & (a2.criterion == CRITERION)].set_index('wiring')
    g4 = pd.read_csv(r119 / 'a4_groups.csv')
    s4 = pd.read_csv(r119 / 'a4_summary.csv').set_index('acq')
    fs = E119.fly_roi_stats('bouton')[CRITERION]
    out, diffs = {}, []
    for acq, X in fly_means().items():
        M = X.mean(axis=0)
        st = fs[acq]
        mine = three_stats(M, CRITERION)
        row = a2.loc[f'fly:{acq}']
        d = {s: float(abs(st[s] - row[s])) for s in THREE}
        d.update({f'{s}_three_stats': float(abs(mine[s] - st[s])) for s in THREE})
        lab = fly_groups(M)
        ids, cnt = np.unique(lab[lab > 0], return_counts=True)
        big = {int(i): int(c) for i, c in zip(ids, cnt) if c >= MIN_GROUP}
        ref = {int(r.group): int(r.n_rois) for r in g4[g4.acq == acq].itertuples()}
        sref = s4.loc[acq]
        ok_groups = (big == ref and int((lab > 0).sum()) == int(sref.n_responsive)
                     and len(ids) == int(sref.n_tuning_groups) and len(big) == int(sref.n_big_groups))
        ok = ok_groups and max(d.values()) <= TOL
        out[acq] = dict(stat_diffs=d, groups_equal=ok_groups, n_rois=int(M.shape[0]), n_responsive=int((lab > 0).sum()),
                        n_groups=int(len(ids)), n_big_groups=len(big), passed=bool(ok),
                        **{s: float(st[s]) for s in THREE})
        if not ok:
            diffs.append(acq)
    res = dict(passed=not diffs, failed=diffs, flies=out, tolerance=TOL)
    _write_json(out_dir() / 'part0.json', res)
    if diffs:
        raise RuntimeError(f'part 0 does not reproduce E119 for {diffs}; see part0.json')
    return res


def null_deal(labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """The responsive ROIs' labels permuted among them (groups of the same sizes), the remainder fixed."""
    lab = labels.copy()
    g = lab > 0
    lab[g] = rng.permutation(lab[g])
    return lab


def part_a(n_null: int = N_NULL, seed: int = SEED, flies: Optional[Dict[str, np.ndarray]] = None) -> dict:
    """The flies' readings, the groups' variance shares, the null and the repeat check (a_*.csv)."""
    flies = flies if flies is not None else fly_means()
    rows, shares, nulls, reps = [], [], [], []
    for fi, (acq, X) in enumerate(sorted(flies.items())):
        M = X.mean(axis=0)
        lab = fly_groups(M)
        sh = variance_shares(M, lab)
        _, _, n = _group_sums(M, lab)
        for (gid, s), size in zip(sh.items(), n):
            shares.append(dict(acq=acq, group=int(gid), n_rois=int(size), share=float(s), remainder=bool(gid == 0)))
        r = readings(M, lab)
        top = sh.drop(0, errors='ignore').sort_values(ascending=False)
        rows.append(dict(acq=acq, n_rois=int(M.shape[0]), n_grouped=int((lab > 0).sum()),
                         n_groups=int(len(np.unique(lab[lab > 0]))),
                         remainder_share=float(sh.get(0, np.nan)), **r, class_min=class_min(r),
                         below=bool(class_min(r) < HALF_WAY),
                         **{f'top{k}_share': float(top.iloc[:k].sum()) for k in LEAVE_OUT}))
        rng = np.random.default_rng([seed, fi])
        for d in range(n_null):
            rn = readings(M, null_deal(lab, rng))
            nulls.append(dict(acq=acq, deal=d, **rn, class_min=class_min(rn)))
        for h in range(X.shape[0]):
            Mtr = np.delete(X, h, axis=0).mean(axis=0)
            lt = fly_groups(Mtr)
            sht = variance_shares(Mtr, lt)
            order = list(sht.drop(0, errors='ignore').sort_values(ascending=False).index)
            rh = readings(X[h], lt, order=order)
            rt = readings(Mtr, lt, order=order)
            reps.append(dict(acq=acq, held_out=h, n_groups=int(len(np.unique(lt[lt > 0]))), **rh,
                             class_min=class_min(rh), below=bool(class_min(rh) < HALF_WAY),
                             train_measured=rt['measured'], train_class_min=class_min(rt)))
    A = pd.DataFrame(rows)
    N = pd.DataFrame(nulls)
    ns = N.groupby('acq').agg(**{f'{k}_mean': (k, 'mean') for k in ('measured', 'equal', 'equal_big') + tuple(f'leave{k}' for k in LEAVE_OUT)},
                              class_min_mean=('class_min', 'mean'), class_min_p05=('class_min', lambda x: x.quantile(0.05)),
                              frac_below=('class_min', lambda x: float((x < HALF_WAY).mean()))).reset_index()
    A = A.merge(ns[['acq', 'frac_below']].rename(columns={'frac_below': 'null_frac_below'}), on='acq')
    A['fly_passes'] = A.below & (A.null_frac_below < NULL_FRAC)
    R = pd.DataFrame(reps)
    out = dict(readings=A, shares=pd.DataFrame(shares), null=N, null_summary=ns, repeat=R)
    return out


# ---------------------------------------------------------------------------
# the model

def model_boutons(wiring: str) -> Tuple[np.ndarray, np.ndarray]:
    """A wiring's boutons of the 28 shared glomeruli, APL on, 'rate' reading (the window's change, E119's
    `model_units`), as (boutons x the 12 megamat odors, glomerulus per bouton)."""
    E119 = _e119()
    u = E119.model_units('bouton', wiring, 0.0)
    return u['on'][list(E119.ODORS)].to_numpy(dtype=float), np.asarray(u['groups'])


def deal_population(U: np.ndarray, gl: np.ndarray, sizes: Sequence[int], rng: np.random.Generator,
                    gloms: Optional[np.ndarray] = None) -> Tuple[np.ndarray, Dict[str, int]]:
    """One deal: each group (of `sizes`) to one of the glomeruli at random with replacement; a glomerulus's weight is the
    ROIs dealt to it, its ROIs that many of its boutons drawn with replacement. Returns (ROIs x odors, weights)."""
    gloms = np.unique(gl) if gloms is None else gloms
    to = rng.integers(0, len(gloms), size=len(sizes))
    w = np.bincount(to, weights=np.asarray(sizes, dtype=float), minlength=len(gloms)).astype(int)
    parts = []
    for g, k in zip(gloms, w):
        if k:
            idx = np.flatnonzero(gl == g)
            parts.append(U[rng.choice(idx, size=k, replace=True)])
    return np.vstack(parts), {str(g): int(k) for g, k in zip(gloms, w) if k}


def equal_population(U: np.ndarray, gl: np.ndarray, rng: np.random.Generator, per: int = EQUAL_PER_GLOM) -> np.ndarray:
    """The glomeruli at equal weight: `per` boutons of each drawn with replacement."""
    return np.vstack([U[rng.choice(np.flatnonzero(gl == g), size=per, replace=True)] for g in np.unique(gl)])


def fly_rate(M: np.ndarray, labels: np.ndarray, criterion: float = CRITERION) -> float:
    """The fly's overall response rate among its grouped ROIs."""
    return float((M[labels > 0] > criterion).mean())


def part_b(n_deals: int = N_DEALS, seed: int = SEED, wirings: Optional[Sequence[str]] = None,
           flies: Optional[Dict[str, np.ndarray]] = None, spread: Optional[dict] = None) -> dict:
    """The deals per fly, wiring and variant ('all' groups, 'big' groups of >= 5 ROIs), the equal-weight and E119
    readings beside (b_*.csv)."""
    E119 = _e119()
    wirings = tuple(wirings or E119.WIRINGS)
    flies = flies if flies is not None else fly_means()
    sp = spread or flies_range()
    deals, eq, leads = [], [], []
    fly_info = {}
    for acq, X in sorted(flies.items()):
        M = X.mean(axis=0)
        lab = fly_groups(M)
        sizes = np.unique(lab[lab > 0], return_counts=True)[1]
        fly_info[acq] = dict(sizes=sizes, rate=fly_rate(M, lab))
    for wi, w in enumerate(wirings):
        U, gl = model_boutons(w)
        gloms = np.unique(gl)
        rng = np.random.default_rng([seed, 1000 + wi])
        reps = []
        for _ in range(EQUAL_DRAWS):
            P = equal_population(U, gl, rng)
            reps.append(three_stats(P, E119.matched_threshold(pd.DataFrame(P, columns=list(E119.ODORS)),
                                                              float(np.mean([f['rate'] for f in fly_info.values()])))))
        eq.append(dict(wiring=w, variant='equal', **pd.DataFrame(reps).mean().to_dict()))
        thr = E119.matched_threshold(pd.DataFrame(U, columns=list(E119.ODORS)),
                                     float(np.mean([f['rate'] for f in fly_info.values()])))
        eq.append(dict(wiring=w, variant='e119_every_bouton', **three_stats(U, thr)))
        for fi, (acq, f) in enumerate(sorted(fly_info.items())):
            for variant, sizes in (('all', f['sizes']), ('big', f['sizes'][f['sizes'] >= MIN_GROUP])):
                rng = np.random.default_rng([seed, wi, fi, 0 if variant == 'all' else 1])
                for d in range(n_deals):
                    P, wts = deal_population(U, gl, sizes, rng, gloms)
                    thr = float(np.quantile(P.ravel(), 1.0 - f['rate']))
                    st = three_stats(P, thr)
                    ok = model_pass(st, sp)
                    deals.append(dict(acq=acq, wiring=w, variant=variant, deal=d, n_rois=int(P.shape[0]),
                                      n_gloms=len(wts), **st, passes=ok))
                    if variant == 'all':
                        topg = max(wts, key=wts.get)
                        leads.append(dict(acq=acq, wiring=w, top_glomerulus=topg,
                                          top_weight_share=wts[topg] / P.shape[0], passes=ok))
        print(f'part B: {w} done', flush=True)
    D = pd.DataFrame(deals)
    L = pd.DataFrame(leads)
    lead = (L.groupby('top_glomerulus').agg(n_deals=('passes', 'size'), n_passing=('passes', 'sum'))
            .assign(share_of_passing=lambda d: d.n_passing / max(1, d.n_passing.sum()),
                    share_of_deals=lambda d: d.n_deals / d.n_deals.sum())
            .sort_values('n_passing', ascending=False).reset_index())
    return dict(deals=D, equal=pd.DataFrame(eq), leads=lead)


def flies_range() -> dict:
    """The flies' range of the three statistics at 2 SDs (E119's a2_spread.csv), and on the grouped ROIs only (beside)."""
    r119 = H.shared_dir() / 'results' / 'E119_stages_per_roi'
    s = pd.read_csv(r119 / 'a2_spread.csv')
    s = s[(s.stage == 'bouton') & (s.criterion == CRITERION)].iloc[0]
    return {k: (float(s[f'{k}_lo']), float(s[f'{k}_hi'])) for k in THREE}


def grouped_range(flies: Optional[Dict[str, np.ndarray]] = None) -> dict:
    """The flies' range of the three statistics on their grouped ROIs only (the remainder out), at 2 SDs."""
    flies = flies if flies is not None else fly_means()
    v = []
    for X in flies.values():
        M = X.mean(axis=0)
        v.append(three_stats(M[fly_groups(M) > 0], CRITERION))
    V = pd.DataFrame(v)
    return {k: (float(V[k].min()), float(V[k].max())) for k in THREE}


def model_pass(st: Dict[str, float], sp: dict) -> bool:
    """A deal passes: the totals' CV inside the flies' range and neither amplitude statistic above the flies' range."""
    lo, hi = sp['totals_cv']
    return bool(lo <= st['totals_cv'] <= hi and st['resp_cv'] <= sp['resp_cv'][1] and st['top1'] <= sp['top1'][1])


def wiring_table(D: pd.DataFrame, sp_grouped: Optional[dict] = None, sp: Optional[dict] = None) -> pd.DataFrame:
    """Per wiring and variant: the fraction of deals passing (pooled over flies; per fly beside), the totals' CV's
    quantiles; with the grouped-ROI range's passing fraction beside."""
    rows = []
    lo, hi = (sp or flies_range())['totals_cv']
    for (w, v), sub in D.groupby(['wiring', 'variant'], sort=False):
        tin = (sub.totals_cv >= lo) & (sub.totals_cv <= hi)
        r = dict(wiring=w, variant=v, frac_pass=float(sub.passes.mean()),
                 frac_totals_in=float(tin.mean()), frac_amp_ok_given_totals_in=float(sub.passes[tin].mean()) if tin.any()
                 else np.nan,
                 totals_cv_median=float(sub.totals_cv.median()), totals_cv_p95=float(sub.totals_cv.quantile(0.95)),
                 totals_cv_max=float(sub.totals_cv.max()), resp_cv_median=float(sub.resp_cv.median()),
                 top1_median=float(sub.top1.median()),
                 n_flies_pass=int((sub.groupby('acq').passes.mean() >= MODEL_FRAC).sum()))
        for acq, s2 in sub.groupby('acq'):
            r[f'frac_pass__{acq}'] = float(s2.passes.mean())
        if sp_grouped is not None:
            r['frac_pass_grouped_range'] = float(np.mean([model_pass(dict(totals_cv=a, resp_cv=b, top1=c), sp_grouped)
                                                          for a, b, c in sub[list(THREE)].to_numpy()]))
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# score

def classes(A: pd.DataFrame, R: pd.DataFrame, W: pd.DataFrame) -> dict:
    """The flies' class, the repeat check, the model's class (on the 'all' variant) and the outcome letter."""
    n_fly = int(A.fly_passes.sum())
    flies_hold = n_fly >= FLY_MIN
    rp = R.groupby('acq').below.sum()
    n_rep = int((rp >= REPEAT_MIN).sum())
    repeat_holds = n_rep >= FLY_MIN
    Wa = W[W.variant == 'all']
    n_w = int((Wa.frac_pass >= MODEL_FRAC).sum())
    model_holds = n_w >= WIRING_MIN
    Wb = W[W.variant == 'big']
    n_wb = int((Wb.frac_pass >= MODEL_FRAC).sum())
    carried = {k: int((A[k] < HALF_WAY).sum()) for k in CLASS_READINGS + ('equal_big', 'leave5')}
    letter = 'c' if not flies_hold else ('d' if not repeat_holds else ('a' if model_holds else 'b'))
    return dict(flies_class=dict(holds=flies_hold, n_flies=n_fly, of=len(A), need=FLY_MIN,
                                 n_flies_below_by_reading=carried),
                repeat_check=dict(holds=repeat_holds, n_flies=n_rep, of=len(rp), need=FLY_MIN, per_fly=rp.to_dict()),
                model_class=dict(holds=model_holds, n_wirings=n_w, of=len(Wa), need=WIRING_MIN, variant='all'),
                model_class_big_groups=dict(holds=n_wb >= WIRING_MIN, n_wirings=n_wb, of=len(Wb)),
                outcome=letter)


def _f(x, spec='.3f'):
    return '-' if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(x, spec)


def comparison_table(v: dict, A: pd.DataFrame, ns: pd.DataFrame, R: pd.DataFrame, W: pd.DataFrame,
                     E: pd.DataFrame, sp: dict, spg: dict) -> str:
    L = ['# E128: the flies\' bouton totals by tuning group', '',
         f"Outcome ({v['outcome']}). Flies' class {'holds' if v['flies_class']['holds'] else 'fails'} "
         f"({v['flies_class']['n_flies']} of {v['flies_class']['of']} flies, need {FLY_MIN}); repeat check "
         f"{'holds' if v['repeat_check']['holds'] else 'fails'} ({v['repeat_check']['n_flies']} of "
         f"{v['repeat_check']['of']}); model's class {'holds' if v['model_class']['holds'] else 'fails'} "
         f"({v['model_class']['n_wirings']} of {v['model_class']['of']} wirings, need {WIRING_MIN}; groups of >= 5 ROIs "
         f"only: {v['model_class_big_groups']['n_wirings']}).", '',
         "Flies below 0.230 by reading: " + ', '.join(f'{k} {n}' for k, n in
                                                     v['flies_class']['n_flies_below_by_reading'].items()) + '.', '',
         f'Half-way point {HALF_WAY}; the flies\' totals\' CV range (E119, 2 SDs) {_f(sp["totals_cv"][0])} to '
         f'{_f(sp["totals_cv"][1])}.', '',
         '## Part A: the flies (totals\' CV over the 12 odors)', '',
         '| fly | groups | measured | equal weight | equal, no remainder | equal, groups >= 5 | leave 1 | leave 2 | leave 3 | leave 5 | '
         'top-3 share | remainder share | null below 0.230 | passes |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in A.itertuples():
        L.append(f'| {r.acq} | {r.n_groups} | {_f(r.measured)} | {_f(r.equal)} | {_f(r.equal_no_rem)} | {_f(r.equal_big)} | '
                 f'{_f(r.leave1)} | {_f(r.leave2)} | {_f(r.leave3)} | {_f(r.leave5)} | {_f(r.top3_share, ".2f")} | '
                 f'{_f(r.remainder_share, ".2f")} | {_f(r.null_frac_below, ".3f")} | {"yes" if r.fly_passes else "no"} |')
    L += ['', 'Null (ROIs dealt into groups of the same sizes, 200 deals per fly): means.', '',
          '| fly | measured | equal | equal, groups >= 5 | leave 1 | leave 3 | smallest of the four, 5th percentile |',
          '|---|---:|---:|---:|---:|---:|---:|']
    for r in ns.itertuples():
        L.append(f'| {r.acq} | {_f(r.measured_mean)} | {_f(r.equal_mean)} | {_f(r.equal_big_mean)} | {_f(r.leave1_mean)} | '
                 f'{_f(r.leave3_mean)} | {_f(r.class_min_p05)} |')
    L += ['', '## The repeat check (groups from two repeats, read on the third)', '',
          '| fly | held out | groups | measured (third) | equal | equal, groups >= 5 | leave 1 | leave 3 | below 0.230 |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in R.itertuples():
        L.append(f'| {r.acq} | {r.held_out} | {r.n_groups} | {_f(r.measured)} | {_f(r.equal)} | {_f(r.equal_big)} | '
                 f'{_f(r.leave1)} | {_f(r.leave3)} | {"yes" if r.below else "no"} |')
    L += ['', '## Part B: the model\'s boutons dealt as the flies\' groups (APL on, rate reading)', '',
          f'Flies\' range: totals\' CV {_f(sp["totals_cv"][0])} to {_f(sp["totals_cv"][1])}, responders\' CV up to '
          f'{_f(sp["resp_cv"][1])}, top-1 % share up to {_f(sp["top1"][1])} (on grouped ROIs only: '
          f'{_f(spg["totals_cv"][0])} to {_f(spg["totals_cv"][1])}, {_f(spg["resp_cv"][1])}, {_f(spg["top1"][1])}).', '',
          '| wiring | variant | deals passing | totals\' CV in range | of those, amplitude not above | '
          'flies with >= 5 % | totals\' CV median | 95th pct | responders\' CV median | '
          'top-1 % median | passing, grouped-ROI range |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in W.itertuples():
        L.append(f'| {r.wiring} | {r.variant} | {_f(r.frac_pass)} | {_f(r.frac_totals_in)} | '
                 f'{_f(r.frac_amp_ok_given_totals_in, ".2f")} | {r.n_flies_pass} | {_f(r.totals_cv_median)} | '
                 f'{_f(r.totals_cv_p95)} | {_f(r.resp_cv_median)} | {_f(r.top1_median)} | '
                 f'{_f(r.frac_pass_grouped_range)} |')
    L += ['', 'Beside: the glomeruli at equal weight, and E119\'s reading (every bouton once).', '',
          '| wiring | variant | totals\' CV | responders\' CV | top-1 % |', '|---|---|---:|---:|---:|']
    for r in E.itertuples():
        L.append(f'| {r.wiring} | {r.variant} | {_f(r.totals_cv)} | {_f(r.resp_cv)} | {_f(r.top1)} |')
    return '\n'.join(L) + '\n'


def figure(A: pd.DataFrame, N: pd.DataFrame, D: pd.DataFrame, sp: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    keys = ['measured', 'equal', 'leave1', 'leave2', 'leave3', 'leave5']
    x = np.arange(len(keys))
    for i, r in enumerate(A.itertuples()):
        ax[0].plot(x + (i - 2) * 0.06, [getattr(r, k) for k in keys], 'o-', ms=4, label=r.acq.split('__megamat')[0])
        nm = N[N.acq == r.acq]
        ax[0].plot(x + (i - 2) * 0.06, [nm[k].mean() for k in keys], ':', color=ax[0].lines[-1].get_color(), lw=1)
    ax[0].axhline(HALF_WAY, color='k', lw=0.8, ls='--')
    ax[0].axhspan(*sp['totals_cv'], color='0.85', zorder=0)
    ax[0].set_xticks(x, ['measured', 'equal\nweight', 'leave 1', 'leave 2', 'leave 3', 'leave 5'])
    ax[0].set_ylabel("totals' CV over the 12 odors")
    ax[0].set_title("Flies' boutons by tuning group\n(dotted: ROIs dealt at random into the same sizes)")
    ax[0].legend(fontsize=7)
    Da = D[D.variant == 'all']
    ws = list(dict.fromkeys(Da.wiring))
    ax[1].boxplot([Da[Da.wiring == w].totals_cv for w in ws], showfliers=False)
    ax[1].set_xticks(np.arange(1, len(ws) + 1), ws)
    ax[1].axhspan(*sp['totals_cv'], color='0.85', zorder=0, label="flies' range")
    ax[1].axhline(HALF_WAY, color='k', lw=0.8, ls='--', label='half-way 0.230')
    ax[1].set_ylabel("totals' CV over the 12 odors")
    ax[1].set_title("Model boutons dealt as each fly's groups\n(APL on, 5 flies x 1000 deals per wiring)")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def score(a: dict, b: dict) -> dict:
    """The classes and the outcome letter (verdict.json, comparison_table.md, figure.png) from parts A and B."""
    od = out_dir()
    sp, spg = flies_range(), grouped_range()
    D = b['deals']
    W = wiring_table(D, spg, sp)
    W.to_csv(od / 'b_wirings.csv', index=False)
    v = classes(a['readings'], a['repeat'], W)
    v.update(half_way=HALF_WAY, flies_range=sp, flies_range_grouped=spg,
             rules_fixed_here=[
                 "equal weight: the total is the mean of the groups' mean responses, the remainder one more group "
                 "(without it beside)",
                 "the flies' class reads the smallest of the four readings (equal weight, leave 1, 2, 3) against 0.230; "
                 "the null's deal counts if the smallest of the same four is below 0.230",
                 'every tuning group counts (any size); the remainder is the unresponsive ROIs',
                 'repeat check: groups and the groups to leave out from the mean of two repeats, readings on the third; '
                 'holds if 4 of 5 flies are below 0.230 on at least 2 of the 3 held-out repeats',
                 "part B: a dealt ROI is one of the glomerulus's boutons drawn with replacement; only the groups are "
                 'dealt (the remainder is not); threshold matched to the fly\'s own overall response rate among its '
                 'grouped ROIs',
                 "part B: a wiring's fraction of passing deals is over the five flies' deals pooled; 'above the flies'' "
                 "is above the top of their range",
             ])
    _write_json(od / 'verdict.json', v)
    (od / 'comparison_table.md').write_text(comparison_table(v, a['readings'], a['null_summary'], a['repeat'], W,
                                                             b['equal'], sp, spg))
    figure(a['readings'], a['null'], D, sp, od / 'figure.png')
    return v


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', default='all', choices=('part0', 'A', 'all'))
    ap.add_argument('--n-null', type=int, default=N_NULL)
    ap.add_argument('--n-deals', type=int, default=N_DEALS)
    args = ap.parse_args(argv)
    t0 = time.time()
    od = out_dir()
    p0 = part0()
    print('part 0 reproduced E119 for every fly', flush=True)
    if args.part == 'part0':
        return
    a = part_a(args.n_null)
    a['readings'].to_csv(od / 'a_readings.csv', index=False)
    a['shares'].to_csv(od / 'a_shares.csv', index=False)
    a['null'].to_csv(od / 'a_null.csv', index=False)
    a['null_summary'].to_csv(od / 'a_null_summary.csv', index=False)
    a['repeat'].to_csv(od / 'a_repeat.csv', index=False)
    print('part A done', flush=True)
    if args.part == 'A':
        return
    b = part_b(args.n_deals, spread=flies_range())
    b['deals'].to_parquet(local_dir() / 'b_deals.parquet')
    summ = (b['deals'].groupby(['acq', 'wiring', 'variant'], sort=False)
            .agg(frac_pass=('passes', 'mean'), totals_cv_median=('totals_cv', 'median'),
                 totals_cv_p95=('totals_cv', lambda x: x.quantile(0.95)), resp_cv_median=('resp_cv', 'median'),
                 top1_median=('top1', 'median'), n_rois=('n_rois', 'first'), n_gloms_median=('n_gloms', 'median'))
            .reset_index())
    summ.to_csv(od / 'b_deals.csv', index=False)
    b['equal'].to_csv(od / 'b_equal.csv', index=False)
    b['leads'].to_csv(od / 'b_leads.csv', index=False)
    v = score(a, b)
    _write_json(od / 'provenance.json', dict(script=H.script_path(__file__), **H.code_commits(),
        generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), seconds=round(time.time() - t0, 1),
        settings=dict(CRITERION=CRITERION, HALF_WAY=HALF_WAY, LEAVE_OUT=LEAVE_OUT, N_NULL=args.n_null,
                      N_DEALS=args.n_deals, SEED=SEED, MIN_GROUP=MIN_GROUP, EQUAL_PER_GLOM=EQUAL_PER_GLOM,
                      EQUAL_DRAWS=EQUAL_DRAWS, NULL_FRAC=NULL_FRAC, MODEL_FRAC=MODEL_FRAC, FLY_MIN=FLY_MIN,
                      WIRING_MIN=WIRING_MIN, REPEAT_MIN=REPEAT_MIN),
        inputs=dict(remy_per_roi=str(H.shared_dir() / 'data' / 'remy_per_roi'),
                    e119=str(H.shared_dir() / 'results' / 'E119_stages_per_roi'),
                    e90_archive=str(_e119().archive_dir())),
        part0=dict(passed=p0['passed'])))
    print(json.dumps({k: v[k] for k in ('outcome', 'flies_class', 'repeat_check', 'model_class')}, default=str), flush=True)
    print(f'done in {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main()
