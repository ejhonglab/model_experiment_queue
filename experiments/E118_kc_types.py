#!/usr/bin/env python3
"""E118: the KC types in the best family and in the recordings: each type's response rate, the two ρ's and the per-pair
APL effect by type and by the sample's type composition (no simulations), then each type's response rate held at
measured ratios (parts B and C, on the entry's rule only).

The entry is item E118 of experiment_priorities.md (proposed by viki 2026-11-19 19:16; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-19 19:07). "The family" is the breadth-2 APL-effect fit (C2, E67) on its eight
wirings (the hemibrain wiring and seven re-drawn within type; E72's part B, `E77_input_side_common_error.members`), read
per wiring and as the population (the mean over wirings of the pair distances per condition, E72's). "The two ρ's" are
ρ on odor-odor distances with APL on against the control flies and with APL off against APL output blocked (TNT) (E106's
`two_rhos`; the KC matrices the mean over recordings of the cosine distance), the model's the cosine distance on raw
spike counts. The KC types are the hemibrain's labels of the model's KCs (αβ 'ab', α′β′ "a'b'", γ 'g').

Parts (CLI `--part`):
- 0, identities: the all-KC reading of each wiring's saved counts reproduces E106's pair distances (raw counts) and its
  population ρ's (0.690 / 0.756 over the 120 pairs); the types' KCs and spikes add up to each wiring's; the weighted
  reading with every weight 1 is the plain reading; E108's 'all' targets are E106's data matrices.
- A (no simulations): A1 by type (each type's KCs alone: response rates and their ratio, share of spikes, the two ρ's
  over the 120, 66, 54 and 91 pairs and against halves A and B of E77's split, the distance of paa from B-myr, the
  APL-effect statistics on raw counts and under E65's calibrated readout with K set on all KCs); A2 the composition (γ
  KCs weighted 0.5, 1, 2, 4, 8 and α′β′ KCs 0 or 1, as if the recorded sample held that share; a weight w multiplies a
  KC's row by sqrt(w), so that cosine similarities count it w times), also against Remy's validation2 somata on E107's
  as-fitted runs (a description); A3 controls (the same weights with the type labels permuted among KCs of equal claw
  count, 20 permutations; random subsets of each type's size; each type's distance of paa from B-myr against its
  response rate over the 48 type x condition x wiring points); A4 what is measured (published per-type response rates,
  only as on record in the queue's notes; E80's saved per-type points read for the two ρ's from E106's archive table);
  the rule for part B and the outcome letter as far as part A decides it.
- B, C: the per-type threshold factors (E80_screen's `type_factors` through E77's `hold_rates`) and the refit; stubs
  only until E114's runs are done.

Outputs: <shared>/results/E118_kc_types/; nothing large (part A reads archives only).
Usage (from the repository root):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E118_kc_types.py --part 0
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E118_kc_types.py --part A
"""
from __future__ import annotations

import argparse
import functools
import json
import socket
import time
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E118_kc_types'
TYPES = ('ab', "a'b'", 'g')
TYPE_LABEL = {'ab': 'αβ', "a'b'": 'α′β′', 'g': 'γ', 'all': 'all KCs'}
G_WEIGHTS = (0.5, 1.0, 2.0, 4.0, 8.0)
APBP_WEIGHTS = (0.0, 1.0)
N_PERM = 20
N_SUBSETS = 20
HALVES = ('A', 'B')
PAA_BMYR = ('paa', 'B-myr')
# the rule for part B (fixed before any reading is scored)
RULE_GAIN, RULE_OTHER, RULE_MAX_W, RULE_WIRINGS, RULE_CONTROL = 0.010, 0.005, 4.0, 6, 0.5


def out_dir() -> Path:
    """<shared>/results/E118_kc_types/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def members() -> List[Tuple[str, str, str]]:
    """The family's eight wirings, [(label, folder, tag)] (E77's `members`)."""
    from experiments import E77_input_side_common_error as E77
    return list(E77.members())


@functools.lru_cache(maxsize=None)
def _family() -> Tuple[Tuple[str, Dict[str, tuple]], ...]:
    from experiments import E72_wiring_population as E72
    out = []
    for lab, folder, tag in members():
        o = E72.load_fit(folder, tag)
        if o is None:
            raise RuntimeError(f'{folder}:{tag}: no saved outputs')
        out.append((lab, o))
    return tuple(out)


def load_family() -> Dict[str, Dict[str, tuple]]:
    """{wiring label: {'on', 'off'}: (responses, spike counts)} of the eight saved fits (E72's `load_fit`)."""
    return dict(_family())


WIRING_PATH = H.shared_dir() / 'data' / 'model_outputs' / 'inputs' / 'wPNKC_bouton_base.parquet'


@functools.lru_cache(maxsize=None)
def _claws_by_kc() -> pd.Series:
    W = pd.read_parquet(WIRING_PATH)
    return pd.Series(W.index.get_level_values(H.KC_ID_LEVEL)).value_counts()


def claw_counts(index: pd.Index) -> pd.Series:
    """Claws per KC (the hemibrain bouton base's wiring, archived inputs/wPNKC_bouton_base.parquet), in `index`'s
    order. The re-draws keep each KC's claw count, so one count serves all eight wirings."""
    c = _claws_by_kc()
    ids = index.get_level_values(H.KC_ID_LEVEL)
    out = pd.Series(c.reindex(ids).to_numpy(), index=index, name='claws')
    if out.isna().any():
        raise ValueError(f'{int(out.isna().sum())} KCs without claws in {WIRING_PATH.name}')
    return out.astype(int)


def kc_types(index: pd.Index) -> np.ndarray:
    """The KC type label of each row of `index` (level 'kc_type')."""
    return np.asarray(index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)


def kc_weights(labels: np.ndarray, g: float = 1.0, apbp: float = 1.0, ab: float = 1.0) -> np.ndarray:
    """Per-KC composition weights: `g` for the γ KCs, `apbp` for the α′β′ KCs, `ab` for the αβ KCs."""
    labels = np.asarray(labels, dtype=object)
    unknown = set(labels) - set(TYPES)
    if unknown:
        raise ValueError(f'unknown KC types {sorted(unknown)}')
    w = np.ones(len(labels), dtype=float)
    w[labels == 'g'] = g
    w[labels == "a'b'"] = apbp
    w[labels == 'ab'] = ab
    return w


def permuted_labels(labels: np.ndarray, claws: np.ndarray, seed: int) -> np.ndarray:
    """The type labels permuted among KCs of equal claw count (each claw count's labels shuffled among its KCs)."""
    labels = np.asarray(labels, dtype=object)
    claws = np.asarray(claws)
    rng = np.random.default_rng(seed)
    out = labels.copy()
    for c in np.unique(claws):
        i = np.flatnonzero(claws == c)
        out[i] = labels[i[rng.permutation(len(i))]]
    return out


def weighted_pair_distances(counts: pd.DataFrame, weights: np.ndarray, kind: str = 'raw', param=None) -> pd.Series:
    """Cosine odor-odor distances after E65's readout (`kind`, `param`) with each KC's row multiplied by sqrt(weight)
    (KCs of weight 0 dropped); weights all 1 give `kcq.harness.readout_pair_distances`."""
    from pyolfsysm import readout as R
    w = np.asarray(weights, dtype=float)
    if len(w) != len(counts):
        raise ValueError(f'{len(w)} weights for {len(counts)} KCs')
    if (w < 0).any():
        raise ValueError('negative weights')
    X = R.readout(H.abbrev_columns(counts), kind, param)
    keep = w > 0
    X = X.loc[keep]
    if not np.all(w[keep] == 1.0):
        X = X.mul(np.sqrt(w[keep]), axis=0)
    return R.pair_distances(X, measure='cosine', order=H.APL16_ODORS)


def weighted_rate(responses: pd.DataFrame, weights: np.ndarray) -> float:
    """The response rate of the weighted sample: sum over KCs of weight x fraction of odors responded to, over the
    summed weight."""
    w = np.asarray(weights, dtype=float)
    frac = (responses.to_numpy(dtype=float) > 0.5).mean(axis=1)
    return float((w * frac).sum() / w.sum())


def wiring_K(o: Mapping[str, tuple]) -> float:
    """E65's calibrated readout's K of one fit, set on all its KCs' APL-on counts (E72's `distances` 'cal')."""
    from experiments.E67_apl_effect_fits import refs
    return float(H.calibrate_readout(o['on'][1], 'saturating', refs()['cal_target'], stat='participation_ratio',
                                     responders_min=2.0))


def distances_for(o: Mapping[str, tuple], weights: np.ndarray, K: float) -> Dict[str, Dict[str, pd.Series]]:
    """{'raw', 'cal'}: {'on', 'off'} weighted pair distances of one fit."""
    return {'raw': {c: weighted_pair_distances(o[c][1], weights) for c in ('on', 'off')},
            'cal': {c: weighted_pair_distances(o[c][1], weights, 'saturating', K) for c in ('on', 'off')}}


@functools.lru_cache(maxsize=None)
def _targets() -> dict:
    from experiments import E106_rho_gap_parts as E106
    from experiments import E108_paa_extrapolation as E108
    dt, et = E108.dist_targets(), E108.effect_targets()
    D = E106.data_matrices()
    return {'dist': {'all': {c: H.canonical_pairs(D[c]) for c in ('on', 'off')},
                     **{h: {c: H.canonical_pairs(dt[h][c]) for c in ('on', 'off')} for h in HALVES}},
            'effect': {'all': H.yang_apl_effect(), **{h: et[h] for h in HALVES}}}


def targets() -> dict:
    """The KC data: 'dist' {'all', 'A', 'B'}: {'on', 'off'} mean matrices (all recordings E106's `data_matrices`, halves
    E108's `dist_targets`); 'effect' {'all', 'A', 'B'}: the per-pair APL effect (all recordings the harness's
    `yang_apl_effect`, halves E108's `effect_targets`)."""
    return _targets()


SUBSET_KW = {'120': {}, '66': dict(odors=H.MEGAMAT12), '54': None, '91': dict(without=PAA_BMYR)}


def _pairs54(index: pd.Index) -> List[Tuple[str, str]]:
    return [p for p in index if p[0] in H.VALIDATION2_4 or p[1] in H.VALIDATION2_4]


def _pair_value(s: pd.Series, a: str, b: str) -> float:
    s = H.canonical_pairs(s) if len(s) == 120 else s
    for k in ((a, b), (b, a)):
        if k in s.index:
            return float(s.loc[k])
    return float('nan')


def score_distances(d: Mapping[str, Mapping[str, pd.Series]], T: Mapping, full: bool = True) -> Dict[str, float]:
    """The readings of one set of distances: the two ρ's over the 120, 66, 54 and 91 pairs (all recordings) and over the
    120 against halves A and B; the distance of paa from B-myr per condition; the mean distance per condition; the
    APL-effect statistics (loss, mean, pattern, amplitude, shape, spread ratio, ρ and Pearson r on the per-pair APL
    effect) on raw counts and under the calibrated readout against all recordings, and the loss, mean and pattern part
    against each half. With `full` False, only the two ρ's over the 120 pairs against all recordings and each half."""
    from experiments import E106_rho_gap_parts as E106
    raw = d['raw']
    row: Dict[str, float] = {}
    for sub, kw in SUBSET_KW.items():
        if sub != '120' and not full:
            continue
        if kw is None:
            kw = dict(pairs=_pairs54(H.canonical_pairs(raw['on']).index))
        r = E106.two_rhos(raw['on'], raw['off'], T['dist']['all']['on'], T['dist']['all']['off'], **kw)
        row[f'rho_on_{sub}'], row[f'rho_off_{sub}'] = r['on'], r['off']
    for h in HALVES:
        r = E106.two_rhos(raw['on'], raw['off'], T['dist'][h]['on'], T['dist'][h]['off'])
        row[f'rho_on_120_half{h}'], row[f'rho_off_120_half{h}'] = r['on'], r['off']
    if not full:
        return row
    for c in ('on', 'off'):
        row[f'paa_bmyr_{c}'] = _pair_value(raw[c], *PAA_BMYR)
        row[f'mean_dist_{c}'] = float(np.nanmean(raw[c].to_numpy(dtype=float)))
    for rd in ('raw', 'cal'):
        on, off = d[rd]['on'], d[rd]['off']
        st = H.apl_effect_stats(on, off, T['effect']['all'])
        ss = H.apl_effect_spread_shape(on, off, T['effect']['all'])
        for k in ('rmse', 'mean', 'rmse_centered', 'spearman', 'pearson'):
            row[f'{rd}_{k}'] = float(st[f'apl_effect_{k}'])
        for k in ('amplitude', 'shape', 'sd_ratio'):
            row[f'{rd}_{k}'] = float(ss[f'apl_effect_{k}'])
        for h in HALVES:
            sh = H.apl_effect_stats(on, off, T['effect'][h])
            for k in ('rmse', 'mean', 'rmse_centered'):
                row[f'{rd}_{k}_half{h}'] = float(sh[f'apl_effect_{k}'])
    return row


def _mean_d(ds: Sequence[Mapping[str, Mapping[str, pd.Series]]]) -> Dict[str, Dict[str, pd.Series]]:
    from experiments import E106_rho_gap_parts as E106
    return {rd: E106.mean_pairs([x[rd] for x in ds]) for rd in ds[0]}


def _spike_share(counts: pd.DataFrame, weights: np.ndarray, base: np.ndarray) -> float:
    s = counts.to_numpy(dtype=float).sum(axis=1)
    return float((np.asarray(weights) * s).sum() / (np.asarray(base) * s).sum())


def reading(fam: Mapping[str, Mapping[str, tuple]], weights: Mapping[str, np.ndarray], Ks: Mapping[str, float],
            T: Mapping, per_wiring: str = 'full', base: Optional[Mapping[str, np.ndarray]] = None) -> pd.DataFrame:
    """Rows per wiring and 'population' (the mean over wirings of the distances): response rates on / off and their
    ratio, the share of spikes per condition (of the weighted sample over the `base` sample's, default all KCs), and
    `score_distances` (`per_wiring` 'full', 'rho' (the ρ's over the 120 only) or 'none' (the population only))."""
    rows, ds = [], []
    for lab, o in fam.items():
        w = weights[lab]
        b = np.ones(len(w)) if base is None else base[lab]
        d = distances_for(o, w, Ks[lab])
        ds.append(d)
        row = dict(wiring=lab, n_kcs=int((w > 0).sum()), weight_sum=float(w.sum()))
        for c in ('on', 'off'):
            row[f'rate_{c}'] = weighted_rate(o[c][0], w)
            row[f'spike_share_{c}'] = _spike_share(o[c][1], w, b)
        row['rate_ratio'] = row['rate_off'] / row['rate_on'] if row['rate_on'] > 0 else float('nan')
        if per_wiring != 'none':
            row.update(score_distances(d, T, full=per_wiring == 'full'))
        rows.append(row)
    per = pd.DataFrame(rows)
    pop = dict(wiring='population', n_kcs=float(per.n_kcs.mean()), weight_sum=float(per.weight_sum.mean()))
    for k in ('rate_on', 'rate_off', 'spike_share_on', 'spike_share_off'):
        pop[k] = float(per[k].mean())
    pop['rate_ratio'] = pop['rate_off'] / pop['rate_on'] if pop['rate_on'] > 0 else float('nan')
    pop.update(score_distances(_mean_d(ds), T))
    return pd.concat([per, pd.DataFrame([pop])], ignore_index=True)


def _labels(fam) -> np.ndarray:
    """The type labels (the same KCs, in the same order, in every wiring; checked)."""
    idx = [o['on'][1].index for o in fam.values()]
    for i in idx[1:]:
        if not i.equals(idx[0]):
            raise RuntimeError('the wirings have different KCs or KC orders')
    return kc_types(idx[0])


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, default=float) + '\n')


def part0(out: Path) -> dict:
    """Part 0's identities (part0_identity.json); raises on a failure."""
    from experiments import E106_rho_gap_parts as E106
    fam = load_family()
    T = targets()
    res: Dict[str, object] = dict(wirings={})
    worst = 0.0
    for (lab, folder, tag), o in zip(members(), fam.values()):
        saved = E106.load_pairs(folder, tag)
        r = {}
        labels = kc_types(o['on'][1].index)
        for c in ('on', 'off'):
            mine = H.canonical_pairs(weighted_pair_distances(o[c][1], np.ones(len(labels))))
            plain = H.canonical_pairs(H.readout_pair_distances(o[c][1], 'raw'))
            sv = H.canonical_pairs(saved[c])
            r[f'max_diff_saved_{c}'] = float((mine - sv).abs().max())
            r[f'max_diff_plain_{c}'] = float((mine - plain).abs().max())
            worst = max(worst, r[f'max_diff_saved_{c}'], r[f'max_diff_plain_{c}'])
            n = {t: int((labels == t).sum()) for t in TYPES}
            spikes = {t: float(o[c][1].to_numpy(dtype=float)[labels == t].sum()) for t in TYPES}
            r[f'n_by_type'] = n
            r[f'types_add_up_{c}'] = bool(sum(n.values()) == len(labels) and
                                          abs(sum(spikes.values()) - float(o[c][1].to_numpy(dtype=float).sum())) < 1e-6)
        res['wirings'][lab] = r
    pop = E106.population(members())
    rhos = E106.two_rhos(pop['on'], pop['off'])
    e106 = pd.read_csv(H.results_dir('E106_rho_gap_parts') / 'conditions.csv')
    e = e106[(e106.family == 'family') & (e106.reading == 'population') & (e106.subset == 'all120')].iloc[0]
    fam_d = [distances_for(o, np.ones(len(o['on'][1])), 1.0) for o in fam.values()]
    mine = score_distances(_mean_d(fam_d), T, full=False)
    res.update(rho_on_e106=float(e.rho_on), rho_off_e106=float(e.rho_off), rho_on_mine=mine['rho_on_120'],
               rho_off_mine=mine['rho_off_120'], rho_on_e106_population=rhos['on'], rho_off_e106_population=rhos['off'])
    dt = __import__('experiments.E108_paa_extrapolation', fromlist=['x']).dist_targets()
    res['e108_all_vs_e106_data_max_diff'] = float(max((H.canonical_pairs(dt['all'][c]) - T['dist']['all'][c]).abs().max()
                                                      for c in ('on', 'off')))
    ok = (worst < 1e-9 and abs(mine['rho_on_120'] - float(e.rho_on)) < 1e-9 and abs(mine['rho_off_120'] - float(e.rho_off)) < 1e-9
          and all(v['types_add_up_on'] and v['types_add_up_off'] for v in res['wirings'].values()))
    res.update(max_abs_diff=worst, passes=bool(ok), host=socket.gethostname(), time=time.strftime('%Y-%m-%d %H:%M'),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _write_json(out / 'part0_identity.json', res)
    print(f"part 0: max |difference| {worst:.2e}; ρ's {mine['rho_on_120']:.4f} / {mine['rho_off_120']:.4f} against E106's "
          f"{float(e.rho_on):.4f} / {float(e.rho_off):.4f}; passes {ok}", flush=True)
    if not ok:
        raise RuntimeError(f'part 0 fails: {res}')
    return res


def type_anatomy(fam: Mapping[str, Mapping[str, tuple]]) -> pd.DataFrame:
    """Per type and all KCs: the number of KCs, the mean claws per KC, and on the real wiring the mean APL→KC weight
    per KC (the connectome's, normalized, and the family's breadth-weighted ones), each relative to all KCs' mean."""
    idx = next(iter(fam.values()))['on'][1].index
    labels = kc_types(idx)
    claws = claw_counts(idx).to_numpy()
    rows = [dict(kc_type=t, n_kcs=int(m.sum()), share_kcs=float(m.mean()), claws_mean=float(claws[m].mean()))
            for t, m in [('all', np.ones(len(labels), bool))] + [(t, labels == t) for t in TYPES]]
    df = pd.DataFrame(rows)
    try:
        from experiments import E10_joint_fit as J
        from experiments import E77_input_screen as S77
        P = S77.c2_parts()
        tmap = dict(zip(idx.get_level_values(H.KC_ID_LEVEL), labels))
        for name, w in (('apl_kc_connectome_rel', P['wnorm']['wAPLKC']),
                        ('apl_kc_breadth_rel', J.breadth_weights(P['wnorm']['wAPLKC'], P['breadth'], S77.BREADTH_K))):
            w = pd.Series(np.asarray(w, dtype=float).ravel(), index=w.index)
            per_kc = w.groupby(level=H.KC_ID_LEVEL).sum()
            ty = per_kc.index.map(tmap)
            m = per_kc.mean()
            df[name] = [1.0] + [float(per_kc[np.asarray(ty == t)].mean() / m) for t in TYPES]
    except Exception as e:  # the weights are a description; their absence leaves the rest of part A
        print(f'APL→KC weights by type not read: {e!r}', flush=True)
    return df


def _type_weights(fam, typ: str, labels=None) -> Dict[str, np.ndarray]:
    labels = _labels(fam) if labels is None else labels
    w = np.ones(len(labels)) if typ == 'all' else (labels == typ).astype(float)
    return {lab: w for lab in fam}


def part_a1(out: Path, fam, T, Ks) -> pd.DataFrame:
    """A1 (by_type.csv): per type and all KCs, per wiring and the population."""
    parts = []
    for typ in ('all',) + TYPES:
        r = reading(fam, _type_weights(fam, typ), Ks, T)
        r.insert(0, 'kc_type', typ)
        parts.append(r)
    df = pd.concat(parts, ignore_index=True)
    df.to_csv(out / 'by_type.csv', index=False)
    return df


def _comp_weights(fam, labels, g, apbp) -> Dict[str, np.ndarray]:
    w = kc_weights(labels, g=g, apbp=apbp)
    return {lab: w for lab in fam}


def part_a2(out: Path, fam, T, Ks) -> pd.DataFrame:
    """A2 (composition.csv): per γ weight and α′β′ weight, per wiring and the population."""
    labels = _labels(fam)
    parts = []
    for apbp in APBP_WEIGHTS:
        for g in G_WEIGHTS:
            r = reading(fam, _comp_weights(fam, labels, g, apbp), Ks, T)
            r.insert(0, 'w_apbp', apbp)
            r.insert(0, 'w_g', g)
            parts.append(r)
    df = pd.concat(parts, ignore_index=True)
    df.to_csv(out / 'composition.csv', index=False)
    return df


def validation2_composition(out: Path) -> pd.DataFrame:
    """A2 against Remy's validation2 somata (E107's as-fitted runs, ρ over her carrying set of pairs and its pairs of paa
    or B-myr, APL on; a description; composition_validation2.csv). Empty when E107's runs or Remy's matrix cannot be
    read here."""
    try:
        from experiments import E106_rho_gap_parts as E106
        from experiments import E107_validation2_panel as E107
        runs = E107.family_runs('validation2', 'fitted')
        remy = E107.remy_v2('table_cos')
        pairs = E107.set_pairs(E107.carrying_set())
    except Exception as e:
        print(f"Remy's validation2 reading skipped: {e!r}", flush=True)
        return pd.DataFrame()
    from pyolfsysm import readout as R
    rows = []
    for apbp in APBP_WEIGHTS:
        for g in G_WEIGHTS:
            per = {}
            for w_lab, o in runs.items():
                s = H.abbrev_columns(o['on'][1])
                w = kc_weights(kc_types(s.index), g=g, apbp=apbp)
                X = s.astype(float).loc[w > 0].mul(np.sqrt(w[w > 0]), axis=0)
                per[w_lab] = H.canonical_pairs(R.pair_distances(X, measure='cosine'))
            pop = pd.concat(list(per.values()), axis=1).dropna().mean(axis=1)
            for sub in ('all', 'pb'):
                v, n = E107.rho(pop, remy, E107.subset_of(pairs, sub))
                rows.append(dict(w_g=g, w_apbp=apbp, remy='table_cos', set=E107.carrying_set(), subset=sub, rho_on=v,
                                 n_pairs=n, paa_bmyr_on=_pair_value(pop, *PAA_BMYR)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'composition_validation2.csv', index=False)
    return df


def _sparseness(a1: pd.DataFrame, subsets: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pts = []
    for _, r in a1[(a1.kc_type != 'all') & (a1.wiring != 'population')].iterrows():
        for c in ('on', 'off'):
            pts.append(dict(kc_type=r.kc_type, condition=c, wiring=r.wiring, rate=r[f'rate_{c}'],
                            paa_bmyr=r[f'paa_bmyr_{c}']))
    P = pd.DataFrame(pts)
    y = P.paa_bmyr.to_numpy(dtype=float)
    fits = []
    X0 = np.column_stack([np.ones(len(P)), P.rate.to_numpy(dtype=float)])
    Xt = np.column_stack([X0] + [(P.kc_type == t).to_numpy(dtype=float) for t in ("a'b'", 'g')])
    for name, X in (('rate', X0), ('rate + type', Xt)):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ b
        r2 = 1 - float((resid ** 2).sum() / ((y - y.mean()) ** 2).sum())
        row = dict(model=name, n_points=len(P), r2=r2, intercept=b[0], slope_rate=b[1])
        if len(b) > 2:
            row.update(type_apbp_vs_ab=b[2], type_g_vs_ab=b[3])
        fits.append(row)
    # the random subsets (no type) on the rate line of the 'rate' fit: their residuals' mean per condition
    if len(subsets):
        b0 = fits[0]
        for c in ('on', 'off'):
            s = subsets[subsets.wiring == 'population']
            res = s[f'paa_bmyr_{c}'] - (b0['intercept'] + b0['slope_rate'] * s[f'rate_{c}'])
            fits.append(dict(model=f'random subsets, residual from the rate line, APL {c}', n_points=len(s),
                             r2=np.nan, intercept=float(res.mean()), slope_rate=np.nan))
    return P, pd.DataFrame(fits)


def part_a3(out: Path, fam, T, Ks) -> Dict[str, pd.DataFrame]:
    """A3: the permuted-label control (composition_permuted.csv), random subsets (random_subsets.csv), the sparseness
    points and their fit (sparseness_points.csv, sparseness_fit.csv; needs by_type.csv)."""
    labels = _labels(fam)
    claws = claw_counts(next(iter(fam.values()))['on'][1].index).to_numpy()
    parts = []
    for seed in range(N_PERM):
        pl = permuted_labels(labels, claws, seed)
        for apbp in APBP_WEIGHTS:
            for g in G_WEIGHTS:
                if g == 1.0 and apbp == 1.0 and seed:
                    continue
                r = reading(fam, _comp_weights(fam, pl, g, apbp), Ks, T, per_wiring='rho')
                r.insert(0, 'w_apbp', apbp)
                r.insert(0, 'w_g', g)
                r.insert(0, 'perm_seed', seed)
                parts.append(r)
        print(f'A3 permutation {seed + 1} of {N_PERM}', flush=True)
    perm = pd.concat(parts, ignore_index=True)
    perm.to_csv(out / 'composition_permuted.csv', index=False)
    rng = np.random.default_rng(118)
    sub_parts = []
    for typ in TYPES:
        n = int((labels == typ).sum())
        for i in range(N_SUBSETS):
            w = np.zeros(len(labels))
            w[rng.choice(len(labels), n, replace=False)] = 1.0
            r = reading(fam, {lab: w for lab in fam}, Ks, T, per_wiring='none')
            r.insert(0, 'draw', i)
            r.insert(0, 'size_of', typ)
            sub_parts.append(r)
    subsets = pd.concat(sub_parts, ignore_index=True)
    subsets.to_csv(out / 'random_subsets.csv', index=False)
    a1 = pd.read_csv(out / 'by_type.csv')
    pts, fit = _sparseness(a1, subsets)
    pts.to_csv(out / 'sparseness_points.csv', index=False)
    fit.to_csv(out / 'sparseness_fit.csv', index=False)
    return dict(perm=perm, subsets=subsets, points=pts, fit=fit)


def published_type_rates() -> pd.DataFrame:
    """Published per-type response rates as on record in the queue's notes (type, rate, ratio to αβ, source); empty
    when none is on record. As of 2026-11-19 none is: E5's per-type settings (α′β′ twice as excitable, APL coupling
    α′β′ x1.6, γ x0.6) are model choices without a cited recording, and E80, E88 and notes/yang-summary.md give no
    per-type response rate of recorded KCs. The entry forbids looking further (nothing published is read)."""
    return pd.DataFrame(columns=['kc_type', 'rate', 'ratio_to_ab', 'source'])


E80_FOLDER = 'E80_remaining_candidates'


def part_a4(out: Path, a1: pd.DataFrame) -> pd.DataFrame:
    """A4 (measured.csv: the family's per-type response rates and ratios to αβ beside the published ones, 'not on hand'
    where none is on record; e80_type_points.csv: E80's per-type points (B4 a type's threshold with APL off, B5 a
    type's APL→KC weights) with the two ρ's from E106's archive table, against C2 on the real wiring)."""
    pub = published_type_rates()
    pop = a1[a1.wiring == 'population'].set_index('kc_type')
    rows = []
    for t in TYPES:
        p = pub[pub.kc_type == t]
        rows.append(dict(kc_type=t, family_rate_on=pop.loc[t, 'rate_on'], family_rate_off=pop.loc[t, 'rate_off'],
                         family_ratio_to_ab_on=pop.loc[t, 'rate_on'] / pop.loc['ab', 'rate_on'],
                         family_ratio_to_ab_off=pop.loc[t, 'rate_off'] / pop.loc['ab', 'rate_off'],
                         published_rate=p.rate.iloc[0] if len(p) else 'not on hand',
                         published_ratio_to_ab=p.ratio_to_ab.iloc[0] if len(p) else 'not on hand',
                         source=p.source.iloc[0] if len(p) else "none on record in the queue's notes (E5, E80, E88, "
                                                              "notes/yang-summary.md)"))
    m = pd.DataFrame(rows)
    m.to_csv(out / 'measured.csv', index=False)
    arch = pd.read_csv(H.results_dir('E106_rho_gap_parts') / 'archive.csv', low_memory=False)
    V = 'spike_counts_cosine'
    cols = {f'{V}__dist_rho_on_all120': 'rho_on_120', f'{V}__dist_rho_off_all120': 'rho_off_120',
            f'{V}__apl_effect_rmse': 'raw_rmse', f'{V}__apl_effect_rmse_centered': 'raw_rmse_centered',
            f'{V}__apl_effect_mean': 'raw_mean'}

    def find(folder, tag):
        hit = arch[(arch.folder == folder) & (arch.tag == tag)]
        if not len(hit):
            hit = arch[arch.copies.astype(str).str.contains(f'{folder}:{tag}', regex=False)]
        return hit
    base = find('E67_apl_effect_fits', 'rate0.068_aplboth_k2_cap3')
    e80 = arch[(arch.folder == E80_FOLDER) & arch.tag.str.contains(r'C2_B4_type|C2_B4_beta_perm|C2_B5_', regex=True)]
    pts = pd.concat([base.assign(point='C2 on the real wiring (base)'), e80.assign(point=e80.tag)], ignore_index=True)
    pts = pts[['point', 'folder', 'tag'] + list(cols)].rename(columns=cols)
    if len(base):
        for k in ('rho_on_120', 'rho_off_120'):
            pts[f'{k}_change'] = pts[k] - float(base[[c for c, v in cols.items() if v == k][0]].iloc[0])
    pts.to_csv(out / 'e80_type_points.csv', index=False)
    return m


def _gains(comp: pd.DataFrame, g: float, apbp: float, base: Optional[pd.DataFrame] = None) -> Dict[str, object]:
    """The changes of the two ρ's at weights (`g`, `apbp`) of `comp` from all KCs (`base`, default `comp`'s own row at
    weights 1; every labelling gives the same all-KC reading)."""
    base = comp if base is None else base
    base = base[(base.w_g == 1.0) & (base.w_apbp == 1.0)].set_index('wiring')
    cur = comp[(comp.w_g == g) & (comp.w_apbp == apbp)].set_index('wiring')
    wir = [w for w in base.index if w != 'population']
    out = {}
    for c in ('on', 'off'):
        out[f'gain_{c}_pop'] = float(cur.loc['population', f'rho_{c}_120'] - base.loc['population', f'rho_{c}_120'])
        for h in HALVES:
            k = f'rho_{c}_120_half{h}'
            out[f'gain_{c}_half{h}'] = float(cur.loc['population', k] - base.loc['population', k])
        if f'rho_{c}_120' in cur.columns and cur.loc[wir, f'rho_{c}_120'].notna().all():
            out[f'n_wirings_gain_{c}'] = int(((cur.loc[wir, f'rho_{c}_120'] - base.loc[wir, f'rho_{c}_120']) >= RULE_GAIN).sum())
            out[f'gain_{c}_wirings'] = (cur.loc[wir, f'rho_{c}_120'] - base.loc[wir, f'rho_{c}_120']).round(4).tolist()
    return out


def rule_b(comp: pd.DataFrame, perm: pd.DataFrame, a1: pd.DataFrame, published: pd.DataFrame) -> dict:
    """The rule for part B: 'outside' (γ's APL-on response rate over αβ's outside the published ratios; None when none
    is on hand) and 'composition' (at a γ weight of RULE_MAX_W or less, one ρ higher by RULE_GAIN or more on
    RULE_WIRINGS of 8 wirings and against each half, the other not lower by more than RULE_OTHER (population, all
    recordings and each half), and the permuted-label control (mean over permutations of its population gain against all
    recordings) gaining less than RULE_CONTROL times as much); 'run_b' if either holds. The settings compared are every γ
    weight up to RULE_MAX_W other than 1, with α′β′ at 0 or 1, and α′β′ at 0 with γ at 1, each against all KCs."""
    pop = a1[a1.wiring == 'population'].set_index('kc_type')
    ratio = float(pop.loc['g', 'rate_on'] / pop.loc['ab', 'rate_on'])
    pr = published[published.kc_type == 'g'].ratio_to_ab.astype(float) if len(published) else pd.Series(dtype=float)
    outside = None if not len(pr) else bool(ratio < pr.min() or ratio > pr.max())
    settings = []
    for apbp in APBP_WEIGHTS:
        for g in G_WEIGHTS:
            if g > RULE_MAX_W or (g == 1.0 and apbp == 1.0):
                continue
            gn = _gains(comp, g, apbp)
            pg = [_gains(perm[perm.perm_seed == s], g, apbp, comp) for s in sorted(perm.perm_seed.unique())]
            row = dict(w_g=g, w_apbp=apbp, **gn)
            passes = False
            for c, other in (('on', 'off'), ('off', 'on')):
                ctrl = float(np.mean([p[f'gain_{c}_pop'] for p in pg]))
                row[f'perm_gain_{c}_pop'] = ctrl
                ok = (gn.get(f'n_wirings_gain_{c}', 0) >= RULE_WIRINGS
                      and all(gn[f'gain_{c}_half{h}'] >= RULE_GAIN for h in HALVES)
                      and gn[f'gain_{other}_pop'] >= -RULE_OTHER
                      and all(gn[f'gain_{other}_half{h}'] >= -RULE_OTHER for h in HALVES)
                      and ctrl < RULE_CONTROL * gn[f'gain_{c}_pop'])
                row[f'passes_{c}'] = bool(ok)
                passes |= ok
            row['passes'] = bool(passes)
            settings.append(row)
    composition = any(s['passes'] for s in settings)
    return dict(gamma_over_ab_on=ratio, published_ratios=pr.tolist(), outside=outside, composition=bool(composition),
                run_b=bool(composition or bool(outside)), settings=settings)


def outcome_a(rule: dict) -> dict:
    """The outcome letter as far as part A decides it: (b) when the composition can move ρ; (a) when it cannot and the
    family's γ rate is inside what is measured; when nothing measured is on hand, (a) for the composition only, with
    the measured comparison open (parts B and C, and (c) to (e), wait for it or for the rule)."""
    if rule['composition']:
        return dict(letter='b', text='the composition can move ρ: the unknown type composition of the recorded somata is '
                                     'an uncertainty of every ρ (size in composition.csv)')
    if rule['outside'] is None:
        return dict(letter='a (composition only)',
                    text="the composition moves nothing by the rule; whether the family's type rates are inside what is "
                         "measured is open (no published per-type rate on hand), so part B does not run on the rule")
    if rule['outside']:
        return dict(letter='c or d (part B decides)', text="the family's γ rate is outside what is measured: part B runs")
    return dict(letter='a', text="the family's type rates are inside what is measured and the composition moves nothing")


def _fmt(x, spec='.3f') -> str:
    try:
        return format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def comparison_table(out: Path, a1: pd.DataFrame, comp: pd.DataFrame, rule: dict, outcome: dict) -> str:
    """comparison_table.md: by type and by composition (the population), against all KCs and the KC data."""
    pop1 = a1[a1.wiring == 'population'].set_index('kc_type')
    popc = comp[comp.wiring == 'population']
    data = dict(rate_on=H.data_target('noise_ctrl_frac'), rate_off=H.data_target('noise_tnt_frac'),
                rate_ratio=H.data_target('noise_rate_ratio_fraction'),
                paa_bmyr_on=_pair_value(targets()['dist']['all']['on'], *PAA_BMYR),
                paa_bmyr_off=_pair_value(targets()['dist']['all']['off'], *PAA_BMYR),
                raw_mean=float(targets()['effect']['all'].mean()))
    cols = [('rate_on', 'response rate, APL on'), ('rate_off', 'response rate, APL off'), ('rate_ratio', 'off / on'),
            ('spike_share_on', 'share of spikes, APL on'), ('rho_on_120', 'ρ on odor-odor distances, APL on (120)'),
            ('rho_off_120', 'ρ on odor-odor distances, APL off (120)'), ('rho_on_91', 'APL on (91, no paa / B-myr)'),
            ('rho_off_91', 'APL off (91)'), ('paa_bmyr_on', 'paa from B-myr, APL on'), ('paa_bmyr_off', 'paa from B-myr, APL off'),
            ('raw_rmse', 'APL-effect loss (raw)'), ('raw_mean', 'mean APL effect (raw)'), ('raw_rmse_centered', 'pattern part (raw)'),
            ('cal_rmse', 'APL-effect loss (calibrated)'), ('cal_rmse_centered', 'pattern part (calibrated)'),
            ('cal_amplitude', 'amplitude part (calibrated)'), ('cal_shape', 'shape part (calibrated)'),
            ('cal_sd_ratio', 'spread ratio (calibrated)'), ('cal_spearman', 'ρ on the per-pair APL effect (calibrated)'),
            ('cal_pearson', 'Pearson r on the per-pair APL effect (calibrated)')]
    heads = [TYPE_LABEL[t] for t in ('all',) + TYPES]
    comps = [(1.0, 1.0), (2.0, 1.0), (4.0, 1.0), (8.0, 1.0), (1.0, 0.0), (4.0, 0.0)]
    heads += [f'γ x{g:g}<br>α′β′ x{a:g}' for g, a in comps] + ['KC data']
    lines = ['# E118 part A: the KC types in the best family (population of eight wirings)', '',
             'Columns: each type alone, then the composition weights (γ KCs counted g times, α′β′ KCs a times), then the '
             'KC data (types not labelled). Raw counts unless "calibrated" (E65\'s saturating readout, K set on all KCs '
             'of each wiring).', '', '| statistic | ' + ' | '.join(heads) + ' |', '|---|' + '---:|' * len(heads)]
    for k, name in cols:
        vals = [_fmt(pop1.loc[t, k]) for t in ('all',) + TYPES]
        for g, a in comps:
            r = popc[(popc.w_g == g) & (popc.w_apbp == a)]
            vals.append(_fmt(r[k].iloc[0]) if len(r) and k != 'spike_share_on' else '')
        dv = data.get(k)
        vals.append(_fmt(dv) if dv is not None else '')
        lines.append(f'| {name} | ' + ' | '.join(vals) + ' |')
    lines += ['', f"Rule for part B: γ's APL-on response rate over αβ's {rule['gamma_over_ab_on']:.3f}; published ratios "
              f"{'not on hand' if rule['outside'] is None else rule['published_ratios']}; the composition moves ρ by the "
              f"rule: {rule['composition']}; part B runs on the rule: {rule['run_b']}.",
              f"Outcome (part A): ({outcome['letter']}) {outcome['text']}.", '',
              'Concentrations: paa is the 10^-4 recording for the flies\' 10^-3, euc and -aPine are at 10^-1.75 and '
              '10^-1.5 for the flies\' 10^-2; every 120-pair and 54-pair number and the distance of paa from B-myr '
              'involve them.']
    txt = '\n'.join(lines) + '\n'
    (out / 'comparison_table.md').write_text(txt)
    return txt


def realism_table(out: Path, a1: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    """realism_table.csv: each type's and each composition's response rates with APL on and off and their ratio beside
    the KC data's, with the bands (APL on 0.045 to 0.085, off / on 1.6 to 2.9; `kcq.apl_effect`) marked. The bands are
    for all-KC samples; the per-type rows are read against them as a description."""
    from kcq import apl_effect as AE
    rows = []
    for _, r in a1[a1.wiring == 'population'].iterrows():
        rows.append(dict(reading=f'type {TYPE_LABEL[r.kc_type]}', rate_on=r.rate_on, rate_off=r.rate_off, ratio=r.rate_ratio))
    for _, r in comp[comp.wiring == 'population'].iterrows():
        rows.append(dict(reading=f'composition γ x{r.w_g:g}, α′β′ x{r.w_apbp:g}', rate_on=r.rate_on, rate_off=r.rate_off,
                         ratio=r.rate_ratio))
    df = pd.DataFrame(rows)
    df['inside_rate_on_band'] = df.rate_on.between(*AE.RATE_ON_BAND)
    df['inside_ratio_band'] = df.ratio.between(*AE.RATE_RATIO_BAND)
    df = pd.concat([df, pd.DataFrame([dict(reading='KC data (noise units)', rate_on=H.data_target('noise_ctrl_frac'),
                                           rate_off=H.data_target('noise_tnt_frac'),
                                           ratio=H.data_target('noise_rate_ratio_fraction'))])], ignore_index=True)
    df.to_csv(out / 'realism_table.csv', index=False)
    return df


def part_a(out: Path) -> dict:
    """Part A: A1 to A4, the rule, the outcome, comparison_table.md, realism_table.csv, provenance_A.json."""
    t0 = time.time()
    fam = load_family()
    T = targets()
    Ks = {lab: wiring_K(o) for lab, o in fam.items()}
    anat = type_anatomy(fam)
    anat.to_csv(out / 'type_anatomy.csv', index=False)
    a1 = part_a1(out, fam, T, Ks)
    print('A1 done', flush=True)
    comp = part_a2(out, fam, T, Ks)
    print('A2 done', flush=True)
    v2 = validation2_composition(out)
    a3 = part_a3(out, fam, T, Ks)
    print('A3 done', flush=True)
    pub = published_type_rates()
    part_a4(out, a1)
    rule = rule_b(comp, a3['perm'], a1, pub)
    outcome = outcome_a(rule)
    _write_json(out / 'rule_b.json', dict(rule=rule, outcome=outcome))
    txt = comparison_table(out, a1, comp, rule, outcome)
    realism_table(out, a1, comp)
    _write_json(out / 'provenance_A.json', dict(
        script='experiments/E118_kc_types.py', part='A', host=socket.gethostname(), time=time.strftime('%Y-%m-%d %H:%M'),
        seconds=time.time() - t0, members=members(), K=Ks, g_weights=G_WEIGHTS, apbp_weights=APBP_WEIGHTS, n_perm=N_PERM,
        n_subsets=N_SUBSETS, rule=dict(gain=RULE_GAIN, other=RULE_OTHER, max_w=RULE_MAX_W, wirings=RULE_WIRINGS,
                                       control=RULE_CONTROL), validation2_rows=int(len(v2)), **H.code_commits()))
    print(txt, flush=True)
    return dict(rule=rule, outcome=outcome)


def type_threshold_factors(index: pd.Index, ratios: Mapping[str, float]) -> pd.Series:
    """Part B: per-KC threshold factors (one per type) aiming the types' APL-on response rates at `ratios` to αβ's."""
    raise NotImplementedError('part B waits for E114\'s runs')


def part_b(out: Path) -> None:
    """Part B (only on the rule): held readings at measured, equal and half-way ratios, with controls."""
    raise NotImplementedError('part B waits for E114\'s runs')


def part_c(out: Path) -> None:
    """Part C (only on 'better' or 'an equal' in part B): the family fitted again with the types' ratios held."""
    raise NotImplementedError('part C waits for part B')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C'))
    args = ap.parse_args()
    out = out_dir()
    if args.part == '0':
        part0(out)
    elif args.part == 'A':
        if not (out / 'part0_identity.json').exists() or not json.loads((out / 'part0_identity.json').read_text())['passes']:
            part0(out)
        part_a(out)
    elif args.part == 'B':
        part_b(out)
    else:
        part_c(out)


if __name__ == '__main__':
    main()
