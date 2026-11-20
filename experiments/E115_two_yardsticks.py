#!/usr/bin/env python3
"""E115: the megamat gap against two sets of flies and the yardstick of the price in response rate ρ (no simulations).

The entry is item E115 of experiment_priorities.md (proposed by viki 2026-11-19 16:58; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 16:47; moved into the list and claimed by viki 17:00). "The family" is the breadth-2
APL-effect fit (C2, E67) on its eight wirings (E72's part B; S = real, P1, P2, M1; T = P3, P4, M2, M3), read as the mean
matrix over the wirings (E106's `population`); "the two ρ's" are ρ on odor-odor distances with APL on against the control
flies and with APL off against APL output blocked (TNT) (E106's `data_matrices`: the mean over Yang's recordings of the
cosine distance). "Remy's somata" are her wild-type megamat KC somata, read per fly: as correlation distance with E90's
`remy_soma_fly_distances` (18 flies) and as cosine distance with E107's `megamat_somata('cosine')` under E107's lock
(from the copy in data/remy_validation2_kc/, E109's `remy_copy`). The same Remy somata stand beside both conditions
(her flies are wild type).

Parts (CLI `--part`):
- 0, the identities (part0.json; parts A and B stop if one fails): the family's two ρ's over the 120 and the 66 megamat
  pairs reproduce E106's conditions.csv; Remy's matrices are read through E107's lock (`require_lock`) unchanged; the
  response rate ρ read from the model-output archive's response matrices reproduces E106's realism.csv for the family
  and E109's runs.csv for its reference and readings on T.
- A, the two sets of flies (66 megamat pairs, per condition): every reader (the family, each wiring, the ORN input read
  four ways, E106's `input_matrices`, of which input_pearson / input_cosine are the family's own input) against Yang's
  flies, against Remy's somata and against the mean rank of the two sets (`two_set_rank`); the correlation of the
  reader's rank residuals against Yang's flies with Remy's flies' rank residuals against Yang's (`residual_correlation`);
  1000 joint resamples of Yang's recordings, Remy's flies and the eight wirings (`resample_part_a`); the pairs by class
  (`pair_classes`) with E61's shares of the gap per class and per odor (`class_shares`); the validation2 odors against
  Remy's validation2 somata (`validation2_table`, E107's runs of the family on that panel); the saved candidates read the
  same way (`candidates_table`). Writes a_*.csv.
- B, the yardstick of response rate ρ: every input change at fixed scales with both response counts re-held on the
  family's wirings and no claim behind it (`CONTROLS`), each paired with its base on the same wiring; per wiring the
  change of response rate ρ with APL on and off, the mean over T's four and over the eight wirings, the same with the
  wirings' per-odor response fractions pooled before ρ is taken, the paired SD, and the size of the input change
  (`input_size`); E109's six inputs and E108's X3 placed in the controls' distribution (`place`). Writes b_*.csv.
- verdict: the three verdicts and the outcome letters (verdict.json, comparison_table.md), the registry rows
  (registry.csv, read by `registry_entries`), provenance.json.

The verdicts (the entry's, fixed before any reading was scored):
- "Shared": the residual correlation of the family is SHARED_MIN (0.30) or more in both conditions and its 5th percentile
  over the resamples is above 0.
- "The megamat level is the flies'": the family's ρ against the two sets' mean rank is within LEVEL_TOL (0.03) of Remy's
  flies' ρ against Yang's, in both conditions.
- "The price is inside the controls": E109's fall in response rate ρ with APL on, in each direction of its split, is
  inside the 5th to 95th percentile of the control changes of comparable size, and the pooled reading falls by no more
  than POOLED_MAX_FALL (0.02).
No class of E97, E108 or E109 is rewritten; no correlation chooses anything; no pair is weighted; the fruity x
small-ester block is a column only.

Readings chosen where the entry is open (fixed here before any reading was scored):
- The residual correlation is Spearman's ρ between the two sets of rank residuals over the 66 pairs (E61's residuals,
  E106's `rank_residuals`; ranks over the pairs both hold), the statistic the planner's probes read (+0.43 / +0.37);
  Pearson's r beside. The verdict reads Remy's somata as correlation distance (E90's 18 flies; her published reading),
  the cosine reading beside.
- Outcome (e): an odor "carries" the class "Yang's flies alone" when it is in at least half of the class's pairs; (e)
  holds when one or two odors carry a class of three or more pairs, in either condition.
- Beside the verdict (no class): the residual correlation with the reader's residuals against one random half of Yang's
  recordings and Remy's residuals against the other half (`residual_correlation_halves`), which removes the share of the
  correlation that Yang's own measurement noise puts into both residuals.
- "The two sets' mean rank": per pair, the mean of its rank among Yang's 66 distances and its rank among Remy's.
- Classes (`pair_classes`): ranks over the 66 pairs; "both sets against the family" when |rank Yang - rank Remy| <=
  CLASS_NEAR (10) and the family's rank is more than CLASS_FAR (15) from both; "Yang's flies alone" when |rank Remy - rank
  family| <= 10 and Yang's rank is more than 15 from both; otherwise "neither". Per condition.
- The size of an input change (`input_size`): the Frobenius norm of the change of the model's ORN input (spike deltas,
  the 16 odors x the input's glomeruli) over the norm of the lab's input. Changes that are not changes of the ORN input
  (E97's claw proportions, E102's antennal-lobe steps) have no size: they enter the distribution of all controls, not
  that of comparable size.
- "Comparable size": a control whose size is within COMPARABLE (0.5 to 2 times) the mean size of the direction's two
  searched inputs. If fewer than MIN_COMPARABLE (10) controls are comparable, the verdict reads against every control
  that has a size and says so.
- A direction's fall is E109's (`direction_values`' d_rate_rho_on: the mean over its two searches' runs on T minus the
  mean over the family's runs on T); a control's change is the same for one input (its runs on T minus the family's on T).
  The pooled reading of a direction pools the per-odor response fractions over its eight runs and the family's over its
  four before ρ is taken.
- E109's control searches (ctlA, ctlB) enter at both sizes they were read at (0.6 and 1); E108's W5 (paa at 10^-5) is a
  wrong-way control, not a permuted or scaled one, and is given beside only; E108's and E97's ceiling-family runs are not
  read.

Concentrations. paa's input is the 10^-4 recording for Yang's flies' 10^-3 (Remy's flies got 10^-4); euc and -aPine are at
10^-1.75 and 10^-1.5 for her flies' 10^-2. Every number of the validation2 odors (step 5) and every 120-pair number of part
B (the response rate ρ is over the 16 odors) involves them; the 66 megamat pairs do not.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): the hemibrain wiring and its re-draws only (E97's male-CNS
permutations, Pratyush's tables, are among the controls read from saved runs); no table of Yang's own Claude's project, no
FlyWire / FAFB / BANC data, no raw-only recording. Nothing is published to claude.ai.

Outputs: <shared>/results/E115_two_yardsticks/ (comparison_table.md, verdict.json, a_*.csv, b_*.csv, registry.csv,
part0.json, provenance.json).
Usage (from the repository root; no simulations):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E115_two_yardsticks.py --part all
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import socket
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats as sstats

from kcq import harness as H
from kcq import apl_effect as AE

EXPERIMENT = 'E115_two_yardsticks'
RESULTS = H.shared_dir() / 'results'
CONDITIONS = ('on', 'off')
GENOTYPE = {'on': 'TNT3c', 'off': 'TNT'}
DATA_LABEL = {'on': 'control flies', 'off': 'APL output blocked (TNT)'}
SHARED_MIN = 0.30
LEVEL_TOL = 0.03
POOLED_MAX_FALL = 0.02
CLASS_NEAR, CLASS_FAR = 10, 15
CARRY_SHARE, CARRY_MIN_PAIRS = 0.5, 3
COMPARABLE = (0.5, 2.0)
MIN_COMPARABLE = 10
N_BOOT, SEED, N_HALVES = 1000, 0, 200
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
ALL_WIRINGS = S_WIRINGS + T_WIRINGS
CLASSES = ('both sets against the family', "Yang's flies alone", 'neither')
REMY_KINDS = ('corr', 'cos')
REMY_LABEL = {'corr': "Remy's somata, correlation distance (18 flies; E90's reader)",
              'cos': "Remy's somata, cosine distance (18 flies; E107's reader, under its lock)"}
FOUR = ('paa', 'B-myr', 'euc', '-aPine')
BLOCK = (('2h', 'IaA', 'pa'), ('2-but', 'eb', 'ep'))   # the fruity x small-ester block (a column only)
E109_DIRECTIONS = {'A': ('A1', 'A2'), 'B': ('B1', 'B2'), 'second': ('C1', 'D1')}
E109_MAIN = ('A1', 'A2', 'B1', 'B2', 'C1', 'D1')
E102_REFIT = {'real': 'rate0.068_aplboth_pnadd1000sc544.444_k2_cap3',
              'P1': 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_pnadd1000sc544.444_k2_cap3',
              'P2': 'wiring-claw_permutation_within_type-seed2_rate0.068_aplboth_pnadd1000sc544.444_k2_cap3',
              'M1': 'wiring-curveball_within_type_matched-seed1_rate0.068_aplboth_pnadd1000sc544.444_k2_cap3'}
F93, F97, F102, F108, F109 = ('E93_missing_glomeruli_input', 'E97_glomerulus_claw_shares',
                              'E102_antennal_lobe_normalization', 'E108_paa_extrapolation', 'E109_input_search_rho')


def _ctl(source, change, kind, folder, tag, base='family', wirings=ALL_WIRINGS, input_file=None, base_file=None):
    return dict(source=source, change=change, kind=kind, folder=folder, tag=tag, base=base, wirings=tuple(wirings),
                input_file=input_file, base_file=base_file)


def _controls() -> List[dict]:
    """The entry's controls (step 7): every input change at fixed scales with both response counts re-held on the
    family's wirings and no claim behind it."""
    i93, i108, i109 = RESULTS / F93 / 'inputs', RESULTS / F108 / 'inputs', RESULTS / F109 / 'inputs'
    out = [_ctl('E93', 'F', "fuller input (E93's F)", F93, 'B_F_{w}', input_file=i93 / 'F.csv', base_file=i93 / 'I0.csv')]
    for lev in (100, 314):
        out.append(_ctl('E102', f'inhadd{lev}', f'antennal-lobe step at pn_inhadd {lev}', F102, f'B_inhadd{lev}_{{w}}_held',
                        base='B_identity_{w}_held'))
    for ch, kind in (('permmean2_s1', 'claw proportions permuted'), ('permmean2_s2', 'claw proportions permuted'),
                     ('invmean2', 'claw proportions inverted'), ('perm_mcnsL_s1', 'claw proportions permuted'),
                     ('perm_mcnsR_s1', 'claw proportions permuted')):
        out.append(_ctl('E97', ch, kind, F97, f'B3_{ch}_{{w}}'))
    for ch in ('P1', 'P2', 'P3', 'P4'):
        out.append(_ctl('E108', ch, "paa's change permuted", F108, f'B_{ch}_{{w}}', input_file=i108 / f'{ch}.csv',
                        base_file=i108 / 'I0.csv'))
    out.append(_ctl('E108', 'S1', 'paa scaled', F108, 'B_S1_{w}', input_file=i108 / 'S1.csv', base_file=i108 / 'I0.csv'))
    for n in E109_MAIN:
        for s in (1, 2, 3, 4):
            out.append(_ctl('E109', f'{n}_perm{s}', 'searched factors permuted', F109, f'{n}_perm{s}_{{w}}', wirings=T_WIRINGS,
                            input_file=i109 / f'{n}_perm{s}.csv', base_file=i108 / 'I0.csv'))
        out.append(_ctl('E109', f'{n}_inv', 'searched factors inverted', F109, f'{n}_inv_{{w}}', wirings=T_WIRINGS,
                        input_file=i109 / f'{n}_inv.csv', base_file=i108 / 'I0.csv'))
    for n in ('ctlA', 'ctlB'):
        for sz in ('0.6', '1'):
            out.append(_ctl('E109', f'{n}_s{sz}', 'control search (2h, 6al)', F109, f'{n}_s{sz}_{{w}}', wirings=T_WIRINGS,
                            input_file=i109 / f'{n}_s{sz}.csv', base_file=i108 / 'I0.csv'))
    return out


def _candidates() -> List[dict]:
    """E109's six searched inputs at 0.6 and E108's X3 (placed in the controls' distribution); E108's W5 beside."""
    i108, i109 = RESULTS / F108 / 'inputs', RESULTS / F109 / 'inputs'
    out = [_ctl('E109', f'{n}_s0.6', 'searched input (E109)', F109, f'{n}_s0.6_{{w}}', wirings=T_WIRINGS,
                input_file=i109 / f'{n}_s0.6.csv', base_file=i108 / 'I0.csv') for n in E109_MAIN]
    out.append(_ctl('E108', 'X3', "paa's extrapolation X3 (E108)", F108, 'B_X3_{w}', input_file=i108 / 'X3.csv',
                    base_file=i108 / 'I0.csv'))
    out.append(_ctl('E108', 'W5', 'paa at 10^-5 (wrong-way control, beside)', F108, 'B_W5_{w}', input_file=i108 / 'W5.csv',
                    base_file=i108 / 'I0.csv'))
    return out


CONTROLS = _controls()
CANDIDATES = _candidates()

_INPUTS: Dict[str, str] = {}


def _note(p) -> Path:
    p = Path(p)
    if p.exists() and str(p) not in _INPUTS:
        _INPUTS[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _e106():
    from experiments import E106_rho_gap_parts as E106
    return E106


# ---------------------------------------------------------------------------
# the two sets of flies: the statistics (made-up matrices in the tests)

def _ranks(x: np.ndarray) -> np.ndarray:
    return sstats.rankdata(np.asarray(x, dtype=float))


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    return float(sstats.spearmanr(a, b)[0])


def _aligned(*xs: pd.Series) -> List[pd.Series]:
    return _e106()._aligned(*xs)


def row_arrays(m: np.ndarray, y: np.ndarray, r: np.ndarray) -> Dict[str, float]:
    """`reader_row` on aligned arrays (one value per pair)."""
    rm, ry, rr = _ranks(m), _ranks(y), _ranks(r)
    res_m, res_r = rm - ry, rr - ry
    two = (ry + rr) / 2
    return dict(rho_yang=_rho(m, y), rho_remy=_rho(m, r), rho_two_sets=_rho(m, two), remy_vs_yang=_rho(r, y),
                resid_corr=_rho(res_m, res_r), resid_corr_pearson=float(np.corrcoef(res_m, res_r)[0, 1]),
                n_pairs=int(len(m)))


def residual_correlation(model: pd.Series, yang: pd.Series, remy: pd.Series, method: str = 'spearman') -> float:
    """The correlation (Spearman's ρ, or 'pearson') over the pairs all three hold between the model's rank residuals
    against `yang` and `remy`'s rank residuals against `yang` (E106's `rank_residuals`: rank in one minus rank in the other,
    ranks over those pairs)."""
    m, y, r = _aligned(model, yang, remy)
    v = row_arrays(m.to_numpy(float), y.to_numpy(float), r.to_numpy(float))
    return v['resid_corr'] if method == 'spearman' else v['resid_corr_pearson']


def two_set_rank(yang: pd.Series, remy: pd.Series) -> pd.Series:
    """Per pair (the pairs both hold), the mean of its rank among `yang`'s distances and its rank among `remy`'s."""
    y, r = _aligned(yang, remy)
    return pd.Series((_ranks(y) + _ranks(r)) / 2, index=y.index, name='two_set_rank')


def pair_classes(model: pd.Series, yang: pd.Series, remy: pd.Series, near: int = CLASS_NEAR, far: int = CLASS_FAR) -> pd.DataFrame:
    """Per pair (the pairs all three hold): the three ranks and the class (CLASSES): 'both sets against the family' when
    |rank Yang - rank Remy| <= near and the family's rank is more than `far` from both; "Yang's flies alone" when |rank
    Remy - rank family| <= near and Yang's rank is more than `far` from both; else 'neither'."""
    m, y, r = _aligned(model, yang, remy)
    rm, ry, rr = _ranks(m), _ranks(y), _ranks(r)
    both = (np.abs(ry - rr) <= near) & (np.abs(rm - ry) > far) & (np.abs(rm - rr) > far)
    alone = (np.abs(rr - rm) <= near) & (np.abs(ry - rm) > far) & (np.abs(ry - rr) > far)
    cls = np.where(both, CLASSES[0], np.where(alone, CLASSES[1], CLASSES[2]))
    return pd.DataFrame(dict(rank_model=rm, rank_yang=ry, rank_remy=rr, cls=cls), index=m.index)


def class_shares(model: pd.Series, yang: pd.Series, classes: pd.Series) -> pd.DataFrame:
    """E61's contributions to 1 - ρ of `model` against `yang` (over the pairs `classes` holds), summed per class (level
    'class', odor '') and per odor within each class (level 'odor'; a pair counts for both its odors): n_pairs,
    contribution, share_of_gap (contribution over 1 - ρ)."""
    E106 = _e106()
    m = E106._canon(model).reindex(classes.index)
    rr = E106.rank_residuals(m, yang)
    cl = classes.reindex(rr.index)
    gap = 1 - H.spearman(rr['model_distance'], rr['data_distance'])
    rows = []
    for c in CLASSES:
        sel = rr[cl == c]
        rows.append(dict(level='class', cls=c, odor='', n_pairs=len(sel), contribution=float(sel['contribution'].sum()),
                         share_of_gap=float(sel['contribution'].sum() / gap) if gap > 0 else np.nan))
        ods = sorted({o for p in sel.index for o in p}, key=lambda o: list(H.APL16_ODORS).index(o) if o in H.APL16_ODORS else 99)
        for o in ods:
            s2 = sel[[o in p for p in sel.index]]
            rows.append(dict(level='odor', cls=c, odor=o, n_pairs=len(s2), contribution=float(s2['contribution'].sum()),
                             share_of_gap=float(s2['contribution'].sum() / gap) if gap > 0 else np.nan))
    return pd.DataFrame(rows)


def carriers(classes: pd.DataFrame, cls: str = CLASSES[1], share: float = CARRY_SHARE) -> List[str]:
    """The odors in at least `share` of the pairs of class `cls` (none when it has fewer than CARRY_MIN_PAIRS pairs)."""
    sel = classes.index[classes['cls'] == cls]
    if len(sel) < CARRY_MIN_PAIRS:
        return []
    cnt = pd.Series([o for p in sel for o in p]).value_counts()
    return [o for o, n in cnt.items() if n / len(sel) >= share]


def reader_row(model: pd.Series, yang: pd.Series, remy: pd.Series) -> Dict[str, float]:
    """One reader's part-A numbers over the pairs all three hold: ρ against Yang's (rho_yang), against Remy's (rho_remy),
    against the two sets' mean rank (rho_two_sets), the residual correlation (resid_corr; Pearson's beside), and Remy's
    flies against Yang's (remy_vs_yang)."""
    m, y, r = _aligned(model, yang, remy)
    return row_arrays(m.to_numpy(float), y.to_numpy(float), r.to_numpy(float))


def residual_correlation_halves(model: pd.Series, yang_rows: pd.DataFrame, remy: pd.Series, n_splits: int = N_HALVES,
                                seed: int = SEED) -> Tuple[float, float, float]:
    """Beside the verdict: the residual correlation with the model's residuals against the mean of one random half of
    Yang's recordings (`yang_rows`, recordings x pairs) and Remy's against the other half's, so that Yang's measurement
    noise does not enter both; (mean, 5th, 95th percentile) over `n_splits` splits."""
    rng = np.random.default_rng(seed)
    E106 = _e106()
    m = E106._canon(model)
    vals = []
    for _ in range(n_splits):
        a, b = E106._split_means(yang_rows, rng)
        mm, ya, yb, rr = _aligned(m, a, b, remy)
        res_m = _ranks(mm) - _ranks(ya)
        res_r = _ranks(rr) - _ranks(yb)
        vals.append(_rho(res_m, res_r))
    v = np.asarray(vals)
    return float(v.mean()), float(np.percentile(v, 5)), float(np.percentile(v, 95))


# ---------------------------------------------------------------------------
# the data and the readers

def megamat_index() -> pd.MultiIndex:
    return pd.MultiIndex.from_tuples(H.pair_list(list(H.MEGAMAT12)), names=['odor_a', 'odor_b'])


def _mm(s: pd.Series) -> pd.Series:
    """A pair Series on the 66 megamat pairs, canonical keys and order (NaN where it lacks one)."""
    return _e106()._canon(s).reindex(megamat_index())


@functools.lru_cache(maxsize=None)
def yang_rows() -> Dict[str, pd.DataFrame]:
    """{'on', 'off'}: Yang's per-recording cosine distances (recordings x pairs) of the control flies and of APL output
    blocked (TNT) (E106's reader)."""
    r = _e106()._recs()
    return {c: r[GENOTYPE[c]] for c in CONDITIONS}


def yang_mean() -> Dict[str, pd.Series]:
    return _e106().data_matrices()


@functools.lru_cache(maxsize=None)
def remy_rows(kind: str) -> pd.DataFrame:
    """Remy's wild-type megamat somata per fly on the 66 pairs (flies x pairs, NaN where a fly lacks an odor): 'corr'
    (E90's `remy_soma_fly_distances`, correlation distance), 'cos' (E107's `megamat_somata('cosine')` under E107's lock,
    from the copy in data/remy_validation2_kc/)."""
    if kind == 'corr':
        from experiments import E90_model_stages as E90
        from pyolfsysm.mb_model import data_root
        _note(Path(data_root) / 'from_remy' / '2024-11-12' / 'xrda_stim_rdm_concat__iscell_good_xid0__correlation.nc')
        rows = E90.remy_soma_fly_distances()
    elif kind == 'cos':
        from experiments import E109_input_search_rho as E109
        with E109.remy_copy() as E107:
            rows, info = E107.megamat_somata('cosine')
            _note(E107.remy_dir() / 'megamat17' / 'kc_soma' / 'xrds_stim_rdm_concat_mm.nc')
        if info['metric'] != 'cosine':
            raise RuntimeError(f"Remy's megamat file gave {info['metric']}, not cosine")
    else:
        raise ValueError(kind)
    E106 = _e106()
    out = pd.DataFrame({f: _mm(E106._canon(rows.loc[f].dropna())) for f in rows.index}).T
    return out


def remy_mean(kind: str) -> pd.Series:
    return remy_rows(kind).mean()


@functools.lru_cache(maxsize=None)
def family() -> Tuple[Dict[str, Dict[str, pd.Series]], Dict[str, pd.Series]]:
    """({wiring: {'on', 'off'}}, the population) of the family (E106's `_readings`)."""
    return _e106()._readings('family')


def readers() -> Dict[str, Dict[str, pd.Series]]:
    """{name: {'on', 'off'}}: the family, each wiring, and the ORN input read four ways (E106's `input_matrices`; the same
    matrix for both conditions)."""
    per, pop = family()
    out = {'family': pop}
    out.update({f'wiring {w}': per[w] for w in ALL_WIRINGS})
    for k, v in _e106().input_matrices().items():
        out[f'input {k}'] = {'on': v, 'off': v}
    return out


def candidate_readers() -> Dict[str, dict]:
    """Step 6: {name: dict(pairs={'on', 'off'}, base={'on', 'off'}, wirings)}: E109's six inputs at 0.6 on T, E108's X3 and
    E97's mean proportions on the eight, E102's refit at pn_inhadd 1000 on S; each with the family on the same wirings
    as its base."""
    E106 = _e106()
    per, _ = family()
    fam = lambda ws: E106.mean_pairs([per[w] for w in ws])
    out = {}
    for n in E109_MAIN:
        out[f'E109 {n} at 0.6'] = dict(pairs=E106.mean_pairs([E106.load_pairs(F109, f'{n}_s0.6_{w}') for w in T_WIRINGS]),
                                       base=fam(T_WIRINGS), wirings='T')
    out['E108 X3'] = dict(pairs=E106.mean_pairs([E106.load_pairs(F108, f'B_X3_{w}') for w in ALL_WIRINGS]),
                          base=fam(ALL_WIRINGS), wirings='eight')
    out["E97 mean proportions (B2)"] = dict(pairs=E106.mean_pairs([E106.load_pairs(F97, f'B2_mean2_{w}') for w in ALL_WIRINGS]),
                                            base=fam(ALL_WIRINGS), wirings='eight')
    out['E102 refit at pn_inhadd 1000'] = dict(pairs=E106.mean_pairs([E106.load_pairs(F102, E102_REFIT[w]) for w in S_WIRINGS]),
                                               base=fam(S_WIRINGS), wirings='S')
    return out


# ---------------------------------------------------------------------------
# part A

def part_a_table() -> pd.DataFrame:
    """Every reader and candidate, per condition and Remy reading: `reader_row` on the 66 pairs."""
    Y = yang_mean()
    rows = []
    rd = readers()
    cand = candidate_readers()
    for kind in REMY_KINDS:
        R = _mm(remy_mean(kind))
        for c in CONDITIONS:
            y = _mm(Y[c])
            for name, m in rd.items():
                rows.append(dict(reader=name, group='reader', remy=kind, condition=c, **reader_row(_mm(m[c]), y, R)))
            for name, d in cand.items():
                rows.append(dict(reader=name, group='candidate', wirings=d['wirings'], remy=kind, condition=c,
                                 **reader_row(_mm(d['pairs'][c]), y, R)))
                rows.append(dict(reader=f'{name}: its base (the family on the same wirings)', group='candidate base',
                                 wirings=d['wirings'], remy=kind, condition=c, **reader_row(_mm(d['base'][c]), y, R)))
    df = pd.DataFrame(rows)
    df['level_minus_flies'] = df['rho_two_sets'] - df['remy_vs_yang']
    return df


def _pair_matrix(series: Sequence[pd.Series]) -> np.ndarray:
    return np.vstack([_mm(s).to_numpy(float) for s in series])


def resample_part_a(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """The joint resamples (step 3): per resample, Yang's recordings (with replacement within each genotype), Remy's flies
    (with replacement, per reading) and the eight wirings (with replacement; the family's population the mean over the
    drawn wirings); each reader's `reader_row` per condition and Remy reading (the wirings read alone and the inputs see
    only the data resampled). One row per resample, reader, condition and Remy reading."""
    rng = np.random.default_rng(seed)
    recs = _e106()._recs()
    per, _ = family()
    W = {c: _pair_matrix([per[w][c] for w in ALL_WIRINGS]) for c in CONDITIONS}
    rd = {k: v for k, v in readers().items() if k != 'family'}
    fixed = {c: {k: _mm(v[c]).to_numpy(float) for k, v in rd.items()} for c in CONDITIONS}
    R = {k: remy_rows(k).to_numpy(float) for k in REMY_KINDS}
    rows = []
    for b in range(n_boot):
        ym = H.resample_recording_means(recs, rng)
        y = {c: _mm(ym[GENOTYPE[c]]).to_numpy(float) for c in CONDITIONS}
        rm = {}
        for k in REMY_KINDS:
            idx = rng.integers(0, len(R[k]), len(R[k]))
            with np.errstate(invalid='ignore'):
                rm[k] = np.nanmean(R[k][idx], axis=0)
        wi = rng.integers(0, len(ALL_WIRINGS), len(ALL_WIRINGS))
        for c in CONDITIONS:
            fam = W[c][wi].mean(axis=0)
            for k in REMY_KINDS:
                ok = ~np.isnan(rm[k])
                rows.append(dict(boot=b, reader='family', remy=k, condition=c, **row_arrays(fam[ok], y[c][ok], rm[k][ok])))
                for name, m in fixed[c].items():
                    rows.append(dict(boot=b, reader=name, remy=k, condition=c, **row_arrays(m[ok], y[c][ok], rm[k][ok])))
    df = pd.DataFrame(rows)
    df['level_minus_flies'] = df['rho_two_sets'] - df['remy_vs_yang']
    return df


BOOT_COLS = ('rho_yang', 'rho_remy', 'rho_two_sets', 'remy_vs_yang', 'resid_corr', 'resid_corr_pearson', 'level_minus_flies')


def with_intervals(point: pd.DataFrame, boot: pd.DataFrame) -> pd.DataFrame:
    """`point` with each BOOT_COLS's 5th and 95th percentile over the resamples (columns <col>_lo, <col>_hi)."""
    q = boot.groupby(['reader', 'remy', 'condition'])[list(BOOT_COLS)].quantile([0.05, 0.95]).unstack()
    q.columns = [f'{c}_{"lo" if p == 0.05 else "hi"}' for c, p in q.columns]
    return point.merge(q.reset_index(), on=['reader', 'remy', 'condition'], how='left')


def halves_table() -> pd.DataFrame:
    """`residual_correlation_halves` for the family and the inputs, per condition and Remy reading."""
    rows = []
    rd = readers()
    for kind in REMY_KINDS:
        R = _mm(remy_mean(kind))
        for c in CONDITIONS:
            for name in ['family'] + [k for k in rd if k.startswith('input')]:
                v, lo, hi = residual_correlation_halves(_mm(rd[name][c]), yang_rows()[c], R)
                rows.append(dict(reader=name, remy=kind, condition=c, resid_corr_halves=v, lo=lo, hi=hi))
    return pd.DataFrame(rows)


def classes_tables() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Step 4 on the family (Remy's correlation reading, the verdict's; the cosine reading beside): per pair the ranks,
    the class, the rank residual and the contribution (a_classes.csv, with 'block' marking the fruity x small-ester
    pairs); the shares per class and per odor (a_class_shares.csv)."""
    E106 = _e106()
    _, pop = family()
    Y = yang_mean()
    pr, sh = [], []
    for kind in REMY_KINDS:
        R = _mm(remy_mean(kind))
        for c in CONDITIONS:
            y, m = _mm(Y[c]), _mm(pop[c])
            cl = pair_classes(m, y, R)
            rr = E106.rank_residuals(m, y)
            t = cl.join(rr[['rank_resid', 'contribution']])
            t['block'] = [(a in BLOCK[0] and b in BLOCK[1]) or (a in BLOCK[1] and b in BLOCK[0]) for a, b in t.index]
            t = t.reset_index().assign(remy=kind, condition=c)
            pr.append(t)
            s = class_shares(m, y, cl['cls'])
            sh.append(s.assign(remy=kind, condition=c))
    return pd.concat(pr, ignore_index=True), pd.concat(sh, ignore_index=True)


def validation2_table() -> pd.DataFrame:
    """Step 5: the family's population on Remy's validation2 panel (E107's runs as fitted on the eight wirings) against
    her validation2 somata (her mean table, correlation and cosine; and the fly mean of the lock's per-fly reading) over
    every pair of the carrying set with paa, B-myr, euc or -aPine ('four'), with paa or B-myr ('pb'), and each odor's
    pairs; and, on the six pairs among the four odors (the only ones Yang's flies also hold), `reader_row` against Yang's
    flies and Remy's validation2 somata. Resampled over Remy's validation2 flies, Yang's recordings and the wirings (1000,
    seed 0) for the flymean readings."""
    from experiments import E107_validation2_panel as E107m
    from experiments import E109_input_search_rho as E109
    E106 = _e106()
    ref = E109.e107_reference(ALL_WIRINGS)
    pop = ref['pop']
    rows = []
    rng = np.random.default_rng(SEED)
    recs = E106._recs()
    with E109.remy_copy() as E107:
        key = E107.carrying_set()
        sp = E107.set_pairs(key)
        sets = {'four': E107.subset_of(sp, 'four'), 'pb': E107.subset_of(sp, 'pb')}
        for o in FOUR:
            sets[f'odor {o}'] = [p for p in sp if o in p]
        six = [p for p in sp if p[0] in FOUR and p[1] in FOUR]
        mats = {'table_corr': E107.remy_v2('table_corr'), 'table_cos': E107.remy_v2('table_cos'),
                'flymean_corr': E107.remy_v2('flymean_corr'), 'flymean_cos': E107.remy_v2('flymean_cos')}
        fly_rows = {'flymean_corr': E107.per_fly('correlation'), 'flymean_cos': E107.per_fly('cosine')}
        min_flies = E107.MIN_FLIES
    for c in CONDITIONS:
        m = H.canonical_pairs(pop[c])
        for rk, R in mats.items():
            for sname, prs in sets.items():
                r, n = E107m.rho(m, R, prs)
                rows.append(dict(part='validation2 pairs', condition=c, remy=rk, pairs=sname, n_pairs=n, rho_remy=r))
            y = E106._canon(yang_mean()[c])
            mm, yy, rr = _aligned(m.loc[E106.select_pairs(m.index, six)], y, R)
            rows.append(dict(part='six pairs of the four odors', condition=c, remy=rk, pairs='six', **row_arrays(
                mm.to_numpy(float), yy.to_numpy(float), rr.to_numpy(float))))
    # resampling of the fly-mean readings (Remy's validation2 flies, Yang's recordings, the wirings)
    per = {w: ref[w] for w in ALL_WIRINGS}
    boot = []
    for b in range(N_BOOT):
        wi = rng.integers(0, len(ALL_WIRINGS), len(ALL_WIRINGS))
        popb = E106.mean_pairs([per[ALL_WIRINGS[i]] for i in wi])
        ym = H.resample_recording_means(recs, rng)
        for rk, FR in fly_rows.items():
            idx = rng.integers(0, len(FR), len(FR))
            Rb = E107m.fly_mean(FR.iloc[idx].reset_index(drop=True), min_flies)
            for c in CONDITIONS:
                m = H.canonical_pairs(popb[c])
                r, n = E107m.rho(m, Rb, sets['four'])
                boot.append(dict(boot=b, condition=c, remy=rk, pairs='four', rho_remy=r))
                mm, yy, rr = _aligned(m.loc[E106.select_pairs(m.index, six)], E106._canon(ym[GENOTYPE[c]]), Rb)
                v = row_arrays(mm.to_numpy(float), yy.to_numpy(float), rr.to_numpy(float))
                boot.append(dict(boot=b, condition=c, remy=rk, pairs='six', **v))
    bt = pd.DataFrame(boot)
    df = pd.DataFrame(rows)
    q = bt.groupby(['condition', 'remy', 'pairs']).quantile([0.05, 0.95]).drop(columns='boot').unstack()
    q.columns = [f'{a}_{"lo" if p == 0.05 else "hi"}' for a, p in q.columns]
    q = q.dropna(axis=1, how='all').reset_index()
    return df.merge(q, on=['condition', 'remy', 'pairs'], how='left')


def validation2_shares() -> pd.DataFrame:
    """Step 5's E61 shares by odor of the family's validation2 population against Remy's validation2 somata (her mean
    table, correlation and cosine) over the carrying set's pairs with any of the four odors."""
    from experiments import E109_input_search_rho as E109
    E106 = _e106()
    pop = E109.e107_reference(ALL_WIRINGS)['pop']
    out = []
    with E109.remy_copy() as E107:
        sp = E107.set_pairs(E107.carrying_set())
        four = E107.subset_of(sp, 'four')
        for rk in ('table_corr', 'table_cos'):
            R = E107.remy_v2(rk)
            for c in CONDITIONS:
                g = E106.gap_shares_by_odor(H.canonical_pairs(pop[c]), R, four)
                out.append(g.reset_index().assign(condition=c, remy=rk))
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------------------
# part B: the yardstick of response rate ρ

def input_size(change: pd.DataFrame, base: pd.DataFrame) -> float:
    """||change - base||_F / ||base||_F over the union of both frames' glomeruli and the base's odors (a glomerulus
    one frame lacks counts as 0 there)."""
    gl = change.index.union(base.index)
    od = list(base.columns)
    c = change.reindex(index=gl, columns=od).fillna(0.0).to_numpy(float)
    b = base.reindex(index=gl, columns=od).fillna(0.0).to_numpy(float)
    return float(np.linalg.norm(c - b) / np.linalg.norm(b))


def _read_input(p: Path) -> pd.DataFrame:
    return pd.read_csv(_note(p), index_col=0)


@functools.lru_cache(maxsize=None)
def fractions(folder: str, tag: str) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the per-odor response fractions (two or more spikes, the archive's response matrices; the 16 odors
    in APL16 order) of a saved run, read as `kcq.apl_effect.archive_opk` reads them."""
    out = {}
    for c, name in (('on', 'on'), ('off', 'acute_off')):
        p = AE.archive_responses(RESULTS / folder, tag, name)
        if p is None:
            raise FileNotFoundError(f'no archived responses for {folder}:{tag} ({name})')
        r = pd.read_parquet(_note(p))
        out[c] = H.abbrev_columns(r.astype(float) > 0.5).mean().reindex(H.APL16_ODORS).astype(float)
    return out


def rate_rho(frac: pd.Series, cond: str) -> float:
    """Response rate ρ: Spearman of per-odor response fractions with the flies' in baseline-noise units (control flies
    for 'on', APL output blocked (TNT) for 'off'; `kcq.apl_effect.frac_data`)."""
    return H.spearman(frac.reindex(H.APL16_ODORS), AE.frac_data(cond, 'noise'))


def pooled_rate_rho(fracs: Sequence[pd.Series], cond: str) -> float:
    """Response rate ρ of the mean of several runs' per-odor response fractions."""
    return rate_rho(pd.concat(list(fracs), axis=1).mean(axis=1), cond)


def percentile_of(x: float, values: Sequence[float]) -> float:
    """The percentage of `values` below `x` (ties counted half)."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if not len(v):
        return np.nan
    return float(100 * ((v < x).sum() + 0.5 * (v == x).sum()) / len(v))


def _family_tag(w: str) -> Tuple[str, str]:
    return {lab: (f, t) for lab, f, t in _e106().family_members('family')}[w]


def _base_fracs(spec: Mapping, w: str) -> Dict[str, pd.Series]:
    if spec['base'] == 'family':
        return fractions(*_family_tag(w))
    return fractions(spec['folder'], spec['base'].format(w=w))


def change_rows(specs: Sequence[Mapping]) -> pd.DataFrame:
    """Per spec and wiring: the run's and its base's response rate ρ with APL on and off and their changes."""
    rows = []
    for s in specs:
        for w in s['wirings']:
            f, f0 = fractions(s['folder'], s['tag'].format(w=w)), _base_fracs(s, w)
            r = dict(source=s['source'], change=s['change'], kind=s['kind'], wiring=w, tag=s['tag'].format(w=w))
            for c in CONDITIONS:
                r[f'rate_rho_{c}'], r[f'base_rate_rho_{c}'] = rate_rho(f[c], c), rate_rho(f0[c], c)
                r[f'd_rate_rho_{c}'] = r[f'rate_rho_{c}'] - r[f'base_rate_rho_{c}']
            rows.append(r)
    return pd.DataFrame(rows)


def change_summary(specs: Sequence[Mapping], rows: pd.DataFrame) -> pd.DataFrame:
    """Per spec: the size of the input change; per condition the mean change over T's four wirings and over the eight
    (where read), the paired SD over its wirings, and the pooled reading (each set's per-odor response fractions averaged
    over the wirings before ρ is taken, minus the base's pooled the same way) over T and over the eight."""
    out = []
    for s in specs:
        g = rows[(rows.source == s['source']) & (rows.change == s['change'])]
        r = dict(source=s['source'], change=s['change'], kind=s['kind'], n_wirings=len(g))
        if s['input_file'] is not None:
            r['size'] = input_size(_read_input(s['input_file']), _read_input(s['base_file']))
        else:
            r['size'] = np.nan
        for c in CONDITIONS:
            for unit, ws in (('T', T_WIRINGS), ('8', ALL_WIRINGS)):
                gg = g[g.wiring.isin(ws)]
                ok = len(gg) == len(ws)
                r[f'd_{c}_{unit}'] = float(gg[f'd_rate_rho_{c}'].mean()) if ok else np.nan
                if ok:
                    pooled = pooled_rate_rho([fractions(s['folder'], s['tag'].format(w=w))[c] for w in ws], c)
                    pooled0 = pooled_rate_rho([_base_fracs(s, w)[c] for w in ws], c)
                    r[f'pooled_d_{c}_{unit}'] = pooled - pooled0
                else:
                    r[f'pooled_d_{c}_{unit}'] = np.nan
            r[f'paired_sd_{c}'] = float(g[f'd_rate_rho_{c}'].std(ddof=1))
        out.append(r)
    return pd.DataFrame(out)


def direction_values(cand_rows: pd.DataFrame, cand_summary: pd.DataFrame) -> Dict[str, dict]:
    """E109's directions (A: A1, A2; B: B1, B2; second split beside: C1, D1): per condition the fall in response rate ρ
    (the mean over the direction's runs on T minus the family's mean on T, E109's `direction_values`), the pooled fall (the
    per-odor response fractions pooled over the direction's eight runs, and the family's over its four), and the mean size
    of the two inputs."""
    out = {}
    for d, names in E109_DIRECTIONS.items():
        chs = [f'{n}_s0.6' for n in names]
        g = cand_rows[(cand_rows.source == 'E109') & cand_rows.change.isin(chs)]
        r = dict(changes=chs, size=float(cand_summary[cand_summary.change.isin(chs)]['size'].mean()))
        for c in CONDITIONS:
            r[f'd_rate_rho_{c}'] = float(g[f'rate_rho_{c}'].mean() - np.mean([rate_rho(fractions(*_family_tag(w))[c], c)
                                                                              for w in T_WIRINGS]))
            runs = [fractions(F109, f'{ch}_{w}')[c] for ch in chs for w in T_WIRINGS]
            base = [fractions(*_family_tag(w))[c] for w in T_WIRINGS]
            r[f'pooled_d_rate_rho_{c}'] = pooled_rate_rho(runs, c) - pooled_rate_rho(base, c)
        out[d] = r
    return out


def price_verdict(directions: Mapping[str, Mapping[str, float]], controls: pd.DataFrame,
                  verdict_directions: Sequence[str] = ('A', 'B')) -> dict:
    """"The price is inside the controls": for each direction of `verdict_directions`, its fall in response rate ρ with
    APL on (d_rate_rho_on) inside the 5th to 95th percentile of the controls' changes (d_on_T) among the controls of
    comparable size (size within COMPARABLE times the direction's size; every control with a size when fewer than
    MIN_COMPARABLE are), and its pooled fall (pooled_d_rate_rho_on) no more than POOLED_MAX_FALL."""
    out = dict(directions={})
    ok_all = True
    sized = controls[controls['size'].notna() & controls['d_on_T'].notna()]
    for d in verdict_directions:
        v = directions[d]
        lo_s, hi_s = COMPARABLE[0] * v['size'], COMPARABLE[1] * v['size']
        comp = sized[(sized['size'] >= lo_s) & (sized['size'] <= hi_s)]
        fallback = len(comp) < MIN_COMPARABLE
        use = sized if fallback else comp
        p5, p95 = (float(np.percentile(use['d_on_T'], q)) for q in (5, 95))
        inside = bool(p5 <= v['d_rate_rho_on'] <= p95)
        pooled_ok = bool(v['pooled_d_rate_rho_on'] >= -POOLED_MAX_FALL)
        out['directions'][d] = dict(fall=v['d_rate_rho_on'], size=v['size'], size_window=[lo_s, hi_s], n_comparable=len(comp),
                                    fallback_all_sized=fallback, n_used=len(use), p5=p5, p95=p95, inside=inside,
                                    percentile=percentile_of(v['d_rate_rho_on'], use['d_on_T']),
                                    pooled_fall=v['pooled_d_rate_rho_on'], pooled_ok=pooled_ok)
        ok_all &= inside and pooled_ok
    out['inside'] = bool(ok_all)
    return out


def place(cands: pd.DataFrame, controls: pd.DataFrame) -> pd.DataFrame:
    """Each candidate's change (mean over T, and over the eight where read) placed among the controls' (percentile, and
    the controls' 5th and 95th percentiles), with APL on and off."""
    rows = []
    for _, r in cands.iterrows():
        for c in CONDITIONS:
            for unit in ('T', '8'):
                x = r[f'd_{c}_{unit}']
                vals = controls[f'd_{c}_{unit}'].dropna()
                if np.isnan(x) or not len(vals):
                    continue
                rows.append(dict(source=r['source'], change=r['change'], condition=c, unit=unit, change_value=x, size=r['size'],
                                 n_controls=len(vals), percentile=percentile_of(x, vals),
                                 controls_p5=float(np.percentile(vals, 5)), controls_p95=float(np.percentile(vals, 95)),
                                 controls_mean=float(vals.mean()), controls_sd=float(vals.std(ddof=1))))
    return pd.DataFrame(rows)


def yardstick_summary(ctl_rows: pd.DataFrame, ctl: pd.DataFrame) -> pd.DataFrame:
    """The yardstick per condition and set of controls (all; each source): the SD, mean, 5th and 95th percentiles over
    the controls of the change in response rate ρ over T (and the eight, and pooled), and the paired SD over the
    controls' wiring-level changes after each control's own mean is removed."""
    rows = []
    groups = [('all', ctl)] + [(s, g) for s, g in ctl.groupby('source', sort=False)]
    for name, g in groups:
        rr = ctl_rows[ctl_rows.source.isin(g.source.unique()) & ctl_rows.change.isin(g.change)]
        for c in CONDITIONS:
            r = dict(controls=name, condition=c, n_controls=len(g))
            for col in (f'd_{c}_T', f'd_{c}_8', f'pooled_d_{c}_T', f'pooled_d_{c}_8'):
                v = g[col].dropna()
                key = col.replace(f'_{c}', '')
                r[f'{key}_n'] = len(v)
                r[f'{key}_mean'] = float(v.mean()) if len(v) else np.nan
                r[f'{key}_sd'] = float(v.std(ddof=1)) if len(v) > 1 else np.nan
                r[f'{key}_p5'] = float(np.percentile(v, 5)) if len(v) else np.nan
                r[f'{key}_p95'] = float(np.percentile(v, 95)) if len(v) else np.nan
            dev = rr[f'd_rate_rho_{c}'] - rr.groupby(['source', 'change'])[f'd_rate_rho_{c}'].transform('mean')
            dof = len(rr) - rr.groupby(['source', 'change']).ngroups
            r['paired_sd_within'] = float(np.sqrt((dev ** 2).sum() / dof)) if dof > 0 else np.nan
            r['wiring_change_sd'] = float(rr[f'd_rate_rho_{c}'].std(ddof=1))
            rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the verdicts

def verdicts(a: pd.DataFrame, classes: pd.DataFrame, price: Mapping) -> dict:
    """The three verdicts and the outcome letters from part A's table with intervals (`with_intervals`; the family's rows
    under Remy's correlation reading), part A's classes and `price_verdict`'s result."""
    fam = a[(a.reader == 'family') & (a.remy == 'corr')].set_index('condition')
    shared = {c: dict(resid_corr=float(fam.loc[c, 'resid_corr']), lo=float(fam.loc[c, 'resid_corr_lo']),
                      hi=float(fam.loc[c, 'resid_corr_hi'])) for c in CONDITIONS}
    is_shared = all(v['resid_corr'] >= SHARED_MIN and v['lo'] > 0 for v in shared.values())
    level = {c: dict(family_two_sets=float(fam.loc[c, 'rho_two_sets']), remy_vs_yang=float(fam.loc[c, 'remy_vs_yang']),
                     difference=float(fam.loc[c, 'rho_two_sets'] - fam.loc[c, 'remy_vs_yang'])) for c in CONDITIONS}
    is_level = all(abs(v['difference']) <= LEVEL_TOL for v in level.values())
    car = {}
    for c in CONDITIONS:
        t = classes[(classes.remy == 'corr') & (classes.condition == c)].set_index(['odor_a', 'odor_b'])
        car[c] = dict(n_alone=int((t['cls'] == CLASSES[1]).sum()), n_both=int((t['cls'] == CLASSES[0]).sum()),
                      carriers=carriers(t))
    e = any(0 < len(v['carriers']) <= 2 for v in car.values())
    letters = []
    if is_shared and is_level:
        letters.append('a')
    if not is_shared:
        letters.append('b')
    letters.append('c' if price['inside'] else 'd')
    if e:
        letters.append('e')
    return dict(shared=dict(holds=bool(is_shared), rule=f'residual correlation >= {SHARED_MIN} in both conditions and its 5th '
                            'percentile over the resamples above 0', values=shared),
                megamat_level=dict(holds=bool(is_level), rule=f'|family vs the two sets - Remy vs Yang| <= {LEVEL_TOL} in both '
                                   'conditions', values=level),
                price=dict(holds=bool(price['inside']), **price),
                yang_alone=dict(holds=bool(e), rule=f'one or two odors in at least {CARRY_SHARE:.0%} of the pairs of the class '
                                f"\"Yang's flies alone\" ({CARRY_MIN_PAIRS} or more pairs), either condition", values=car),
                outcomes=letters,
                outcome_note=('' if (is_shared and is_level) or not is_shared else
                              "shared, but the megamat level is not the flies': neither (a) nor (b)"))


# ---------------------------------------------------------------------------
# part 0

def part0() -> dict:
    """The identities (part0.json); raises if one fails."""
    E106 = _e106()
    res = {}
    # 1. the family's two ρ's over the 120 and the 66 pairs against E106's conditions.csv
    cond = pd.read_csv(_note(RESULTS / E106.EXPERIMENT / 'conditions.csv'))
    _, pop = family()
    D = yang_mean()
    diffs = []
    for s in ('all120', 'megamat66'):
        r = cond[(cond.family == 'family') & (cond.reading == 'population') & (cond.subset == s)].iloc[0]
        for c in CONDITIONS:
            diffs.append(abs(E106.rho(pop[c], D[c], s) - float(r[f'rho_{c}'])))
    res['family_rhos_vs_e106'] = dict(max_abs_diff=float(max(diffs)), ok=bool(max(diffs) < 1e-12))
    # 2. Remy's matrices through E107's lock
    from experiments import E109_input_search_rho as E109
    with E109.remy_copy() as E107:
        lock = E107.require_lock()
        res['e107_lock'] = dict(ok=True, sha256=lock.get('sha256'), files_checked=True)
    rc, rs = remy_rows('corr'), remy_rows('cos')
    res['remy_rows'] = dict(corr_flies=int(len(rc)), cos_flies=int(len(rs)), corr_pairs=int(rc.shape[1]),
                            cos_pairs=int(rs.shape[1]), ok=bool(rc.shape[1] == 66 and rs.shape[1] == 66 and len(rc) == 18))
    v = {c: H.spearman(_mm(remy_mean('corr')), _mm(D[c])) for c in CONDITIONS}
    reg = {'on': H.data_target('rho_ctrl_remy_soma'), 'off': H.data_target('rho_tnt_remy_soma')}
    res['remy_vs_yang_vs_registry'] = dict(ours=v, registry=reg, max_abs_diff=float(max(abs(v[c] - reg[c]) for c in CONDITIONS)),
                                           note='beside: the registry reads kc_remy_combo from shared12_pairs.csv')
    # 3. response rate ρ from the archive against E106's realism.csv and E109's runs.csv
    real = pd.read_csv(_note(RESULTS / E106.EXPERIMENT / 'realism.csv'))
    real = real[(real.family == 'family') & real.folder.notna()]
    d3 = []
    for _, r in real.iterrows():
        f = fractions(r.folder, r.tag)
        d3 += [abs(rate_rho(f[c], c) - float(r[f'rate_rho_{c}'])) for c in CONDITIONS]
    res['rate_rho_family_vs_e106'] = dict(n=len(d3), max_abs_diff=float(max(d3)), ok=bool(max(d3) < 1e-12))
    runs = pd.read_csv(_note(RESULTS / F109 / 'runs.csv'), low_memory=False)
    d4 = []
    for _, r in runs.iterrows():
        if r.change == 'I0':
            f = fractions(*_family_tag(r.wiring))
        elif r.change == 'E108ref':
            continue
        else:
            f = fractions(F109, r.tag)
        d4 += [abs(rate_rho(f[c], c) - float(r[f'rate_rho_{c}'])) for c in CONDITIONS]
    res['rate_rho_e109_runs'] = dict(n=len(d4), max_abs_diff=float(max(d4)), ok=bool(max(d4) < 1e-12))
    # 4. E109's direction values of response rate ρ against its gate.json
    g = json.loads(_note(RESULTS / F109 / 'gate.json').read_text())['gates']['six:I0:s0.6']
    crow = change_rows([c for c in CANDIDATES if c['source'] == 'E109'])
    csum = change_summary([c for c in CANDIDATES if c['source'] == 'E109'], crow)
    dv = direction_values(crow, csum)
    d5 = [abs(dv[d][f'd_rate_rho_{c}'] - g['directions'][d][f'd_rate_rho_{c}']) for d in ('A', 'B') for c in CONDITIONS]
    d5 += [abs(dv['second'][f'd_rate_rho_{c}'] - g['second'][f'd_rate_rho_{c}']) for c in CONDITIONS]
    res['e109_directions_vs_gate'] = dict(max_abs_diff=float(max(d5)), ok=bool(max(d5) < 1e-12))
    # 5. the lab's input as E93 and E108 saved it
    i93, i108 = _read_input(RESULTS / F93 / 'inputs' / 'I0.csv'), _read_input(RESULTS / F108 / 'inputs' / 'I0.csv')
    res['lab_input_e93_e108'] = dict(max_abs_diff=float(np.abs(i93.to_numpy(float) - i108.reindex_like(i93).to_numpy(float)).max()),
                                     ok=bool(np.allclose(i93.to_numpy(float), i108.reindex_like(i93).to_numpy(float), rtol=0, atol=1e-9)))
    res['all_ok'] = bool(all(v.get('ok', True) for v in res.values() if isinstance(v, dict)))
    (out_dir() / 'part0.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    if not res['all_ok']:
        raise RuntimeError(f'part 0 failed: {json.dumps(res, default=float)}')
    return res


def realism_table() -> pd.DataFrame:
    """The family's realism quantities per wiring and their mean (E106's realism.csv, which reads the archive), with the
    bands (`kcq.apl_effect`) and the flies' values beside: the APL-on response rate (band 0.045 to 0.085; control flies
    0.068), the APL-off / APL-on ratio (1.6 to 2.9), response rate ρ with APL on (0.37 or more) and the odors-per-KC
    distances (bands in baseline-noise units)."""
    E106 = _e106()
    r = pd.read_csv(_note(RESULTS / E106.EXPERIMENT / 'realism.csv'))
    r = r[(r.family == 'family') & r.folder.notna()].copy()
    cols = ['rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv']
    mean = r[cols].mean().to_dict()
    r = pd.concat([r[['reading'] + cols], pd.DataFrame([dict(reading='mean of the eight', **mean)])], ignore_index=True)
    bands = dict(reading='band', rate_on=f'{AE.RATE_ON_BAND[0]}-{AE.RATE_ON_BAND[1]}',
                 ratio=f'{AE.RATE_RATIO_BAND[0]}-{AE.RATE_RATIO_BAND[1]}', rate_rho_on=f'>= {AE.RATE_RHO_MIN:.2f}',
                 opk_on_tv=f"<= {AE.opk_band('on', 'noise'):.3f}", opk_off_tv=f"<= {AE.opk_band('off', 'noise'):.3f}")
    ratio, _ = AE.data_rate_ratio('noise')
    flies = dict(reading='flies (baseline-noise units; rates the mean of the per-odor response fractions)',
                 rate_on=float(AE.frac_data('on', 'noise').mean()), rate_off=float(AE.frac_data('off', 'noise').mean()),
                 ratio=ratio)
    return pd.concat([r, pd.DataFrame([bands, flies])], ignore_index=True)


# ---------------------------------------------------------------------------
# the parts

def part_a() -> Dict[str, pd.DataFrame]:
    out = out_dir()
    point = part_a_table()
    boot = resample_part_a()
    a = with_intervals(point, boot)
    a.to_csv(out / 'a_two_sets.csv', index=False)
    boot.to_csv(out / 'a_resamples.csv.gz', index=False)
    hv = halves_table()
    hv.to_csv(out / 'a_halves.csv', index=False)
    cl, sh = classes_tables()
    cl.to_csv(out / 'a_classes.csv', index=False)
    sh.to_csv(out / 'a_class_shares.csv', index=False)
    v2 = validation2_table()
    v2.to_csv(out / 'a_validation2.csv', index=False)
    v2s = validation2_shares()
    v2s.to_csv(out / 'a_validation2_shares.csv', index=False)
    realism_table().to_csv(out / 'a_realism.csv', index=False)
    return dict(a=a, halves=hv, classes=cl, shares=sh, v2=v2, v2s=v2s)


def part_b() -> Dict[str, pd.DataFrame]:
    out = out_dir()
    crow = change_rows(CONTROLS)
    csum = change_summary(CONTROLS, crow)
    krow = change_rows(CANDIDATES)
    ksum = change_summary(CANDIDATES, krow)
    pd.concat([crow.assign(role='control'), krow.assign(role='candidate')]).to_csv(out / 'b_wirings.csv', index=False)
    pd.concat([csum.assign(role='control'), ksum.assign(role='candidate')]).to_csv(out / 'b_changes.csv', index=False)
    ys = yardstick_summary(crow, csum)
    ys.to_csv(out / 'b_yardstick.csv', index=False)
    pl = place(ksum, csum)
    pl.to_csv(out / 'b_placed.csv', index=False)
    dv = direction_values(krow, ksum)
    pv = price_verdict(dv, csum)
    (out / 'b_directions.json').write_text(json.dumps(dict(directions=dv, price=pv), indent=1, default=float) + '\n')
    return dict(rows=crow, controls=csum, cand_rows=krow, cands=ksum, yardstick=ys, placed=pl, directions=dv, price=pv)


def registry_rows(a: pd.DataFrame, ys: pd.DataFrame) -> pd.DataFrame:
    """The registry rows of both yardsticks (registry.csv)."""
    fam = a[(a.reader == 'family') & (a.remy == 'corr')].set_index('condition')
    rows = []
    fdesc = 'the best family (the breadth-2 APL-effect fit on its eight wirings, E72; raw counts)'
    for c in CONDITIONS:
        r = fam.loc[c]
        rows.append(dict(name=f'family_rho_{c}_two_sets_66', value=r['rho_two_sets'], lo=r['rho_two_sets_lo'],
                         hi=r['rho_two_sets_hi'], definition=f'ρ on odor-odor distances of {fdesc}, APL {c}, against the mean '
                         f"rank of two sets of flies (Yang's {DATA_LABEL[c]} and Remy's wild-type megamat KC somata, correlation "
                         'distance, 18 flies) over the 66 megamat pairs; interval the 5th to 95th percentiles over 1000 joint '
                         "resamples of Yang's recordings, Remy's flies and the eight wirings, seed 0 (E115)"))
        rows.append(dict(name=f'family_resid_corr_remy_{c}_66', value=r['resid_corr'], lo=r['resid_corr_lo'], hi=r['resid_corr_hi'],
                         definition=f"Spearman ρ over the 66 megamat pairs between the rank residuals of {fdesc}, APL {c}, "
                         f"against Yang's {DATA_LABEL[c]} and the rank residuals of Remy's wild-type somata (correlation "
                         'distance) against the same flies: how far the family\'s megamat errors are shared with a second set '
                         'of wild-type flies; interval as family_rho_<on|off>_two_sets_66 (E115)'))
    y = ys[ys.controls == 'all'].set_index('condition')
    for c in CONDITIONS:
        r = y.loc[c]
        what = (f'the {int(r["d_T_n"])} control input changes at fixed scales with both response counts re-held on the '
                "family's wirings (E93's fuller input, E102's steps at pn_inhadd 100 and 314, E97's permuted and inverted claw "
                "proportions, E108's permuted and scaled paa, E109's 30 permuted and inverted factors and its control searches) "
                f"of the change in response rate ρ with APL {c} (the mean over T's four wirings P3, P4, M2, M3, each against "
                "the family's own run on the wiring)")
        rows.append(dict(name=f'rate_rho_{c}_control_change_mean', value=r['d_T_mean'], lo=r['d_T_p5'], hi=r['d_T_p95'],
                         definition=f'mean over {what}; interval the 5th to 95th percentiles of the changes (E115)'))
        rows.append(dict(name=f'rate_rho_{c}_control_change_sd', value=r['d_T_sd'], lo=np.nan, hi=np.nan,
                         definition=f'SD (ddof 1) over {what} (E115)'))
        rows.append(dict(name=f'rate_rho_{c}_control_paired_sd', value=r['paired_sd_within'], lo=np.nan, hi=np.nan,
                         definition=f'paired SD over wirings of the same controls\' changes in response rate ρ with APL {c}: '
                         "the root mean square of each wiring's change about its control's mean over wirings (pooled over the "
                         'controls, ddof the number of controls) (E115)'))
    return pd.DataFrame(rows)


def registry_entries() -> List[Tuple[str, float, str, Optional[Sequence[float]]]]:
    """The two yardsticks for kcq/data_targets.py, read from results/E115_two_yardsticks/registry.csv: [(name, value,
    definition, 90 % interval or None)] (`family_rho_<on|off>_two_sets_66`, `family_resid_corr_remy_<on|off>_66`,
    `rate_rho_<on|off>_control_change_mean`, `_control_change_sd`, `rate_rho_<on|off>_control_paired_sd`)."""
    df = pd.read_csv(H.results_dir(EXPERIMENT) / 'registry.csv')
    out = []
    for _, r in df.iterrows():
        ci = None if pd.isna(r['lo']) else (float(r['lo']), float(r['hi']))
        out.append((r['name'], float(r['value']), r['definition'], ci))
    return out


# ---------------------------------------------------------------------------
# the comparison table

def _f(x, spec='.3f') -> str:
    try:
        return '–' if x is None or np.isnan(float(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def comparison_table(A: Mapping[str, pd.DataFrame], B: Mapping, V: Mapping) -> str:
    a = A['a']
    L = ['# E115: the megamat gap against two sets of flies and the yardstick of the price in response rate ρ', '',
         'No simulations. Raw counts (cosine on spike counts); Yang\'s flies the mean over recordings of the cosine distance; '
         "Remy's somata her wild-type megamat KC somata (correlation distance, 18 flies: the verdict's reading; cosine beside). "
         'Intervals: 5th to 95th percentiles over 1000 joint resamples (Yang\'s recordings, Remy\'s flies, the eight wirings).', '',
         '## Part A: the 66 megamat pairs', '',
         '| reader | Remy reading | APL | ρ vs Yang\'s flies | ρ vs Remy\'s flies | ρ vs the two sets\' mean rank | Remy vs Yang | '
         'residual correlation (ρ) | residual correlation (Pearson) |', '|---|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in a[a.group == 'reader'].iterrows():
        if r.reader.startswith('wiring') and r.remy == 'cos':
            continue
        iv = lambda k: f"{_f(r[k])} [{_f(r.get(k + '_lo'))}, {_f(r.get(k + '_hi'))}]" if (k + '_lo') in r else _f(r[k])
        L.append(f"| {r.reader} | {r.remy} | {r.condition} | {iv('rho_yang')} | {iv('rho_remy')} | {iv('rho_two_sets')} | "
                 f"{iv('remy_vs_yang')} | {iv('resid_corr')} | {_f(r.resid_corr_pearson)} |")
    L += ['', '### Beside: residual correlation with Yang\'s recordings split into halves (model vs half 1, Remy vs half 2; 200 splits)', '',
          '| reader | Remy reading | APL | mean | 5th | 95th |', '|---|---|---|---:|---:|---:|']
    for _, r in A['halves'].iterrows():
        L.append(f'| {r.reader} | {r.remy} | {r.condition} | {_f(r.resid_corr_halves)} | {_f(r.lo)} | {_f(r.hi)} |')
    L += ['', '### The pairs by class (family, Remy correlation reading)', '',
          '| APL | class | pairs | share of 1 − ρ |', '|---|---|---:|---:|']
    sh = A['shares']
    for _, r in sh[(sh.level == 'class') & (sh.remy == 'corr')].iterrows():
        L.append(f'| {r.condition} | {r.cls} | {r.n_pairs} | {_f(r.share_of_gap)} |')
    L += ['', '### Saved candidates (66 pairs, Remy correlation reading)', '',
          '| reader | wirings | APL | ρ vs Yang | ρ vs Remy | ρ vs two sets | residual correlation |', '|---|---|---|---:|---:|---:|---:|']
    for _, r in a[a.group.str.startswith('candidate') & (a.remy == 'corr')].iterrows():
        L.append(f'| {r.reader} | {r.wirings} | {r.condition} | {_f(r.rho_yang)} | {_f(r.rho_remy)} | {_f(r.rho_two_sets)} | '
                 f'{_f(r.resid_corr)} |')
    L += ['', '### The validation2 odors (the family on Remy\'s validation2 panel)', '',
          '| part | APL | Remy matrix | pairs | n | ρ vs Remy | ρ vs Yang | resid. corr. |', '|---|---|---|---|---:|---:|---:|---:|']
    for _, r in A['v2'].iterrows():
        L.append(f"| {r.part} | {r.condition} | {r.remy} | {r.pairs} | {_f(r.n_pairs, '.0f')} | {_f(r.rho_remy)} | "
                 f"{_f(r.get('rho_yang'))} | {_f(r.get('resid_corr'))} |")
    L += ['', '## Part B: the yardstick of response rate ρ', '',
          '| controls | APL | n | mean change (T) | SD (T) | 5th | 95th | pooled mean (T) | pooled SD (T) | paired SD |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in B['yardstick'].iterrows():
        L.append(f'| {r.controls} | {r.condition} | {r.d_T_n} | {_f(r.d_T_mean, "+.3f")} | {_f(r.d_T_sd)} | {_f(r.d_T_p5, "+.3f")} | '
                 f'{_f(r.d_T_p95, "+.3f")} | {_f(r.pooled_d_T_mean, "+.3f")} | {_f(r.pooled_d_T_sd)} | {_f(r.paired_sd_within)} |')
    L += ['', '### The candidates placed among the controls (mean over T\'s four wirings)', '',
          '| candidate | APL | change | size | percentile among controls | controls 5th | 95th |', '|---|---|---:|---:|---:|---:|---:|']
    for _, r in B['placed'][B['placed'].unit == 'T'].iterrows():
        L.append(f'| {r.source} {r.change} | {r.condition} | {_f(r.change_value, "+.3f")} | {_f(r["size"])} | {_f(r.percentile, ".0f")} | '
                 f'{_f(r.controls_p5, "+.3f")} | {_f(r.controls_p95, "+.3f")} |')
    L += ['', '### E109\'s directions', '', '| direction | size | fall, APL on | pooled fall, APL on | fall, APL off | pooled, APL off |',
          '|---|---:|---:|---:|---:|---:|']
    for d, v in B['directions'].items():
        L.append(f"| {d} | {_f(v['size'])} | {_f(v['d_rate_rho_on'], '+.4f')} | {_f(v['pooled_d_rate_rho_on'], '+.4f')} | "
                 f"{_f(v['d_rate_rho_off'], '+.4f')} | {_f(v['pooled_d_rate_rho_off'], '+.4f')} |")
    L += ['', '## Verdicts', '']
    for k in ('shared', 'megamat_level', 'price', 'yang_alone'):
        L.append(f"- {k}: {'holds' if V[k]['holds'] else 'does not hold'}")
    for d, v in V['price']['directions'].items():
        L.append(f"  - direction {d}: fall {_f(v['fall'], '+.4f')} vs controls' 5th–95th [{_f(v['p5'], '+.4f')}, {_f(v['p95'], '+.4f')}] "
                 f"(n {v['n_used']}{', every sized control' if v['fallback_all_sized'] else ', comparable size'}; percentile "
                 f"{_f(v['percentile'], '.0f')}); pooled fall {_f(v['pooled_fall'], '+.4f')}")
    L.append(f"- outcomes: {', '.join(V['outcomes'])}" + (f" ({V['outcome_note']})" if V['outcome_note'] else ''))
    return '\n'.join(L) + '\n'


def provenance() -> dict:
    return dict(experiment=EXPERIMENT, script='experiments/E115_two_yardsticks.py', host=socket.gethostname(),
                when=datetime.now().strftime('%Y-%m-%d %H:%M'), n_boot=N_BOOT, seed=SEED, n_halves=N_HALVES,
                rules=dict(shared_min=SHARED_MIN, level_tol=LEVEL_TOL, pooled_max_fall=POOLED_MAX_FALL, class_near=CLASS_NEAR,
                           class_far=CLASS_FAR, comparable=COMPARABLE, min_comparable=MIN_COMPARABLE, carry_share=CARRY_SHARE),
                approx_conc_note=getattr(H, 'APPROX_CONC_NOTE', ''), inputs=dict(_INPUTS), **H.code_commits())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('0', 'A', 'B', 'all'), default='all')
    args = ap.parse_args()
    out = out_dir()
    print('part 0:', json.dumps(part0(), default=float))
    if args.part == '0':
        return
    A = part_a() if args.part in ('A', 'all') else None
    B = part_b() if args.part in ('B', 'all') else None
    if A is not None and B is not None:
        V = verdicts(A['a'], A['classes'], B['price'])
        (out / 'verdict.json').write_text(json.dumps(V, indent=1, default=float) + '\n')
        registry_rows(A['a'], B['yardstick']).to_csv(out / 'registry.csv', index=False)
        (out / 'comparison_table.md').write_text(comparison_table(A, B, V))
        print(json.dumps(V, indent=1, default=float))
    (out / 'provenance.json').write_text(json.dumps(provenance(), indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
