#!/usr/bin/env python3
"""E106: the gap in ρ on odor-odor distances read by what could close it: the two conditions, the odors, the odor part
and the pair part, and the yardsticks for "not worse" (no simulations).

The entry is item E106 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08; moved into the list 2026-11-17 16:20). "The two ρ's" are ρ on
odor-odor distances with APL on against the control flies and with APL off against APL output blocked (TNT), each KC
matrix the mean over recordings of the cosine distance (`kcq.harness.yang_mean_rdm`), the model's the cosine distance on
spike counts (raw counts, the saved model_pair_distances_long.csv). "The family" is the best family, the breadth-2
APL-effect fit (C2, E67) on its eight wirings (the hemibrain wiring and seven re-drawn; E72's part B), read as the mean
matrix over the wirings per condition (E72's population); "the ceiling family" its four wirings with E69's firing
ceiling (E72's part G: E69's fit on the real wiring, part G's on P1, P2 and M1).

Steps (the entry's design):
1. every archived variant: the two ρ's over the 120 pairs and three subsets (the 66 megamat pairs, the 54 with paa,
   -aPine, euc or B-myr, the 91 without paa and B-myr), written into results/apl_effect/by_variant.csv by
   `kcq.apl_effect` (columns `spike_counts_cosine__dist_rho_<on|off>_<subset>`); this script reads them back
   (archive.csv, one row per distinct variant; archive_summary.csv);
2. the family and the ceiling family: the split between the conditions (each condition rebuilt from the other
   condition's distances and an APL effect, the model's or the flies'; conditions.csv), E61's rank residuals by pair and
   by odor (residuals_pairs.csv, residuals_odors.csv), the odor part and the pair part of each matrix with both hybrids
   (odor_pair_parts.csv), each odor left out and each odor's pairs (and the 29 pairs of paa and B-myr) put at the data's
   values (in residuals_odors.csv and odors.csv);
3. the model's ORN input read the same way, as Pearson and as cosine distance (input.csv, and the input's rows in the
   step-2 tables): its two ρ's by subset beside the family's, so that what the KC stage adds over its input is one column;
4. the yardsticks, as registry entries (`registry_entries`, read by kcq/data_targets.py) and in yardsticks.csv: the
   ceilings by subset of pairs, one recording against the mean of the others, Remy's wild-type KC somata against both
   matrices (the registry's), the SD over wirings, the resampling SD and the split-half difference of the family's ρ's,
   the paired SD of a near-null change (E93's fuller input, E102's weakest step) and the agreement of two populations of
   E85's model flies whose inputs come from disjoint ORN flies;
5. the trade-off table (tradeoff.csv): the two ρ's against the APL-effect loss, its parts and the bands' quantities over
   the archive; the inputs searched for the pattern part (E77's, read in E83's part A2; E83's arm I; E86's candidates)
   read for the two ρ's on their held-out wirings (searched_inputs.csv);
6. robustness (robustness.csv): against the distances without the trial noise (E81) and against the six typical
   recordings with APL output blocked (TNT) (E76: all eight but the two with the lowest mean distance).
probes.csv sets each of the planning probes beside the script's value.

The public interface other experiments import (E107, E108; stable signatures, pair Series as kcq.harness makes them, an
explicit list of pairs or odors): `select_pairs`, `rho_pairs` (ρ on odor-odor distances over chosen pairs), `two_rhos`
(the two ρ's at once), `rank_residuals` and `gap_shares_by_odor` (E61's shares of the gap by odor).

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no table of Yang's own Claude's connectome project, no
FlyWire / FAFB / BANC data and no raw-only recording is read. Archived variants whose wirings were built from such
tables (`kcq.apl_effect.HELD_VARIANTS`: E97's points of 2026-11-15, and E4's two FAFB wirings) get no new reading and
are left out of every table here.

"At the data's values" puts a pair's model distance at the flies' distance ('raw') or at the model's own distance of the
flies' percentile ('quantile': the flies' rank of the pair mapped onto the model's distribution); the hybrids of the odor
part and the pair part add the parts as they are ('raw') or with the donated part rescaled to the SD of the part it
replaces ('matched').

Outputs: <shared>/results/E106_rho_gap_parts/ (the CSVs above, comparison_table.md, provenance.json).
Usage: ./venv/bin/python -u experiments/E106_rho_gap_parts.py
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats as sstats

from kcq import harness as H
from kcq import apl_effect as AE

EXPERIMENT = 'E106_rho_gap_parts'
RESULTS = H.shared_dir() / 'results'
CONDITIONS = ('on', 'off')
GENOTYPE = {'on': 'TNT3c', 'off': 'TNT'}
DATA_LABEL = {'on': 'control flies', 'off': 'APL output blocked (TNT)'}
SUBSETS = tuple(AE.PAIR_SUBSETS)
SUBSET_N = {'all120': 120, 'megamat66': 66, 'v2odors54': 54, 'no_paa_bmyr91': 91}
SUBSET_SHORT = {'all120': '120', 'megamat66': '66', 'v2odors54': '54', 'no_paa_bmyr91': '91'}
PAA_BMYR = ('paa', 'B-myr')
ODORS = list(H.APL16_ODORS)
V = AE.V
N_BOOT, N_SPLITS = 1000, 500
# the near-null changes read on the same wirings (the entry's yardstick): (name, folder, variant tag by wiring, base)
# with base 'family' (the family's own fit on the wiring) or a tag pattern of the same folder
NEAR_NULL = (
    ('e93_fuller_input', "E93's fuller input F, held reading", 'E93_missing_glomeruli_input', 'B_F_{w}', 'family'),
    ('e102_weakest_step', "E102's weakest step (inhibition added at level 100, held reading) against its identity step",
     'E102_antennal_lobe_normalization', 'B_inhadd100_{w}_held', 'B_identity_{w}_held'),
)
# the inputs searched for the pattern part, read on held-out wirings (the four T wirings, P3, P4, M2, M3: E77's and E83's
# searches were on S, real P1 P2 M1; E86's candidates likewise): (label, folder, runs table, part, change)
SEARCHED = (
    ("E77's searched input, half A, size 0.6 (E83 part A2)", 'E83_cross_wiring_fit', 'runs.csv', 'A2', 'A2_halfA_s0.6'),
    ("E77's searched input, half A, size 1 (E83 part A2)", 'E83_cross_wiring_fit', 'runs.csv', 'A2', 'A2_halfA_s1'),
    ("E77's searched input, half B, size 0.6 (E83 part A2)", 'E83_cross_wiring_fit', 'runs.csv', 'A2', 'A2_halfB_s0.6'),
    ("E77's searched input, half B, size 1 (E83 part A2)", 'E83_cross_wiring_fit', 'runs.csv', 'A2', 'A2_halfB_s1'),
    ("E83 arm I, half A, size 0.6", 'E83_cross_wiring_fit', 'runs.csv', 'DI', 'DI_halfA_s0.6'),
    ("E83 arm I, half B, size 0.6", 'E83_cross_wiring_fit', 'runs.csv', 'DI', 'DI_halfB_s0.6'),
    ("E83 arm I, half A, size 1", 'E83_cross_wiring_fit', 'runs.csv', 'DI', 'DI_halfA_s1'),
    ("E83 arm I, half B, size 1", 'E83_cross_wiring_fit', 'runs.csv', 'DI', 'DI_halfB_s1'),
    ("E86 candidate A", 'E86_reduced_input_candidate', 'runs.csv', 'A2', 'cand_A'),
    ("E86 candidate B", 'E86_reduced_input_candidate', 'runs.csv', 'A2', 'cand_B'),
    ("E86 C1, half A", 'E86_reduced_input_candidate', 'runs.csv', 'Cread', 'C1_A'),
    ("E86 C1, half B", 'E86_reduced_input_candidate', 'runs.csv', 'Cread', 'C1_B'),
    ("E86 C2, half C", 'E86_reduced_input_candidate', 'runs.csv', 'Cread', 'C2_C'),
    ("E86 C2, half D", 'E86_reduced_input_candidate', 'runs.csv', 'Cread', 'C2_D'),
)
HELD_OUT = ('P3', 'P4', 'M2', 'M3')
E85_STAGES = ('B1S', 'B1T')
# the planning probes (the entry's "What the probes found"), for probes.csv: key -> (what, probe value)
PROBES = {
    'n_distinct': ('distinct archived variants', 2031),
    'max_on': ('largest ρ with APL on in the archive (none above)', 0.73),
    'max_off': ('largest ρ with APL off in the archive (none above)', 0.79),
    'n_both_070_077': ('variants with 0.70 with APL on and 0.77 with APL off together', 0),
    'max_on_outside': ('largest ρ with APL on outside the bands', 0.71),
    'max_off_outside': ('largest ρ with APL off outside the bands', 0.77),
    'n_inside': ('variants inside the bands', 1608),
    'corr_loss_on_inside': ('rank correlation of ρ with APL on with the APL-effect loss, inside the bands', -0.34),
    'corr_loss_off_inside': ('rank correlation of ρ with APL off with the APL-effect loss, inside the bands', -0.36),
    'family_on': ('the family, ρ with APL on', 0.69),
    'family_off': ('the family, ρ with APL off', 0.76),
    'cond_on_flies_off_family_effect': ("APL on rebuilt from the flies' APL-off distances and the family's APL effect", 0.85),
    'cond_on_family_off_flies_effect': ("APL on rebuilt from the family's APL-off distances and the flies' APL effect", 0.62),
    'cond_resid_on_off': ("ρ between the two conditions' rank residuals", 0.79),
    'family_vs_orn_on': ('the family against the ORN matrix, APL on', 0.93),
    'family_vs_orn_off': ('the family against the ORN matrix, APL off', 0.95),
    'flies_vs_orn_on': ('the control flies against the ORN matrix', 0.65),
    'flies_vs_orn_off': ('APL output blocked (TNT) against the ORN matrix', 0.73),
    'resid_family_orn_on': ("the family's rank residuals against the ORN input's, APL on", 0.89),
    'resid_family_orn_off': ("the family's rank residuals against the ORN input's, APL off", 0.88),
    'input_cosine_on': ("the cosine distance of the model's own input vectors against the control flies", 0.70),
    'input_cosine_off': ("the cosine distance of the model's own input vectors against APL output blocked (TNT)", 0.79),
    'family_91_on': ('the family on the 91 pairs without paa and B-myr, APL on', 0.79),
    'family_91_off': ('the family on the 91 pairs without paa and B-myr, APL off', 0.87),
    'orn_91_on': ('the ORN input on the 91 pairs, against the control flies', 0.76),
    'orn_91_off': ('the ORN input on the 91 pairs, against APL output blocked (TNT)', 0.86),
    'at_data_29_on': ("the family with the 29 pairs of paa and B-myr at the data's values, APL on", 0.83),
    'at_data_29_off': ("the family with the 29 pairs of paa and B-myr at the data's values, APL off", 0.89),
    'gain91_halves_min': ('gain without paa and B-myr against a half of the recordings, lower end (5th percentile over '
                          'halves, the smaller condition)', 0.08),
    'gain91_halves_max': ('gain without paa and B-myr against a half of the recordings, upper end (95th percentile, the '
                          'larger condition)', 0.12),
    'gain91_recordings_positive': ('recordings (of 19) against which the gain without paa and B-myr is positive', 19),
    'hyb_flies_odor_on': ("the flies' odor part with the family's pair part, APL on", 0.74),
    'hyb_flies_odor_off': ("the flies' odor part with the family's pair part, APL off", 0.79),
    'hyb_family_odor_on': ("the family's odor part with the flies' pair part, APL on", 0.89),
    'hyb_family_odor_off': ("the family's odor part with the flies' pair part, APL off", 0.96),
    'per_odor_worth_on': ("what per-odor handles are worth at most (the first hybrid minus the family), APL on", 0.05),
    'per_odor_worth_off': ("what per-odor handles are worth at most (the first hybrid minus the family), APL off", 0.03),
    'family_66_on': ('the family on the 66 megamat pairs, APL on', 0.77),
    'family_66_off': ('the family on the 66 megamat pairs, APL off', 0.82),
    'orn_66_on': ('the ORN input on the 66 megamat pairs, against the control flies', 0.78),
    'orn_66_off': ('the ORN input on the 66 megamat pairs, against APL output blocked (TNT)', 0.85),
    'remy_66_on': ("Remy's wild-type KC somata against the control flies (66 pairs)", 0.85),
    'remy_66_off': ("Remy's wild-type KC somata against APL output blocked (TNT) (66 pairs)", 0.87),
    'sd_wirings_on': ('SD over the eight wirings of ρ with APL on', 0.02),
    'sd_wirings_off': ('SD over the eight wirings of ρ with APL off', 0.02),
    'boot_sd_on': ("SD of the family's ρ with APL on over resamples of the recordings", 0.02),
    'boot_sd_off': ("SD of the family's ρ with APL off over resamples of the recordings", 0.02),
    'split_diff_min': ("difference of the family's ρ between two halves of the recordings, lower end (mean absolute "
                       "difference, the smaller condition)", 0.04),
    'split_diff_max': ("difference of the family's ρ between two halves of the recordings, upper end (mean absolute "
                       "difference, the larger condition)", 0.06),
    'paired_sd_min': ('paired SD of a near-null change over the wirings, smallest (of both changes and conditions)', 0.004),
    'paired_sd_max': ('paired SD of a near-null change over the wirings, largest', 0.009),
    'e85_agreement_66_on': ("two populations of E85's model flies from disjoint ORN flies, 66 megamat pairs, APL on", 0.96),
    'e85_agreement_66_off': ("two populations of E85's model flies from disjoint ORN flies, 66 megamat pairs, APL off", 0.96),
    'searched_on_min': ('change in ρ with APL on of the searched inputs on held-out wirings (population of the four), '
                        'smallest (size 0.6 and E86)', -0.08),
    'searched_on_max': ('change in ρ with APL on of the searched inputs on held-out wirings (population of the four), '
                        'largest (size 0.6 and E86)', 0.00),
    'noise_free_on': ('the family against the noise-free control flies', 0.67),
    'noise_free_off': ('the family against the noise-free APL output blocked (TNT)', 0.75),
}

# other readings of a probe, set beside it in probes.csv: key -> [(reading, value key in `probe_values`)]
PROBE_ALTS = {
    'max_off': [('inside the bands', 'max_off_inside')], 'max_on': [('inside the bands', 'max_on_inside')],
    'at_data_29_on': [('quantile', 'at_data_29_on_quantile')], 'at_data_29_off': [('quantile', 'at_data_29_off_quantile')],
    **{f'hyb_{w}_odor_{c}': [(r, f'hyb_{w}_odor_{c}_{r}') for r in ('matched', 'quantile', 'best')]
       for w in ('flies', 'family') for c in CONDITIONS},
    **{f'per_odor_worth_{c}': [(r, f'per_odor_worth_{c}_{r}') for r in ('matched', 'quantile', 'best', 'optimal')]
       for c in CONDITIONS},
    'gain91_halves_min': [('smallest of all halves', 'gain91_halves_min_abs')],
    'gain91_halves_max': [('largest of all halves', 'gain91_halves_max_abs')],
    'split_diff_min': [('root mean square', 'split_rms_min')], 'split_diff_max': [('root mean square', 'split_rms_max')],
    'searched_on_min': [('all searched inputs, size 1 included', 'searched_on_min_all')],
    'searched_on_max': [('all searched inputs, size 1 included', 'searched_on_max_all')],
    'orn_66_on': [("the 66-pair reference table's ORN column (registry rho_ctrl_orn_66)", 'orn_66_on_ref66')],
    'orn_66_off': [("the 66-pair reference table's ORN column (registry rho_tnt_orn_66)", 'orn_66_off_ref66')],
}

def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


_INPUTS: Dict[str, str] = {}


def _note(p: Path) -> Path:
    """Records an input file's sha256 for the provenance."""
    p = Path(p)
    if str(p) not in _INPUTS:
        _INPUTS[str(p)] = _sha(p)
    return p


# ---------------------------------------------------------------------------
# pairs, subsets and the splits' algebra (made-up matrices in the tests)

_ORDER = {o: i for i, o in enumerate(H.APL16_ODORS)}


def _canon(x: pd.Series) -> pd.Series:
    """A pair Series in canonical pair order (`kcq.harness.canonical_pairs`) when it holds every pair of its odors,
    else (a subset of pairs) with each key put in canonical orientation and the order kept."""
    try:
        return H.canonical_pairs(x)
    except KeyError:
        keys = [(a, b) if _ORDER.get(a, len(_ORDER)) <= _ORDER.get(b, len(_ORDER)) else (b, a) for a, b in x.index]
        return pd.Series(x.to_numpy(), index=pd.MultiIndex.from_tuples(keys, names=['odor_a', 'odor_b']), name=x.name)


def subset_pairs(index: pd.Index, subset: str) -> pd.Index:
    """The pairs of `index` in `subset` (`kcq.apl_effect.PAIR_SUBSETS`), in `index`'s order."""
    keep = AE.PAIR_SUBSETS[subset]
    return index[np.array([bool(keep(*p)) for p in index], dtype=bool)]


def _sub(s: pd.Series, subset: str) -> pd.Series:
    return s.loc[subset_pairs(s.index, subset)]


def rho(model: pd.Series, data: pd.Series, subset: str = 'all120') -> float:
    """ρ on odor-odor distances of `model` against `data` over the pairs of `subset` (pairs either lacks dropped)."""
    return H.spearman(_sub(_canon(model), subset), _canon(data))


# --- the public interface other experiments build on (E107, E108): pair Series as kcq.harness makes them (index
# (odor_a, odor_b), either orientation), and an explicit list of pairs or odors

def select_pairs(index: pd.Index, pairs: Optional[Sequence[Tuple[str, str]]] = None, *, odors: Optional[Sequence[str]] = None,
                 without: Optional[Sequence[str]] = None) -> pd.Index:
    """The pairs of a pair index to read: those in `pairs` (a list of (a, b), either orientation), and / or those with
    both odors in `odors`, and / or those with neither odor in `without` (all of `index` when none is given), in
    `index`'s order. The index is taken in canonical orientation (`kcq.harness.canonical_pairs`' keys)."""
    idx = _canon(pd.Series(np.zeros(len(index)), index=index)).index
    keep = np.ones(len(idx), dtype=bool)
    if pairs is not None:
        want = {tuple(sorted(p, key=lambda o: _ORDER.get(o, len(_ORDER)))) for p in pairs}
        keep &= np.array([p in want for p in idx], dtype=bool)
    if odors is not None:
        keep &= np.array([p[0] in odors and p[1] in odors for p in idx], dtype=bool)
    if without is not None:
        keep &= np.array([p[0] not in without and p[1] not in without for p in idx], dtype=bool)
    return idx[keep]


def rho_pairs(model: pd.Series, data: pd.Series, pairs: Optional[Sequence[Tuple[str, str]]] = None, *,
              odors: Optional[Sequence[str]] = None, without: Optional[Sequence[str]] = None) -> float:
    """ρ on odor-odor distances (Spearman over pairs) of the pair Series `model` against the pair Series `data`, over the
    pairs `select_pairs` keeps (`pairs`, `odors`, `without`; all pairs both have when none is given); pairs either side
    lacks or has NaN for are dropped. E.g. the 66 megamat pairs: `odors=kcq.harness.MEGAMAT12`; the 91 without paa and
    B-myr: `without=('paa', 'B-myr')`."""
    m = _canon(model)
    return H.spearman(m.loc[select_pairs(m.index, pairs, odors=odors, without=without)], _canon(data))


def two_rhos(on: pd.Series, off: pd.Series, data_on: Optional[pd.Series] = None, data_off: Optional[pd.Series] = None,
             pairs: Optional[Sequence[Tuple[str, str]]] = None, *, odors: Optional[Sequence[str]] = None,
             without: Optional[Sequence[str]] = None) -> Dict[str, float]:
    """The two ρ's on odor-odor distances, {'on': APL on against `data_on`, 'off': APL off against `data_off`}, over the
    pairs `select_pairs` keeps (`rho_pairs`). `data_on` / `data_off` default to the control flies' and APL output blocked
    (TNT)'s mean matrices (`data_matrices`: cosine, mean over recordings)."""
    if data_on is None or data_off is None:
        D = data_matrices()
        data_on = D['on'] if data_on is None else data_on
        data_off = D['off'] if data_off is None else data_off
    kw = dict(pairs=pairs, odors=odors, without=without)
    return {'on': rho_pairs(on, data_on, **kw), 'off': rho_pairs(off, data_off, **kw)}


def gap_shares_by_odor(model: pd.Series, data: pd.Series, pairs: Optional[Sequence[Tuple[str, str]]] = None, *,
                       odors: Optional[Sequence[str]] = None, without: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """E61's reading of the gap by odor, over the pairs `select_pairs` keeps: ρ = 1 - 6 sum(d^2) / (n (n^2 - 1)) with d a
    pair's rank residual (`rank_residuals`: its rank in `model` minus its rank in `data`, ranks over those pairs), so each
    pair's contribution to 1 - ρ is 6 d^2 / (n (n^2 - 1)) (exact without ties). Per odor (index; the odors of those pairs):
    `n_pairs`, `mean_abs_resid`, `mean_signed_resid` (over the odor's pairs), `share_of_gap` (the odor's pairs'
    contributions over 1 - ρ; each pair counts for both its odors, so the shares add to 2), `rho` (over all the pairs
    kept), `loo_rho` (ρ over them without the odor's pairs) and `loo_delta` (loo_rho - rho)."""
    m = _canon(model)
    m = m.loc[select_pairs(m.index, pairs, odors=odors, without=without)]
    r = rank_residuals(m, data)
    full = H.spearman(r['model_distance'], r['data_distance'])
    ods = list(dict.fromkeys([o for p in r.index for o in p]))
    ods = sorted(ods, key=lambda o: _ORDER.get(o, len(_ORDER)))
    rows = {}
    for o in ods:
        idx = [p for p in r.index if o in p]
        rest = [p for p in r.index if o not in p]
        loo = H.spearman(r.loc[rest, 'model_distance'], r.loc[rest, 'data_distance']) if len(rest) > 2 else np.nan
        rows[o] = dict(n_pairs=len(idx), mean_abs_resid=float(r.loc[idx, 'rank_resid'].abs().mean()),
                       mean_signed_resid=float(r.loc[idx, 'rank_resid'].mean()),
                       share_of_gap=float(r.loc[idx, 'contribution'].sum() / (1 - full)) if full < 1 else np.nan,
                       rho=full, loo_rho=loo, loo_delta=loo - full)
    out = pd.DataFrame.from_dict(rows, orient='index')
    out.index.name = 'odor'
    return out


def rank_residuals(model: pd.Series, data: pd.Series) -> pd.DataFrame:
    """Per pair (E61): the ranks over the pairs both have (ties at the average rank), the rank residual d (model minus
    data) and the contribution 6 d^2 / (n (n^2 - 1)) to 1 - ρ (exact without ties)."""
    m, d = H.align_pairs(_canon(model), _canon(data))
    rm, rd = sstats.rankdata(m.to_numpy(dtype=float)), sstats.rankdata(d.to_numpy(dtype=float))
    n = len(m)
    dd = rm - rd
    return pd.DataFrame(dict(model_distance=m.to_numpy(dtype=float), data_distance=d.to_numpy(dtype=float), rank_model=rm,
        rank_data=rd, rank_resid=dd, contribution=6 * dd ** 2 / (n * (n ** 2 - 1))), index=m.index)


def _aligned(*xs: pd.Series) -> List[pd.Series]:
    xs = [_canon(x) for x in xs]
    idx = xs[0].index
    for x in xs[1:]:
        idx = idx.intersection(x.index, sort=False)
    ok = np.ones(len(idx), dtype=bool)
    for x in xs:
        ok &= x.reindex(idx).notna().to_numpy()
    idx = idx[ok]
    return [x.reindex(idx) for x in xs]


def condition_split(m_on: pd.Series, m_off: pd.Series, d_on: pd.Series, d_off: pd.Series,
                    subset: str = 'all120') -> Dict[str, float]:
    """The split between the conditions, over the pairs of `subset`: each condition's ρ; APL on rebuilt from the flies'
    APL-off distances plus the model's per-pair APL effect (m_on - m_off) and from the model's APL-off distances plus the
    flies' effect (d_on - d_off), each against the control flies, and APL off rebuilt likewise (d_on - model effect,
    m_on - flies' effect) against APL output blocked (TNT); each model condition against the other condition's flies;
    the flies' two matrices against each other; ρ between the two conditions' rank residuals."""
    m_on, m_off, d_on, d_off = _aligned(m_on, m_off, d_on, d_off)
    eff_m, eff_d = m_on - m_off, d_on - d_off
    r = lambda a, b: rho(a, b, subset)
    res_on = rank_residuals(_sub(m_on, subset), _sub(d_on, subset))['rank_resid']
    res_off = rank_residuals(_sub(m_off, subset), _sub(d_off, subset))['rank_resid']
    return dict(rho_on=r(m_on, d_on), rho_off=r(m_off, d_off),
        on_from_flies_off_model_effect=r(d_off + eff_m, d_on), on_from_model_off_flies_effect=r(m_off + eff_d, d_on),
        off_from_flies_on_model_effect=r(d_on - eff_m, d_off), off_from_model_on_flies_effect=r(m_on - eff_d, d_off),
        model_off_vs_flies_on=r(m_off, d_on), model_on_vs_flies_off=r(m_on, d_off), flies_on_vs_flies_off=r(d_on, d_off),
        model_on_vs_model_off=r(m_on, m_off), resid_on_vs_off=H.spearman(res_on, res_off.reindex(res_on.index)))


def odor_pair_parts(pairs: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """A pair Series split into its odor part (the least-squares fit of value(a, b) = s_a + s_b over the pairs,
    `kcq.harness.odor_shares`; the design spans the constant, so the odor part carries the mean) and its pair part (the
    rest, orthogonal to every odor's column): (odor part per pair, pair part per pair, s per odor)."""
    s, _ = H.odor_shares(pairs)
    X, odors = H.pair_design(pairs.index, list(s.index))
    fit = pd.Series(X @ s.to_numpy(dtype=float), index=pairs.index, name='odor_part')
    return fit, (pairs - fit).rename('pair_part'), s


def _hybrid(odor_part: pd.Series, pair_part: pd.Series, own_part: Optional[pd.Series] = None, which: str = 'pair') -> pd.Series:
    """odor part + pair part; with `own_part` the donated one (`which`) rescaled to `own_part`'s SD over the pairs."""
    if own_part is None:
        return odor_part + pair_part
    if which == 'pair':
        sd = float(pair_part.std(ddof=0))
        return odor_part + (pair_part * float(own_part.std(ddof=0)) / sd if sd > 0 else pair_part)
    c = odor_part.mean()
    sd = float(odor_part.std(ddof=0))
    return c + (odor_part - c) * (float(own_part.std(ddof=0)) / sd if sd > 0 else 1.0) + pair_part


def quantile_map(x: pd.Series, onto: pd.Series) -> pd.Series:
    """`x` mapped onto the distribution of `onto` (same pairs): each pair's value the quantile of `onto` at the pair's
    percentile in `x` (average rank over n, `np.quantile`'s linear rule); keeps `x`'s ranks."""
    q = (sstats.rankdata(x.to_numpy(dtype=float)) - 1) / max(len(x) - 1, 1)
    return pd.Series(np.quantile(onto.to_numpy(dtype=float), q), index=x.index)


LAMBDAS = np.round(np.linspace(0.0, 4.0, 81), 3)


def _best_scale(base: pd.Series, part: pd.Series, data: pd.Series) -> Tuple[float, float]:
    """(the largest ρ against `data` of base + λ part over λ in `LAMBDAS`, that λ)."""
    vals = [H.spearman(base + lam * part, data) for lam in LAMBDAS]
    i = int(np.nanargmax(vals))
    return float(vals[i]), float(LAMBDAS[i])


def optimal_odor_part(pair_part: pd.Series, data: pd.Series) -> pd.Series:
    """The odor part (a value s_a + s_b per pair) that, added to `pair_part` (orthogonal to the odor design, as
    `odor_pair_parts` returns it), maximizes the Pearson correlation of the sum with the ranks of `data`: the projection of
    the centred ranks onto the odor design, scaled by |p|^2 / (p . y) (p the pair part, y the centred ranks; the maximum of
    (t |y_X| + p . y) / sqrt(t^2 + |p|^2) over the projection's length t). The Pearson-optimal per-odor handle, a stand-in
    for the Spearman-optimal one."""
    y = pd.Series(sstats.rankdata(data.reindex(pair_part.index).to_numpy(dtype=float)), index=pair_part.index)
    y = y - y.mean()
    yx, _, _ = odor_pair_parts(y)
    p = pair_part - pair_part.mean()
    py = float(p @ y)
    scale = float(p @ p) / py if py > 0 else 1.0
    return yx * scale


def odor_pair_hybrids(model: pd.Series, data: pd.Series) -> Dict[str, float]:
    """Each matrix split into its odor part and pair part (`odor_pair_parts`, over the 120 pairs): ρ against `data` of
    the model, of the flies' odor part plus the model's pair part (`flies_odor_model_pair`), of the model's odor part plus
    the flies' pair part (`model_odor_flies_pair`), of the model's odor part alone and of its pair part alone; the same
    hybrids with the donated part rescaled to the SD of the part it replaces (`_matched`), with the flies' matrix mapped
    onto the model's distribution first (`_quantile`, `quantile_map`), with the donated part at the scale that gives the
    largest ρ (`_best`, its scale `_best_scale`), and the model's pair part with the Pearson-optimal odor part
    (`model_pair_optimal_odor`, `optimal_odor_part`: what a per-odor handle could reach); ρ over the 16 odors between the
    two sides' odor parts (s) and over the pairs between their pair parts; the share of each side's variance over pairs
    that its odor part explains, and the parts' SDs."""
    m, d = _aligned(model, data)
    om, pm, sm = odor_pair_parts(m)
    od, pd_, sd = odor_pair_parts(d)
    odq, pdq, _ = odor_pair_parts(quantile_map(d, m))
    r = lambda x: H.spearman(x, d)
    best_fo, lam_fo = _best_scale(om.mean() + pm, od - od.mean(), d)
    best_mo, lam_mo = _best_scale(om, pd_, d)
    out = dict(model=r(m), flies_odor_model_pair=r(od + pm), model_odor_flies_pair=r(om + pd_),
        flies_odor_model_pair_matched=r(_hybrid(od, pm, own_part=om, which='odor')),
        model_odor_flies_pair_matched=r(_hybrid(om, pd_, own_part=pm, which='pair')),
        flies_odor_model_pair_quantile=r(odq + pm), model_odor_flies_pair_quantile=r(om + pdq),
        flies_odor_model_pair_best=best_fo, flies_odor_model_pair_best_scale=lam_fo,
        model_odor_flies_pair_best=best_mo, model_odor_flies_pair_best_scale=lam_mo,
        model_pair_optimal_odor=r(optimal_odor_part(pm, d) + pm),
        model_odor_part_only=r(om), model_pair_part_only=r(pm), odor_shares_rho=H.spearman(sm, sd.reindex(sm.index)),
        pair_parts_rho=H.spearman(pm, pd_), model_odor_r2=float(1 - pm.var(ddof=0) / m.var(ddof=0)),
        data_odor_r2=float(1 - pd_.var(ddof=0) / d.var(ddof=0)), model_pair_sd=float(pm.std(ddof=0)),
        data_pair_sd=float(pd_.std(ddof=0)), model_odor_sd=float((om - om.mean()).std(ddof=0)),
        data_odor_sd=float((od - od.mean()).std(ddof=0)))
    return out


def leave_one_odor_out(model: pd.Series, data: pd.Series) -> pd.Series:
    """Per odor: ρ over the pairs without the odor."""
    m, d = _aligned(model, data)
    odors = list(dict.fromkeys([o for p in m.index for o in p]))
    return pd.Series({o: H.spearman(m.loc[[p for p in m.index if o not in p]], d) for o in odors}, name='loo_rho')


def at_data_values(model: pd.Series, data: pd.Series, odors: Sequence[str], how: str = 'raw') -> pd.Series:
    """The model's distances with every pair involving one of `odors` put at the data's value: 'raw' the flies' distance
    itself, 'quantile' the model's own distance at the flies' percentile of the pair (its average rank among the flies'
    pairs, over n, mapped through the model's empirical quantiles, `np.quantile`'s linear rule)."""
    m, d = _aligned(model, data)
    hit = np.array([p[0] in odors or p[1] in odors for p in m.index], dtype=bool)
    out = m.copy()
    if how == 'raw':
        out[hit] = d[hit].to_numpy(dtype=float)
    elif how == 'quantile':
        out[hit] = quantile_map(d, m)[hit].to_numpy(dtype=float)
    else:
        raise ValueError(f'{how=} not raw or quantile')
    return out


# ---------------------------------------------------------------------------
# reading the archive

@functools.lru_cache(maxsize=None)
def _folder_pairs(folder: str) -> Dict[str, Tuple[pd.Series, pd.Series]]:
    p = RESULTS / folder / 'model_pair_distances_long.csv'
    if p.exists():
        _note(p)
    return AE.model_pairs(RESULTS / folder)


def load_pairs(folder: str, tag: str) -> Dict[str, pd.Series]:
    """{'on', 'off'}: a saved variant's pair distances (raw counts, cosine; results/<folder>/model_pair_distances_long.csv,
    as kcq.apl_effect reads them). Refuses the held variants."""
    if AE.held(folder, tag):
        raise ValueError(f'{folder}:{tag} is held (kcq.apl_effect.HELD_VARIANTS): no new reading')
    pr = _folder_pairs(folder).get(tag)
    if pr is None:
        raise KeyError(f'no saved distances for {folder}:{tag}')
    return {'on': pr[0], 'off': pr[1]}


@functools.lru_cache(maxsize=None)
def _members(name: str) -> Tuple[Tuple[str, str, str], ...]:
    from experiments import E72_wiring_population as E72
    _note(RESULTS / E72.EXPERIMENT / 'fits.csv')
    if name == 'family':
        return tuple(E72.new_fit_members())
    if name == 'ceiling':
        return tuple(E72.new_fit_members(prefix_tag=E72.G_SUFFIX, real_folder=E72.G_REAL[0], real_tag=E72.G_REAL[1],
                                         draws=E72.G_DRAWS))
    raise ValueError(f'{name=} not family or ceiling')


def family_members(name: str = 'family') -> List[Tuple[str, str, str]]:
    """[(wiring label, results folder, tag)] of the family ('family', eight wirings: `E72_wiring_population.
    new_fit_members`) or the ceiling family ('ceiling', four: E72's part G)."""
    return list(_members(name))


def mean_pairs(ds: Sequence[Mapping[str, pd.Series]]) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the mean over `ds` of their pair distances per condition (E72's population; pairs all have)."""
    out = {}
    for c in CONDITIONS:
        df = pd.concat([H.canonical_pairs(d[c]) for d in ds], axis=1)
        out[c] = df.dropna().mean(axis=1)
    return out


def population(members: Sequence[Tuple[str, str, str]]) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the mean over the members of their pair distances per condition (E72's population)."""
    return mean_pairs([load_pairs(f, t) for _, f, t in members])


@functools.lru_cache(maxsize=None)
def _recs() -> Dict[str, pd.DataFrame]:
    _note(H.yang_data_dir() / 'computed' / 'rdm_per_recording_long.csv')
    return {g: H.yang_recording_rdms(g, 'cosine') for g in ('TNT3c', 'TNT')}


def data_matrices() -> Dict[str, pd.Series]:
    """{'on': the control flies' mean matrix, 'off': APL output blocked (TNT)'s} (`kcq.harness.yang_mean_rdm`: the mean
    over recordings of the cosine distance)."""
    r = _recs()
    return {'on': r['TNT3c'].mean(), 'off': r['TNT'].mean()}


@functools.lru_cache(maxsize=None)
def _inputs() -> Tuple[Tuple[str, pd.Series], ...]:
    from kcq import data_targets as DT
    deltas = H.apl16_orn_deltas()
    out = (('orn_pearson', H.orn16_consensus_pairs()), ('orn_cosine', H.canonical_pairs(DT.orn_cosine_pairs())),
           ('input_pearson', H.pair_distances(deltas, 'pearson')), ('input_cosine', H.pair_distances(deltas, 'cosine')))
    _note(DT.INPUTS['orn']); _note(DT.INPUTS['orn_vectors'])
    return out


INPUT_LABEL = {'orn_pearson': "the ORN reference (1 - Pearson r of the olfdb fly-mean vectors, 36 glomeruli; the registry's)",
               'orn_cosine': '1 - cosine similarity of the same fly-mean vectors',
               'input_pearson': "1 - Pearson r of the model's own input (its ORN spike deltas, 36 glomeruli)",
               'input_cosine': "1 - cosine similarity of the model's own input (its ORN spike deltas)"}


def input_matrices() -> Dict[str, pd.Series]:
    """The model's ORN input read as distances: {name: pair Series} for the registry's ORN reference (`orn_pearson`, 1 -
    Pearson of the olfdb fly-mean vectors), the cosine of those vectors (`orn_cosine`), and 1 - Pearson and 1 - cosine of
    the model's own input (`input_pearson`, `input_cosine`: `kcq.harness.apl16_orn_deltas`, the spike deltas)."""
    return dict(_inputs())


STAT_COLS = [f'{V}__{k}' for k in ('apl_effect_rmse', 'apl_effect_rmse_centered', 'apl_effect_mean', 'apl_effect_spearman',
                                   'apl_effect_pearson', 'apl_effect_sd_model')]


def rho_col(cond: str, subset: str = 'all120') -> str:
    return f'{V}__dist_rho_{cond}_{subset}'


def archive(by_variant: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """by_variant.csv's rows with the two ρ's of step 1, the held variants out, one row per distinct variant: rows with
    the same APL-effect statistics and the same ρ's (to 12 significant digits) are the same saved distances in several
    folders or under several tags and are counted once (the first kept; `n_copies`, `copies`); rows with neither ρ (no
    120-pair reading) are left out. `inside` marks the variants inside the realism bands (no `realism_flag`)."""
    if by_variant is None:
        p = _note(AE.OUT / 'by_variant.csv')
        by_variant = pd.read_csv(p, low_memory=False)
    df = by_variant.copy()
    df = df[[not AE.held(f, t) for f, t in zip(df.folder, df.tag)]]
    cols = [rho_col(c, s) for c in CONDITIONS for s in SUBSETS]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f'by_variant.csv lacks {missing}: run `./venv/bin/python -m kcq.apl_effect` first')
    df = df[df[rho_col('on')].notna() | df[rho_col('off')].notna()]
    key = df[STAT_COLS + cols].astype(float).apply(lambda col: col.map(lambda x: f'{x:.12g}'))
    df = df.assign(_key=key.agg('|'.join, axis=1))
    df['realism_flag'] = df['realism_flag'].fillna('').astype(str)
    g = df.groupby('_key', sort=False)
    first = g.head(1).copy()
    first['n_copies'] = g['tag'].transform('size').loc[first.index].astype(int)
    first['copies'] = g.apply(lambda x: ' ; '.join(f'{a}:{b}' for a, b in zip(x.folder, x.tag)), include_groups=False
                              ).reindex(first['_key']).to_numpy()
    first['inside'] = first['realism_flag'] == ''
    return first.drop(columns='_key').reset_index(drop=True)


# ---------------------------------------------------------------------------
# the tables

def _readings(name: str) -> Tuple[Dict[str, Dict[str, pd.Series]], Dict[str, pd.Series]]:
    """({wiring label: {'on', 'off'}}, the population) of a family."""
    mem = family_members(name)
    per = {lab: load_pairs(f, t) for lab, f, t in mem}
    return per, mean_pairs(list(per.values()))


def family_tables(name: str = 'family') -> Dict[str, pd.DataFrame]:
    """Step 2 for one family: 'conditions' (`condition_split` per subset, the population and each wiring), 'pairs' (E61's
    rank residuals per pair of the population and of the ORN input), 'odors' (per odor: mean absolute rank residual, its
    share of 1 - ρ, ρ without the odor, ρ with its pairs at the data's values), 'parts' (`odor_pair_hybrids`, the
    population and each wiring), 'sets' (the 29 pairs of paa and B-myr together: left out, at the data's values)."""
    D = data_matrices()
    orn = input_matrices()['orn_pearson']
    per, pop = _readings(name)
    readings = {'population': pop, **per}
    cond_rows, part_rows = [], []
    for lab, m in readings.items():
        for s in SUBSETS:
            cond_rows.append(dict(family=name, reading=lab, subset=s, n_pairs=len(subset_pairs(H.canonical_pairs(m['on']).index, s)),
                                  **condition_split(m['on'], m['off'], D['on'], D['off'], s)))
        for c in CONDITIONS:
            part_rows.append(dict(family=name, reading=lab, condition=c, **odor_pair_hybrids(m[c], D[c])))
    pair_rows, odor_rows, set_rows = [], [], []
    for c in CONDITIONS:
        r = rank_residuals(pop[c], D[c])
        ro = rank_residuals(orn, D[c]).reindex(r.index)
        rho_full = H.spearman(pop[c], D[c])
        for (a, b), x in r.iterrows():
            pair_rows.append(dict(family=name, condition=c, odor_a=a, odor_b=b, pair=f'{a} / {b}',
                megamat=(a in H.MEGAMAT12 and b in H.MEGAMAT12), with_paa_bmyr=(a in PAA_BMYR or b in PAA_BMYR), **x.to_dict(),
                orn_rank_resid=float(ro.loc[(a, b), 'rank_resid']), orn_distance=float(ro.loc[(a, b), 'model_distance'])))
        gs = gap_shares_by_odor(pop[c], D[c])
        go = gap_shares_by_odor(orn, D[c])
        for o in ODORS:
            odor_rows.append(dict(family=name, condition=c, odor=o, **gs.loc[o].to_dict(),
                mean_abs_resid_orn=float(go.loc[o, 'mean_abs_resid']), share_of_gap_orn=float(go.loc[o, 'share_of_gap']),
                loo_rho_orn=float(go.loc[o, 'loo_rho']),
                at_data_raw=H.spearman(at_data_values(pop[c], D[c], [o], 'raw'), D[c]),
                at_data_quantile=H.spearman(at_data_values(pop[c], D[c], [o], 'quantile'), D[c])))
        for lab, odors in (('paa and B-myr', PAA_BMYR), ('the four validation2 odors', tuple(H.VALIDATION2_4))):
            keep = [p for p in r.index if p[0] not in odors and p[1] not in odors]
            set_rows.append(dict(family=name, condition=c, odors=lab, n_pairs=len(r) - len(keep), rho=rho_full,
                rho_without=H.spearman(pop[c].reindex(keep), D[c]), rho_without_orn=H.spearman(orn.reindex(keep), D[c]),
                at_data_raw=H.spearman(at_data_values(pop[c], D[c], odors, 'raw'), D[c]),
                at_data_quantile=H.spearman(at_data_values(pop[c], D[c], odors, 'quantile'), D[c])))
    return dict(conditions=pd.DataFrame(cond_rows), parts=pd.DataFrame(part_rows), pairs=pd.DataFrame(pair_rows),
                odors=pd.DataFrame(odor_rows), sets=pd.DataFrame(set_rows))


def input_table() -> Dict[str, pd.DataFrame]:
    """Step 3: 'rhos' (each input reading's two ρ's by subset beside the family's and the ceiling family's, and the KC
    stage's gain over the input, family minus input), 'follow' (the family's matrices against each input matrix, and its
    rank residuals against the input's), 'parts' (`odor_pair_hybrids` and `condition_split` of each input reading, the
    input standing for both conditions)."""
    D = data_matrices()
    fam = {n: _readings(n)[1] for n in ('family', 'ceiling')}
    rows, follow, parts = [], [], []
    for name, x in input_matrices().items():
        for s in SUBSETS:
            row = dict(input=name, what=INPUT_LABEL[name], subset=s)
            for c in CONDITIONS:
                row[f'rho_{c}'] = rho(x, D[c], s)
                for fn, f in fam.items():
                    row[f'{fn}_rho_{c}'] = rho(f[c], D[c], s)
                    row[f'{fn}_gain_{c}'] = row[f'{fn}_rho_{c}'] - row[f'rho_{c}']
            rows.append(row)
        for c in CONDITIONS:
            for fn, f in fam.items():
                rf = rank_residuals(f[c], D[c])['rank_resid']
                ri = rank_residuals(x, D[c])['rank_resid'].reindex(rf.index)
                follow.append(dict(input=name, family=fn, condition=c, family_vs_input=H.spearman(f[c], x),
                                   resid_family_vs_input=H.spearman(rf, ri)))
            parts.append(dict(input=name, condition=c, **odor_pair_hybrids(x, D[c])))
    flies = [dict(input=name, flies_on_vs_input=H.spearman(D['on'], x), flies_off_vs_input=H.spearman(D['off'], x))
             for name, x in input_matrices().items()]
    return dict(rhos=pd.DataFrame(rows), follow=pd.DataFrame(follow).merge(pd.DataFrame(flies), on='input'),
                parts=pd.DataFrame(parts))


TRADE_COLS = {
    'apl_effect_loss': f'{V}__apl_effect_rmse', 'pattern_part': f'{V}__apl_effect_rmse_centered',
    'mean_apl_effect': f'{V}__apl_effect_mean', 'abs_mean_error': None, 'amplitude_part': f'{V}__apl_effect_amplitude',
    'shape_part': f'{V}__apl_effect_shape', 'spread_ratio': f'{V}__apl_effect_sd_ratio',
    'rho_per_pair_apl_effect': f'{V}__apl_effect_spearman', 'rate_on': 'rate_on', 'rate_ratio': None,
    'rate_rho_on': 'rate_rho_on', 'opk_on_tv': 'opk_on_tv', 'opk_off_tv': 'opk_off_tv'}


def _trade_values(a: pd.DataFrame, k: str) -> pd.Series:
    if k == 'abs_mean_error':
        return a[f'{V}__apl_effect_mean_error'].abs()
    if k == 'rate_ratio':
        return a['rate_off'] / a['rate_on']
    return a[TRADE_COLS[k]]


def tradeoff_table(arch: pd.DataFrame) -> pd.DataFrame:
    """Step 5: rank correlations over the archive's distinct variants (inside the bands and all) between each of the two
    ρ's (120 pairs and the subsets) and the APL-effect loss, its parts and the bands' quantities."""
    rows = []
    for scope, a in (('inside the bands', arch[arch.inside]), ('all', arch), ('outside the bands', arch[~arch.inside])):
        for c in CONDITIONS:
            for s in SUBSETS:
                x = a[rho_col(c, s)].astype(float)
                for k in TRADE_COLS:
                    y = _trade_values(a, k).astype(float)
                    ok = x.notna() & y.notna()
                    rows.append(dict(scope=scope, condition=c, subset=s, against=k, n=int(ok.sum()),
                        rank_corr=float(sstats.spearmanr(x[ok], y[ok])[0]) if ok.sum() > 2 else np.nan))
    return pd.DataFrame(rows)


def archive_summary(arch: pd.DataFrame) -> pd.DataFrame:
    """Step 1's summary: per scope (all, inside, outside the bands) and subset, the number of distinct variants, the
    largest ρ with APL on and with APL off and which variants hold them, how many reach given levels together, and
    outcome (c)'s count: variants whose two ρ's on the 66 megamat pairs are both 0.03 or more above the ORN input's."""
    D = data_matrices()
    orn = input_matrices()['orn_pearson']
    rows = []
    for scope, a in (('all', arch), ('inside the bands', arch[arch.inside]), ('outside the bands', arch[~arch.inside])):
        for s in SUBSETS:
            on, off = a[rho_col('on', s)].astype(float), a[rho_col('off', s)].astype(float)
            i_on, i_off = on.idxmax(), off.idxmax()
            row = dict(scope=scope, subset=s, n=len(a), max_on=float(on.max()), max_on_variant=f"{a.folder[i_on]}:{a.tag[i_on]}",
                       max_off=float(off.max()), max_off_variant=f"{a.folder[i_off]}:{a.tag[i_off]}",
                       median_on=float(on.median()), median_off=float(off.median()),
                       orn_on=rho(orn, D['on'], s), orn_off=rho(orn, D['off'], s))
            row['n_above_orn_both'] = int(((on > row['orn_on']) & (off > row['orn_off'])).sum())
            row['n_above_orn_both_003'] = int(((on >= row['orn_on'] + 0.03) & (off >= row['orn_off'] + 0.03)).sum())
            if s == 'all120':
                row['n_both_070_077'] = int(((on >= 0.70) & (off >= 0.77)).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def outcome_c(arch: pd.DataFrame) -> pd.DataFrame:
    """Outcome (c)'s candidates: variants whose ρ's on the 66 megamat pairs both exceed their input's by 0.03 or more,
    against the ORN reference and against the cosine of the model's own input (the larger of each input's ρ's)."""
    D = data_matrices()
    ins = input_matrices()
    base = {c: max(rho(ins[n], D[c], 'megamat66') for n in ('orn_pearson', 'input_pearson', 'orn_cosine', 'input_cosine'))
            for c in CONDITIONS}
    on, off = arch[rho_col('on', 'megamat66')], arch[rho_col('off', 'megamat66')]
    hit = (on >= base['on'] + 0.03) & (off >= base['off'] + 0.03)
    out = arch.loc[hit, ['exp', 'folder', 'tag', rho_col('on', 'megamat66'), rho_col('off', 'megamat66'), rho_col('on'),
                         rho_col('off'), f'{V}__apl_effect_rmse', 'realism_flag']].copy()
    out['input_on_66_max'], out['input_off_66_max'] = base['on'], base['off']
    return out


def _runs_tags(folder: str, table: str, part: str, change: str) -> Dict[str, str]:
    p = _note(RESULTS / folder / table)
    r = pd.read_csv(p, low_memory=False)
    r = r[(r['part'].astype(str) == part) & (r['change'].astype(str) == change)]
    return dict(zip(r['wiring'].astype(str), r['tag'].astype(str)))


def searched_inputs_table() -> pd.DataFrame:
    """Step 5: the two ρ's of the inputs searched for the pattern part (E77's, E83's arm I, E86's candidates) on the
    held-out wirings (P3, P4, M2, M3), per wiring and as the population over them, each paired with the family's own fit
    on the same wirings (the lab's input; the held readings reproduce those fits exactly, E83's identity)."""
    D = data_matrices()
    fam = {lab: (f, t) for lab, f, t in family_members('family')}
    rows = []
    for label, folder, table, part, change in SEARCHED:
        tags = _runs_tags(folder, table, part, change)
        ws = [w for w in HELD_OUT if w in tags and tags[w] in _folder_pairs(folder)]
        if not ws:
            continue
        var = {w: load_pairs(folder, tags[w]) for w in ws}
        base = {w: load_pairs(*fam[w]) for w in ws}
        for w, reading in list((w, w) for w in ws) + [('population', 'population')]:
            v = mean_pairs([var[x] for x in ws]) if w == 'population' else var[w]
            b = mean_pairs([base[x] for x in ws]) if w == 'population' else base[w]
            row = dict(input=label, folder=folder, change=change, wiring=reading, n_wirings=len(ws) if w == 'population' else 1)
            for c in CONDITIONS:
                for s in ('all120', 'megamat66', 'no_paa_bmyr91'):
                    row[f'rho_{c}_{SUBSET_SHORT[s]}'] = rho(v[c], D[c], s)
                    row[f'base_rho_{c}_{SUBSET_SHORT[s]}'] = rho(b[c], D[c], s)
                    row[f'change_{c}_{SUBSET_SHORT[s]}'] = row[f'rho_{c}_{SUBSET_SHORT[s]}'] - row[f'base_rho_{c}_{SUBSET_SHORT[s]}']
            st = H.apl_effect_stats(v['on'], v['off'])
            sb = H.apl_effect_stats(b['on'], b['off'])
            row.update(pattern=st['apl_effect_rmse_centered'], base_pattern=sb['apl_effect_rmse_centered'],
                       loss=st['apl_effect_rmse'], base_loss=sb['apl_effect_rmse'])
            rows.append(row)
    return pd.DataFrame(rows)


def noise_free_data() -> Dict[str, pd.Series]:
    """The flies' noise-free mean matrices (E81's noise_free_pairs.csv, Yang's z units: distances from products of
    different repeats), {'on', 'off'}."""
    p = _note(RESULTS / 'E81_repeat_structure' / 'noise_free_pairs.csv')
    t = pd.read_csv(p)
    t = t[t['units'] == 'z'].set_index(['odor_a', 'odor_b'])
    return {'on': H.canonical_pairs(t['ctrl_noise_free']), 'off': H.canonical_pairs(t['tnt_noise_free'])}


def six_typical_tnt() -> Tuple[pd.Series, List[tuple]]:
    """The mean matrix of the six typical recordings with APL output blocked (TNT): all eight but the two with the lowest
    mean distance (E76's `strongest_recordings`, in Yang's z), and those two."""
    tnt = _recs()['TNT']
    two = list(tnt.mean(axis=1).sort_values().index[:2])
    return tnt.drop(index=two).mean(), two


def robustness_table() -> pd.DataFrame:
    """Step 6: the family's and the ceiling family's ρ's (and the ORN input's) against the noise-free distances (E81) and,
    with APL off, against the six typical recordings with APL output blocked (TNT) (E76), by subset."""
    D = data_matrices()
    nf = noise_free_data()
    six, two = six_typical_tnt()
    targets = {'as measured': D, 'noise-free (E81)': nf, 'six typical TNT recordings (E76)': {'on': D['on'], 'off': six}}
    readings = {n: _readings(n)[1] for n in ('family', 'ceiling')}
    orn = input_matrices()['orn_pearson']
    readings['ORN input'] = {'on': orn, 'off': orn}
    rows = []
    for tname, T in targets.items():
        for rname, m in readings.items():
            for s in SUBSETS:
                rows.append(dict(target=tname, reading=rname, subset=s, rho_on=rho(m['on'], T['on'], s),
                                 rho_off=rho(m['off'], T['off'], s), excluded_tnt=' + '.join('_'.join(map(str, r)) for r in two)
                                 if tname.startswith('six') else ''))
    return pd.DataFrame(rows)


def realism_table(arch_full: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """The realism of the two families (for the write-up's Realism section): each member wiring's row of by_variant.csv
    (response rates with APL on and off, their ratio, the response rate ρ, the odors-per-KC distances, the flag; the KC
    data's side in baseline-noise units, as the bands read it) and the mean over the members; the mean pair distances of
    the population (raw counts) beside the KC data's; the KC data's values and the bands."""
    bv = pd.read_csv(AE.OUT / 'by_variant.csv', low_memory=False) if arch_full is None else arch_full
    bv = bv.set_index(['folder', 'tag'])
    cols = ['rate_on', 'rate_off', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv']
    rows = []
    for name in ('family', 'ceiling'):
        mem = family_members(name)
        sub = []
        for lab, f, t in mem:
            r = bv.loc[(f, t)]
            d = dict(family=name, reading=lab, folder=f, tag=t, **{c: float(r[c]) for c in cols},
                     realism_flag=r['realism_flag'] if isinstance(r['realism_flag'], str) else '')
            d['ratio'] = d['rate_off'] / d['rate_on']
            sub.append(d)
        rows += sub
        m = {c: float(np.mean([x[c] for x in sub])) for c in cols + ['ratio']}
        pop = _readings(name)[1]
        rows.append(dict(family=name, reading='mean over the wirings', folder='', tag='', **m,
                         realism_flag='; '.join(sorted({x['realism_flag'] for x in sub if x['realism_flag']})),
                         mean_dist_on=float(pop['on'].mean()), mean_dist_off=float(pop['off'].mean())))
    D = data_matrices()
    rows.append(dict(family='KC data', reading='control flies / APL output blocked (TNT), baseline-noise units',
                     rate_on=H.data_target('noise_ctrl_frac'), rate_off=H.data_target('noise_tnt_frac'),
                     ratio=H.data_target('noise_rate_ratio_fraction'), opk_on_tv=AE.opk_band('on', 'noise'),
                     opk_off_tv=AE.opk_band('off', 'noise'), mean_dist_on=float(D['on'].mean()), mean_dist_off=float(D['off'].mean()),
                     realism_flag=f'bands: APL-on rate {AE.RATE_ON_BAND[0]}-{AE.RATE_ON_BAND[1]}, ratio {AE.RATE_RATIO_BAND[0]}-'
                                  f'{AE.RATE_RATIO_BAND[1]}, response rate ρ with APL on at least {AE.RATE_RHO_MIN:.2f}, '
                                  'odors per KC at most the opk_*_tv given'))
    return pd.DataFrame(rows)


def _split_means(rows: pd.DataFrame, rng: np.random.Generator) -> Tuple[pd.Series, pd.Series]:
    idx = rng.permutation(len(rows))
    return rows.iloc[idx[:len(rows) // 2]].mean(), rows.iloc[idx[len(rows) // 2:]].mean()


def gain_without_paa_bmyr(pop: Mapping[str, pd.Series], seed: int = 0) -> pd.DataFrame:
    """The family's gain on the 91 pairs without paa and B-myr over the 120 (ρ_91 - ρ_120) against each half of
    `N_SPLITS` random splits of each genotype's recordings and against each single recording."""
    rng = np.random.default_rng(seed)
    rows = []
    for c in CONDITIONS:
        R = _recs()[GENOTYPE[c]]
        for k in range(N_SPLITS):
            for h, half in zip('AB', _split_means(R, rng)):
                rows.append(dict(condition=c, kind='half', index=k, half=h, rho_120=rho(pop[c], half),
                                 rho_91=rho(pop[c], half, 'no_paa_bmyr91')))
        for rec, row in R.iterrows():
            rows.append(dict(condition=c, kind='recording', index='_'.join(map(str, rec)), half='',
                             rho_120=rho(pop[c], row), rho_91=rho(pop[c], row, 'no_paa_bmyr91')))
    out = pd.DataFrame(rows)
    out['gain'] = out['rho_91'] - out['rho_120']
    return out


def _e85_populations() -> Tuple[List[str], Dict[str, Dict[str, Dict[str, pd.Series]]]]:
    """E85's 'repro' model flies held on the eight wirings (stages B1S, B1T): (flies, {fly: {wiring: {'on', 'off'}}})."""
    folder = 'E85_orn_fly_population'
    runs = pd.read_csv(_note(RESULTS / folder / 'runs.csv'), low_memory=False)
    runs = runs[runs['stage'].isin(E85_STAGES) & (runs['variant'] == 'repro')]
    flies = sorted(runs['fly'].astype(str).unique(), key=lambda f: int(f[1:]))
    by = {f: {} for f in flies}
    for _, r in runs.iterrows():
        by[str(r['fly'])][str(r['wiring'])] = load_pairs(folder, str(r['tag']))
    return flies, by


def e85_agreement() -> pd.DataFrame:
    """Two populations of E85's model flies whose inputs come from disjoint ORN flies: every split of the nine ORN flies
    into four and five, each population the mean over its flies on all eight wirings; ρ between the two populations'
    mean matrices per condition and subset (the input differs between flies only on the 12 megamat odors)."""
    flies, by = _e85_populations()
    rows = []
    for a in itertools.combinations(flies, 4):
        b = [f for f in flies if f not in a]
        pa = mean_pairs([d for f in a for d in by[f].values()])
        pb = mean_pairs([d for f in b for d in by[f].values()])
        for c in CONDITIONS:
            for s in SUBSETS:
                rows.append(dict(flies_a=' '.join(a), flies_b=' '.join(b), condition=c, subset=s,
                                 rho=H.spearman(_sub(pa[c], s), pb[c])))
    return pd.DataFrame(rows)


def near_null_changes() -> pd.DataFrame:
    """The paired change of the two ρ's for each of `NEAR_NULL` on each wiring it was read on (the variant minus its base
    on the same wiring), per subset."""
    D = data_matrices()
    fam = {lab: (f, t) for lab, f, t in family_members('family')}
    rows = []
    for key, what, folder, pat, base in NEAR_NULL:
        have = _folder_pairs(folder)
        for w in fam:
            tag = pat.format(w=w)
            if tag not in have:
                continue
            v = load_pairs(folder, tag)
            b = load_pairs(*fam[w]) if base == 'family' else load_pairs(folder, base.format(w=w))
            for s in SUBSETS:
                rows.append(dict(change=key, what=what, wiring=w, subset=s, **{f'change_{c}': rho(v[c], D[c], s) - rho(b[c], D[c], s)
                                                                                for c in CONDITIONS}))
    return pd.DataFrame(rows)


def family_resampling(pop: Mapping[str, pd.Series], seed: int = 0) -> Dict[str, Dict[str, np.ndarray]]:
    """The family's ρ's against `N_BOOT` resamples of the recordings (with replacement within each genotype,
    `kcq.harness.resample_recording_means`) and their absolute differences between the halves of `N_SPLITS` random
    splits: {'boot' | 'split_diff': {'<cond>_<subset>': values}}."""
    rng = np.random.default_rng(seed)
    recs = _recs()
    boot = {f'{c}_{s}': [] for c in CONDITIONS for s in SUBSETS}
    for _ in range(N_BOOT):
        m = H.resample_recording_means(recs, rng)
        for c in CONDITIONS:
            for s in SUBSETS:
                boot[f'{c}_{s}'].append(rho(pop[c], m[GENOTYPE[c]], s))
    diff = {f'{c}_{s}': [] for c in CONDITIONS for s in SUBSETS}
    for _ in range(N_SPLITS):
        for c in CONDITIONS:
            a, b = _split_means(recs[GENOTYPE[c]], rng)
            for s in SUBSETS:
                diff[f'{c}_{s}'].append(abs(rho(pop[c], a, s) - rho(pop[c], b, s)))
    return dict(boot={k: np.asarray(v) for k, v in boot.items()}, split_diff={k: np.asarray(v) for k, v in diff.items()})


@functools.lru_cache(maxsize=None)
def _family_rs() -> Dict[str, Dict[str, np.ndarray]]:
    return family_resampling(_readings('family')[1])


def recording_loo(rows: pd.DataFrame, subset: str) -> np.ndarray:
    """Per recording: ρ of its pair distances with the mean of the other recordings of its genotype, over `subset`."""
    cols = subset_pairs(rows.columns, subset)
    R = rows.loc[:, cols]
    return np.array([H.spearman(R.iloc[i], R.drop(index=R.index[i]).mean()) for i in range(len(R))])


@functools.lru_cache(maxsize=None)
def _registry_values() -> Tuple[Tuple[str, float, str, Optional[Tuple[float, float]]], ...]:
    out = []
    recs = _recs()
    for g, name in (('TNT3c', 'ctrl'), ('TNT', 'tnt')):
        lab = DATA_LABEL['on' if g == 'TNT3c' else 'off']
        for s in ('all120', 'megamat66', 'v2odors54', 'no_paa_bmyr91'):
            n = SUBSET_SHORT[s]
            R = recs[g].loc[:, subset_pairs(recs[g].columns, s)]
            sh = H.split_half_reliability(R, n_splits=200, seed=0)
            if s in ('v2odors54', 'no_paa_bmyr91'):
                out.append((f'{name}_ceiling_{n}', sh['ceiling'], f'measurement ceiling on odor-odor distances of the {lab} '
                            f'matrix on the {SUBSET_N[s]} pairs {_subset_words(s)} (sqrt of the Spearman-Brown-corrected '
                            'split-half ρ, 200 splits, seed 0; as ctrl_ceiling)', None))
                out.append((f'{name}_split_half_{n}', sh['split_half'], f'ρ between the mean matrices of random halves of the '
                            f'{lab} recordings, {SUBSET_N[s]} pairs {_subset_words(s)} (200 splits)', (sh['split_half_lo'],
                            sh['split_half_hi'])))
            if s != 'megamat66':
                loo = recording_loo(recs[g], s)
                out.append((f'{name}_rec_rho_loo_{n}', float(np.mean(loo)), f"mean over the {lab} recordings of one "
                            f"recording's ρ on odor-odor distances with the mean of the other recordings of its genotype "
                            f"({SUBSET_N[s]} pairs {_subset_words(s)}; interval the 5th to 95th percentiles over "
                            'recordings; ctrl_rec_rho_loo / tnt_rec_rho_loo are the 66 megamat pairs)',
                            tuple(np.percentile(loo, [5, 95]))))
    per, pop = _readings('family')
    D = data_matrices()
    fam = 'the best family (the breadth-2 APL-effect fit on its eight wirings, E72; raw counts)'
    for c in CONDITIONS:
        for s in SUBSETS:
            v = np.array([rho(m[c], D[c], s) for m in per.values()])
            out.append((f'family_rho_{c}_sd_wirings_{SUBSET_SHORT[s]}', float(np.std(v, ddof=1)), f'SD (ddof 1) over the eight '
                        f"wirings of {fam} of each wiring's ρ on odor-odor distances, APL {c} against the {DATA_LABEL[c]}, "
                        f'{SUBSET_N[s]} pairs {_subset_words(s)} (E106)', None))
    rs = _family_rs()
    for c in CONDITIONS:
        for s in ('all120', 'megamat66', 'no_paa_bmyr91'):
            k, n = f'{c}_{s}', SUBSET_SHORT[s]
            b = rs['boot'][k]
            out.append((f'family_rho_{c}_boot_sd_{n}', float(np.std(b, ddof=1)), f'SD over {N_BOOT} resamples of the recordings '
                        f'(with replacement within each genotype, seed 0) of the ρ on odor-odor distances of the population '
                        f'of {fam}, APL {c} against the {DATA_LABEL[c]}, {SUBSET_N[s]} pairs {_subset_words(s)}; interval '
                        'the 5th to 95th percentiles of the ρ (E106)', tuple(np.percentile(b, [5, 95]))))
            d = rs['split_diff'][k]
            out.append((f'family_rho_{c}_split_diff_{n}', float(np.mean(d)), f'mean over {N_SPLITS} random splits of the '
                        f'{DATA_LABEL[c]} recordings into halves of the absolute difference between the ρ on odor-odor '
                        f'distances of the population of {fam} against one half and against the other, APL {c}, '
                        f'{SUBSET_N[s]} pairs {_subset_words(s)}; interval the 5th to 95th percentiles (E106)',
                        tuple(np.percentile(d, [5, 95]))))
    nn = near_null_changes()
    for key, what, *_ in NEAR_NULL:
        g = nn[(nn.change == key) & (nn.subset == 'all120')]
        for c in CONDITIONS:
            v = g[f'change_{c}'].to_numpy(dtype=float)
            out.append((f'paired_rho_{c}_sd_{key}', float(np.std(v, ddof=1)), f'paired SD of a near-null change: SD (ddof 1) '
                        f'over the {len(v)} wirings it was read on of the change in ρ on odor-odor distances (120 pairs, APL '
                        f'{c} against the {DATA_LABEL[c]}) of {what}, each wiring against its base on the same wiring (E106)',
                        None))
            out.append((f'paired_rho_{c}_mean_{key}', float(np.mean(v)), f'mean over the same wirings of the same change (E106)',
                        None))
    ag = e85_agreement()
    for c in CONDITIONS:
        for s in ('megamat66', 'all120'):
            v = ag[(ag.condition == c) & (ag.subset == s)]['rho'].to_numpy(dtype=float)
            out.append((f'e85_population_agreement_{c}_{SUBSET_SHORT[s]}', float(np.mean(v)), "ρ on odor-odor distances between "
                        "two populations of E85's model flies (the 'repro' flies, held reading, on the eight wirings) whose "
                        'inputs come from disjoint ORN flies (four and five of the nine megamat ORN flies), APL '
                        f'{c}, {SUBSET_N[s]} pairs {_subset_words(s)}; mean over the {len(v)} splits, interval the 5th to 95th '
                        'percentiles (E106)', tuple(np.percentile(v, [5, 95]))))
    return tuple(out)


def _subset_words(s: str) -> str:
    return {'all120': 'of the panel', 'megamat66': 'of the 12 megamat odors',
            'v2odors54': 'with paa, -aPine, euc or B-myr', 'no_paa_bmyr91': 'without paa and B-myr'}[s]


# the registry's input keys (kcq/data_targets.py INPUTS) each entry reads, by name prefix
REGISTRY_INPUTS = {'ctrl_': ['kc_rdms'], 'tnt_': ['kc_rdms'],
                   'family_': ['kc_rdms', 'e106_e67_pairs', 'e106_e72_pairs', 'e106_e72_fits'],
                   'paired_': ['kc_rdms', 'e106_e67_pairs', 'e106_e72_pairs', 'e106_e72_fits', 'e106_e93_pairs', 'e106_e102_pairs'],
                   'e85_': ['e106_e85_pairs', 'e106_e85_runs']}


def registry_inputs(name: str) -> List[str]:
    """The registry's input keys an entry of `registry_entries` reads."""
    for k, v in REGISTRY_INPUTS.items():
        if name.startswith(k):
            return v
    raise KeyError(name)


def registry_entries() -> List[Tuple[str, float, str, Optional[Sequence[float]]]]:
    """Step 4's yardsticks for kcq/data_targets.py: [(name, value, definition, 90 % interval or None)]: the ceilings and
    split halves on the 54 and 91 pairs (`<ctrl|tnt>_ceiling_<54|91>`; the 120 and the 66 are ctrl_ceiling and
    ctrl_ceiling_66), one recording against the mean of the others (`<ctrl|tnt>_rec_rho_loo_<120|54|91>`; the 66 are
    ctrl_rec_rho_loo), the family's SD over wirings (`family_rho_<on|off>_sd_wirings_<n>`), its resampling SD and
    split-half difference (`family_rho_<on|off>_boot_sd_<n>`, `_split_diff_<n>`), the paired SD and mean of the near-null
    changes (`paired_rho_<on|off>_<sd|mean>_<change>`) and the agreement of E85's populations
    (`e85_population_agreement_<on|off>_<66|120>`). Remy's somata against both matrices are the registry's
    rho_ctrl_remy_soma and rho_tnt_remy_soma."""
    return [tuple(x) for x in _registry_values()]


def yardsticks_table() -> pd.DataFrame:
    """yardsticks.csv: every registry entry of `registry_entries` and the registry's existing entries it sits beside."""
    rows = [dict(name=n, value=v, lo=np.nan if ci is None else ci[0], hi=np.nan if ci is None else ci[1], source='E106 (new entry)',
                 definition=d) for n, v, d, ci in registry_entries()]
    for n in ('ctrl_ceiling', 'tnt_ceiling', 'ctrl_ceiling_66', 'tnt_ceiling_66', 'ctrl_split_half_rho', 'tnt_split_half_rho',
              'ctrl_split_half_66', 'tnt_split_half_66', 'ctrl_rec_rho_loo', 'tnt_rec_rho_loo', 'rho_ctrl_remy_soma',
              'rho_tnt_remy_soma', 'rho_ctrl_orn', 'rho_tnt_orn', 'rho_ctrl_orn_66', 'rho_tnt_orn_66', 'rho_ctrl_tnt',
              'rho_ctrl_tnt_66'):
        try:
            v = H.data_target(n)
        except KeyError:
            continue
        try:
            lo, hi = H.data_target_ci(n)
        except (KeyError, TypeError, ValueError):
            lo, hi = np.nan, np.nan
        rows.append(dict(name=n, value=v, lo=lo, hi=hi, source='registry (existing)', definition=''))
    return pd.DataFrame(rows)


def _f(x, spec='.2f') -> str:
    return '–' if x is None or (isinstance(x, float) and np.isnan(x)) else format(x, spec)


def probe_values(T: Mapping[str, pd.DataFrame], arch: pd.DataFrame) -> Dict[str, float]:
    """The script's value of each of `PROBES` (and of `PROBE_ALTS`' other readings) from the tables."""
    summ = T['archive_summary'].set_index(['scope', 'subset'])
    s_all, s_out = summ.loc[('all', 'all120')], summ.loc[('outside the bands', 'all120')]
    s_in = summ.loc[('inside the bands', 'all120')]
    tr = T['tradeoff']
    tr_in = tr[(tr.scope == 'inside the bands') & (tr.subset == 'all120') & (tr.against == 'apl_effect_loss')].set_index('condition')
    cond = T['conditions']
    cp = cond[(cond.family == 'family') & (cond.reading == 'population') & (cond.subset == 'all120')].iloc[0]
    fol = T['input_follow']
    fo = fol[(fol.input == 'orn_pearson') & (fol.family == 'family')].set_index('condition')
    rh = T['input_rhos'].set_index(['input', 'subset'])
    sets = T['sets']
    st = sets[(sets.family == 'family') & (sets.odors == 'paa and B-myr')].set_index('condition')
    gain = T['gain91']
    gh = gain[gain.kind == 'half'].groupby('condition')['gain']
    parts = T['parts']
    pp = parts[(parts.family == 'family') & (parts.reading == 'population')].set_index('condition')
    ys = T['yardsticks'].set_index('name')['value']
    nn = T['near_null']
    nn120 = nn[nn.subset == 'all120']
    paired = [float(np.std(g[f'change_{c}'], ddof=1)) for _, g in nn120.groupby('change') for c in CONDITIONS]
    srch = T['searched_inputs']
    sp = srch[srch.wiring == 'population']
    sp06 = sp[~sp.change.astype(str).str.endswith('_s1')]
    rob = T['robustness']
    rnf = rob[(rob.target == 'noise-free (E81)') & (rob.reading == 'family') & (rob.subset == 'all120')].iloc[0]
    rs = _family_rs()
    rms = {c: float(np.sqrt(np.mean(rs['split_diff'][f'{c}_all120'] ** 2))) for c in CONDITIONS}
    mad = {c: float(np.mean(rs['split_diff'][f'{c}_all120'])) for c in CONDITIONS}
    v = dict(
        n_distinct=len(arch), max_on=s_all.max_on, max_off=s_all.max_off, n_both_070_077=s_all.n_both_070_077,
        max_on_inside=s_in.max_on, max_off_inside=s_in.max_off,
        max_on_outside=s_out.max_on, max_off_outside=s_out.max_off, n_inside=int(arch.inside.sum()),
        corr_loss_on_inside=tr_in.loc['on', 'rank_corr'], corr_loss_off_inside=tr_in.loc['off', 'rank_corr'],
        family_on=cp.rho_on, family_off=cp.rho_off, cond_on_flies_off_family_effect=cp.on_from_flies_off_model_effect,
        cond_on_family_off_flies_effect=cp.on_from_model_off_flies_effect, cond_resid_on_off=cp.resid_on_vs_off,
        family_vs_orn_on=fo.loc['on', 'family_vs_input'], family_vs_orn_off=fo.loc['off', 'family_vs_input'],
        flies_vs_orn_on=fo.loc['on', 'flies_on_vs_input'], flies_vs_orn_off=fo.loc['on', 'flies_off_vs_input'],
        resid_family_orn_on=fo.loc['on', 'resid_family_vs_input'], resid_family_orn_off=fo.loc['off', 'resid_family_vs_input'],
        input_cosine_on=rh.loc[('input_cosine', 'all120'), 'rho_on'], input_cosine_off=rh.loc[('input_cosine', 'all120'), 'rho_off'],
        family_91_on=rh.loc[('orn_pearson', 'no_paa_bmyr91'), 'family_rho_on'],
        family_91_off=rh.loc[('orn_pearson', 'no_paa_bmyr91'), 'family_rho_off'],
        orn_91_on=rh.loc[('orn_pearson', 'no_paa_bmyr91'), 'rho_on'], orn_91_off=rh.loc[('orn_pearson', 'no_paa_bmyr91'), 'rho_off'],
        at_data_29_on=st.loc['on', 'at_data_raw'], at_data_29_off=st.loc['off', 'at_data_raw'],
        at_data_29_on_quantile=st.loc['on', 'at_data_quantile'], at_data_29_off_quantile=st.loc['off', 'at_data_quantile'],
        gain91_halves_min=float(gh.quantile(0.05).min()), gain91_halves_max=float(gh.quantile(0.95).max()),
        gain91_halves_min_abs=float(gh.min().min()), gain91_halves_max_abs=float(gh.max().max()),
        gain91_recordings_positive=int((gain[gain.kind == 'recording']['gain'] > 0).sum()),
        family_66_on=rh.loc[('orn_pearson', 'megamat66'), 'family_rho_on'], family_66_off=rh.loc[('orn_pearson', 'megamat66'), 'family_rho_off'],
        orn_66_on=rh.loc[('orn_pearson', 'megamat66'), 'rho_on'], orn_66_off=rh.loc[('orn_pearson', 'megamat66'), 'rho_off'],
        orn_66_on_ref66=H.data_target('rho_ctrl_orn_66'), orn_66_off_ref66=H.data_target('rho_tnt_orn_66'),
        remy_66_on=H.data_target('rho_ctrl_remy_soma'), remy_66_off=H.data_target('rho_tnt_remy_soma'),
        sd_wirings_on=ys['family_rho_on_sd_wirings_120'], sd_wirings_off=ys['family_rho_off_sd_wirings_120'],
        boot_sd_on=ys['family_rho_on_boot_sd_120'], boot_sd_off=ys['family_rho_off_boot_sd_120'],
        split_diff_min=min(mad.values()), split_diff_max=max(mad.values()),
        split_rms_min=min(rms.values()), split_rms_max=max(rms.values()),
        paired_sd_min=min(paired), paired_sd_max=max(paired),
        e85_agreement_66_on=ys['e85_population_agreement_on_66'], e85_agreement_66_off=ys['e85_population_agreement_off_66'],
        searched_on_min=float(sp06['change_on_120'].min()), searched_on_max=float(sp06['change_on_120'].max()),
        searched_on_min_all=float(sp['change_on_120'].min()), searched_on_max_all=float(sp['change_on_120'].max()),
        noise_free_on=rnf.rho_on, noise_free_off=rnf.rho_off)
    for c in CONDITIONS:
        v[f'hyb_flies_odor_{c}'] = pp.loc[c, 'flies_odor_model_pair']
        v[f'hyb_family_odor_{c}'] = pp.loc[c, 'model_odor_flies_pair']
        for r in ('matched', 'quantile', 'best'):
            v[f'hyb_flies_odor_{c}_{r}'] = pp.loc[c, f'flies_odor_model_pair_{r}']
            v[f'hyb_family_odor_{c}_{r}'] = pp.loc[c, f'model_odor_flies_pair_{r}']
        v[f'per_odor_worth_{c}'] = pp.loc[c, 'flies_odor_model_pair'] - pp.loc[c, 'model']
        for r in ('matched', 'quantile', 'best'):
            v[f'per_odor_worth_{c}_{r}'] = pp.loc[c, f'flies_odor_model_pair_{r}'] - pp.loc[c, 'model']
        v[f'per_odor_worth_{c}_optimal'] = pp.loc[c, 'model_pair_optimal_odor'] - pp.loc[c, 'model']
    return {k: float(x) for k, x in v.items()}


COUNT_PROBES = ('n_distinct', 'n_inside', 'n_both_070_077', 'gain91_recordings_positive')


def probes_table(T: Mapping[str, pd.DataFrame], arch: pd.DataFrame) -> pd.DataFrame:
    """probes.csv: each planning probe beside the script's value (reading 'main reading'), and the other
    readings of `PROBE_ALTS` beside it; `agree` when they differ by at most 0.01 (counts: exactly; the archive's maxima:
    the script's is at most the probe's bound)."""
    got = probe_values(T, arch)
    rows = []

    def ok(k, s, p):
        if k in COUNT_PROBES:
            return int(round(s)) == int(p)
        if k in ('max_on', 'max_off'):
            return s <= p + 0.005
        return abs(s - p) <= 0.01 + 1e-9

    for k, (what, p) in PROBES.items():
        for reading, key in [('main reading', k)] + PROBE_ALTS.get(k, []):
            s = got[key]
            rows.append(dict(probe=k, what=what, reading=reading, probe_value=p, script_value=s, difference=s - p,
                             agree=bool(ok(k, s, p))))
    return pd.DataFrame(rows)


def comparison_table(T: Mapping[str, pd.DataFrame], probes: pd.DataFrame, out: Path) -> None:
    """comparison_table.md: the family beside the ORN input and the flies' ceilings, and the probes."""
    D = data_matrices()
    rh = T['input_rhos'].set_index(['input', 'subset'])
    L = ['# E106: the gap in ρ on odor-odor distances read by what could close it', '',
         'ρ on odor-odor distances, APL on (against the control flies) / APL off (against APL output blocked (TNT)); raw '
         'counts; the family is the mean matrix over its eight wirings.', '',
         '| statistic | ceiling (KC data) | ORN input | family | ceiling family | family minus input |',
         '|---|---:|---:|---:|---:|---:|']
    for s in SUBSETS:
        n = SUBSET_SHORT[s]
        cl = {c: H.data_target({'all120': f'{g}_ceiling', 'megamat66': f'{g}_ceiling_66'}.get(s, f'{g}_ceiling_{n}'))
              if s in ('all120', 'megamat66') else T['yardsticks'].set_index('name')['value'].get(f'{g}_ceiling_{n}', np.nan)
              for c, g in (('on', 'ctrl'), ('off', 'tnt'))}
        r = rh.loc[('orn_pearson', s)]
        L.append(f"| ρ on odor-odor distances, {SUBSET_N[s]} pairs {_subset_words(s)} | {_f(cl['on'])} / {_f(cl['off'])} | "
                 f"{_f(r.rho_on)} / {_f(r.rho_off)} | {_f(r.family_rho_on)} / {_f(r.family_rho_off)} | "
                 f"{_f(r.ceiling_rho_on)} / {_f(r.ceiling_rho_off)} | {_f(r.family_gain_on, '+.2f')} / {_f(r.family_gain_off, '+.2f')} |")
    cond = T['conditions']
    cp = cond[(cond.reading == 'population') & (cond.subset == 'all120')].set_index('family')
    L += ['', '## The split between the conditions (population, 120 pairs)', '',
          '| statistic | family | ceiling family |', '|---|---:|---:|']
    for k, lab in (('rho_on', 'ρ on odor-odor distances, APL on'), ('rho_off', 'ρ on odor-odor distances, APL off'),
                   ('on_from_flies_off_model_effect', "APL on from the flies' APL-off distances + the model's APL effect"),
                   ('on_from_model_off_flies_effect', "APL on from the model's APL-off distances + the flies' APL effect"),
                   ('off_from_flies_on_model_effect', "APL off from the flies' APL-on distances - the model's APL effect"),
                   ('off_from_model_on_flies_effect', "APL off from the model's APL-on distances - the flies' APL effect"),
                   ('resid_on_vs_off', "ρ between the two conditions' rank residuals")):
        L.append(f"| {lab} | {_f(cp.loc['family', k])} | {_f(cp.loc['ceiling', k])} |")
    pp = T['parts'][(T['parts'].reading == 'population') & (T['parts'].family == 'family')].set_index('condition')
    L += ['', '## The odor part and the pair part (the family, population; APL on / APL off)', '',
          '| reading | raw | matched | quantile | best scale |', '|---|---:|---:|---:|---:|']
    for w, lab in (('flies_odor_model_pair', "the flies' odor part + the family's pair part"),
                   ('model_odor_flies_pair', "the family's odor part + the flies' pair part")):
        cells = [f"{_f(pp.loc['on', w + sfx])} / {_f(pp.loc['off', w + sfx])}" for sfx in ('', '_matched', '_quantile', '_best')]
        L.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    L.append(f"| the family's pair part + the Pearson-optimal odor part (in sample, 16 parameters) | "
             f"{_f(pp.loc['on', 'model_pair_optimal_odor'])} / {_f(pp.loc['off', 'model_pair_optimal_odor'])} | | | |")
    L.append(f"| the family | {_f(pp.loc['on', 'model'])} / {_f(pp.loc['off', 'model'])} | | | |")
    st = T['sets'][T['sets'].family == 'family'].set_index(['odors', 'condition'])
    L += ['', '## paa and B-myr (the family, population; APL on / APL off)', '',
          '| reading | paa and B-myr (29 pairs) | the four validation2 odors (54 pairs) |', '|---|---:|---:|']
    for k, lab in (('rho', 'the family, 120 pairs'), ('rho_without', 'the family without those pairs'),
                   ('rho_without_orn', 'the ORN input without those pairs'), ('at_data_raw', "those pairs at the flies' values (raw)"),
                   ('at_data_quantile', "those pairs at the flies' values (quantile)")):
        L.append(f"| {lab} | {_f(st.loc[('paa and B-myr', 'on'), k])} / {_f(st.loc[('paa and B-myr', 'off'), k])} | "
                 f"{_f(st.loc[('the four validation2 odors', 'on'), k])} / {_f(st.loc[('the four validation2 odors', 'off'), k])} |")
    rob = T['robustness'][T['robustness'].subset == 'all120']
    L += ['', '## Robustness (120 pairs; APL on / APL off)', '', '| target | family | ceiling family | ORN input |', '|---|---:|---:|---:|']
    for tname, g in rob.groupby('target', sort=False):
        g = g.set_index('reading')
        L.append(f'| {tname} | ' + ' | '.join(f"{_f(g.loc[r, 'rho_on'])} / {_f(g.loc[r, 'rho_off'])}" for r in ('family', 'ceiling', 'ORN input')) + ' |')
    sr = T['searched_inputs'][T['searched_inputs'].wiring == 'population']
    L += ['', '## The searched inputs on the held-out wirings (population of P3, P4, M2, M3; change against the lab input)', '',
          '| input | ρ, APL on / off | change, APL on / off | change on the 66 megamat pairs | pattern part (lab input) |',
          '|---|---:|---:|---:|---:|']
    for _, r in sr.iterrows():
        L.append(f'| {r.input} | {_f(r.rho_on_120)} / {_f(r.rho_off_120)} | {_f(r.change_on_120, "+.3f")} / {_f(r.change_off_120, "+.3f")} | '
                 f'{_f(r.change_on_66, "+.3f")} / {_f(r.change_off_66, "+.3f")} | {_f(r.pattern, ".3f")} ({_f(r.base_pattern, ".3f")}) |')
    tr = T['tradeoff'][(T['tradeoff'].subset == 'all120')]
    L += ['', '## The trade-off over the archive (rank correlation over distinct variants; APL on / APL off)', '',
          '| against | inside the bands | all |', '|---|---:|---:|']
    for k in TRADE_COLS:
        g = tr[tr.against == k].set_index(['scope', 'condition'])['rank_corr']
        L.append(f"| {k} | {_f(g[('inside the bands', 'on')])} / {_f(g[('inside the bands', 'off')])} | "
                 f"{_f(g[('all', 'on')])} / {_f(g[('all', 'off')])} |")
    re_ = T['realism']
    L += ['', '## Realism of the two families (by_variant.csv; the KC data in baseline-noise units)', '',
          '| reading | response rate, APL on / off | ratio | response rate ρ, APL on / off | odors per KC, APL on / off | '
          'mean distance, APL on / off | flag |', '|---|---:|---:|---:|---:|---:|---|']
    for _, r in re_.iterrows():
        L.append(f'| {r.family}: {r.reading} | {_f(r.rate_on, ".3f")} / {_f(r.rate_off, ".3f")} | {_f(r.ratio)} | '
                 f'{_f(r.get("rate_rho_on"))} / {_f(r.get("rate_rho_off"))} | {_f(r.opk_on_tv)} / {_f(r.opk_off_tv)} | '
                 f'{_f(r.get("mean_dist_on"), ".3f")} / {_f(r.get("mean_dist_off"), ".3f")} | {r.realism_flag if isinstance(r.realism_flag, str) else ""} |')
    L += ['', '## The planning probes against the script', '', '| probe | reading | probe value | script | agree |',
          '|---|---|---:|---:|---|']
    for _, r in probes.iterrows():
        cnt = r.probe in COUNT_PROBES
        L.append(f'| {r.what} | {r.reading} | {_f(r.probe_value, ".0f" if cnt else ".3f")} | '
                 f'{_f(r.script_value, ".0f" if cnt else ".3f")} | {"yes" if r.agree else "no"} |')
    (out / 'comparison_table.md').write_text('\n'.join(L) + '\n')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()
    out = args.out or out_dir()
    out.mkdir(parents=True, exist_ok=True)
    arch = archive()
    T: Dict[str, pd.DataFrame] = {}
    T['archive_summary'] = archive_summary(arch)
    T['tradeoff'] = tradeoff_table(arch)
    T['outcome_c'] = outcome_c(arch)
    fams = [family_tables(n) for n in ('family', 'ceiling')]
    for k in ('conditions', 'parts', 'pairs', 'odors', 'sets'):
        T[k] = pd.concat([f[k] for f in fams], ignore_index=True)
    it = input_table()
    T['input_rhos'], T['input_follow'], T['input_parts'] = it['rhos'], it['follow'], it['parts']
    T['gain91'] = gain_without_paa_bmyr(_readings('family')[1])
    T['searched_inputs'] = searched_inputs_table()
    T['robustness'] = robustness_table()
    T['near_null'] = near_null_changes()
    T['e85_agreement'] = e85_agreement()
    T['yardsticks'] = yardsticks_table()
    T['realism'] = realism_table()
    probes = probes_table(T, arch)
    keep = [c for c in arch.columns if c in ('exp', 'folder', 'tag', 'n_copies', 'copies', 'inside', 'realism_flag', 'rate_on',
            'rate_off', 'rate_rho_on', 'opk_on_tv', 'opk_off_tv') or c in STAT_COLS or '__dist_rho_' in c
            or c.endswith(('_amplitude', '_shape', '_sd_ratio', '_mean_error'))]
    arch[keep].to_csv(out / 'archive.csv', index=False)
    files = dict(archive_summary='archive_summary.csv', tradeoff='tradeoff.csv', outcome_c='outcome_c.csv',
                 conditions='conditions.csv', parts='odor_pair_parts.csv', pairs='residuals_pairs.csv',
                 odors='residuals_odors.csv', sets='odors.csv', input_rhos='input.csv', input_follow='input_follow.csv',
                 input_parts='input_parts.csv', gain91='gain_without_paa_bmyr.csv', searched_inputs='searched_inputs.csv',
                 robustness='robustness.csv', near_null='near_null_changes.csv', e85_agreement='e85_agreement.csv',
                 yardsticks='yardsticks.csv', realism='realism.csv')
    for k, f in files.items():
        T[k].to_csv(out / f, index=False)
    probes.to_csv(out / 'probes.csv', index=False)
    comparison_table(T, probes, out)
    six, two = six_typical_tnt()
    (out / 'provenance.json').write_text(json.dumps(dict(
        generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), script='experiments/E106_rho_gap_parts.py',
        host=subprocess.run(['hostname'], capture_output=True, text=True).stdout.strip(), **H.code_commits(),
        family={n: [list(m) for m in family_members(n)] for n in ('family', 'ceiling')},
        data=dict(on='control flies (TNT3c), cosine, mean over recordings', off='APL output blocked (TNT), cosine, mean over recordings',
                  six_typical_tnt_excluded=['_'.join(map(str, r)) for r in two]),
        held_variants=[[f, p.pattern] for f, p in AE.HELD_VARIANTS], n_boot=N_BOOT, n_splits=N_SPLITS,
        approx_conc_note=H.APPROX_CONC_NOTE, outputs=sorted(list(files.values()) + ['archive.csv', 'probes.csv']),
        inputs=_INPUTS), indent=1, ensure_ascii=False) + '\n')
    print(probes.to_string(index=False, max_colwidth=70), flush=True)
    print(f'-> {out}', flush=True)


if __name__ == '__main__':
    main()
