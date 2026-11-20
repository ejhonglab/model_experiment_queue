#!/usr/bin/env python3
"""E116: an input of paa and B-myr that owes nothing to Yang's KC recordings: searched against Remy's wild-type
validation2 somata and read against Yang's flies, with Hallem's measured paa beside.

The entry is item E116 of experiment_priorities.md (proposed by viki 2026-11-19 16:58; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 16:47; moved into the list 17:00 at viki's user's request; claimed by viki 17:21).
"The family" is C2's (the breadth-2 APL-effect fit, E67) on E72 part B's eight wirings at their own fits (E83's
`saved_fit`; S, searched: real, P1, P2, M1; T, held out: P3, P4, M2, M3). "The two ρ's" are ρ on odor-odor distances
with APL on against the control flies and with APL off against APL output blocked (TNT) (E108's `dist_targets`: all
recordings, and halves A / B of E77's `recording_halves(0)`), on raw counts (cosine distance on spike counts). A thin
script on E109's pieces (its held reading, its validation2 reader, its gate and classes), E107's lock and readers of
Remy's validation2 files (from the copy in data/remy_validation2_kc/), E83's pseudo-odor pool and E38's Hallem reader.

E115's outcome is (c) (the price inside the controls), so step 5's model-only penalty covers the two searched odors only.

Concentrations. paa's input is the lab's recording at 10^-4: what Remy's flies got, and a stand-in for Yang's flies'
10^-3, so every number against Yang's flies involves it; Hallem's paa is at 10^-2 (part H: the as-measured arm is an
upper bound, the profile-only arm carries the class). euc and -aPine are at 10^-1.75 and 10^-1.5 for Yang's flies'
10^-2 and are in every 120-pair statistic (`H.APPROX_CONC_NOTE` in every provenance file).

Parts (CLI `--part`), in the entry's order of work:
- 0 (simulations, about 30 min): the identities (`part0`): E107's as-fitted validation2 runs reproduce on all eight
  wirings (`fitted_on`, 0 differences); the held reading of the lab's input (E85's `held_reading`) reproduces the saved
  fits on all eight (0 differences, both factors 1); the validation2 input's paa and B-myr columns are the lab input's;
  the pseudo-odor runs of paa and B-myr at factor 1 on S reproduce E107's validation2 columns; the objective at the
  lab's input (the incumbents of the search's start) equals E107's ρ (`E107.rho`) of its saved runs against each half
  of Remy's flies; Remy's halves (`remy_halves`) written to halves.json. Writes part0.json ('all_exact').
- H (simulations, about 2 h): Hallem's paa with no search (`hallem_inputs`), the held reading on all eight wirings:
  H_meas (Hallem's paa on the 16 shared glomeruli as measured, the lab's elsewhere), H_prof (the same rescaled to the
  lab's sum over those 16: profile only; carries the class), H_dm6 (the lab's paa with DM6 alone Hallem's), H_perm1 to
  H_perm4 (H_prof's 16 values permuted over the 16 glomeruli, seeds 1 to 4: the controls). Rows in runs_H.csv.
- search (`--search NAME ...`; resumable per generation; one lane of four workers, one per wiring of S, 4 threads each):
  Ra, Rb (halves a, b of Remy's first split), Rc, Rd (her second split); paa then B-myr, each glomerulus's factor in
  [0.5, 2] on the lab's response (72 factors), a (1 + 32) evolution strategy, 12 generations per odor (24 in all), sigma
  from 0.3; candidates simulated as pseudo-odors on S's four wirings at their fits (E83's InputPool, the lab input's
  breadth held) and put into each wiring's as-fitted validation2 run (E107's) in place of the odor's column; the
  objective, maximized: ρ on odor-odor distances of S's population with APL on over Remy's 136 pairs against the
  searched half's mean over flies, minus E83's rate penalty on the searched odor's response fractions (APL on and APL
  off, each wiring against its own fit's, beyond 15 %, 10 per unit). It reads no recording of Yang's. Writes
  search_history_<name>.csv, search_candidates_<name>.csv, input_factors_<name>.csv (E109's layout), inputs/,
  search_runs_<name>.csv, provenance_search_<name>.json.
- read (simulations; resumable): on T, both response counts re-held (E109's held reading): each finished search's input
  applied to the lab's input at 0.6 and at 1 of its log size; its factors permuted over each odor's responding glomeruli
  (E83's `permuted_factors`, seeds 1 to 4) at 0.6 and inverted (-0.6); the combined input (step 11: the best search of
  part R by its own last objective on its searched half, times E109's B2 factors, at 0.6). Rows in runs_R.csv.
- remy (simulations; resumable): step 8, the searched inputs at 0.6 and 1 applied to the validation2 input's paa and
  B-myr on T as fitted (E109's `fitted_on`); rows in runs_V.csv.
- score (no simulations): readings.csv (E108's `unit_scores` per change, unit and readout), remy_readings.csv, the gates
  (gate.json, gate.csv: part H per arm, part R), entries.csv / entries_null.csv (step 10), beside.csv (step 11),
  comparison_table.md.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.

Readings chosen where the entry is open (fixed here, before any run):
- Remy's halves (`remy_halves`): her validation2 flies (E107's lock: per-fly cosine, the lock's drop list; 20 flies),
  split by E77's rule (the flies sorted, `numpy.random.default_rng(seed)`, up to 1000 permutations, the first
  floor(n/2) one half) with the first permutation accepted under which each of the 136 pairs is held by the lock's
  least number of flies per half (`min_per_half`, 1) in each half; split 1 seed 0 (halves a, b), split 2 seed 1 (c, d,
  as E86's second split). A half's matrix is the mean over its flies of each pair (the flies holding it).
- The objective's model matrix is the mean over S's wirings of the cosine distance on spike counts (E106's
  `mean_pairs`, as E107 reads the family); ρ by E107's `rho` over the 136.
- "24 generations": 12 per odor, so that the search has E109's 24 generations in all (the entry's compute estimate).
- The rate penalty's reference is each wiring's fit (its as-fitted validation2 run, whose paa and B-myr columns are the
  fit's); no pattern-part penalty (it would read Yang's recordings).
- Part R's gate (the entry's, read as E109's `gate`): each search of Remy's first split (Ra, Rb) is a direction; T's
  population change from the family at the lab's input against Yang's all recordings and against A and B: a gain of
  0.010 with APL on or 0.015 with APL off on T's population and on at least 3 of the 4 wirings, against each of the
  three; the other ρ not lower by more than 0.005 against any of them; the second split (the mean of Rc and Rd against
  all recordings) of the same sign; a gain against Remy's unseen half: each first-split search's T population (as
  fitted, APL on, 136 pairs) at least 0.010 above the family's against the half it never saw, and the second split's
  mean above 0. Controls: the mean of the 16 permuted (four searches x four seeds) and of the four inverted against all
  recordings each below half the candidate's mean gain (the four searches against all recordings). The other measures
  against all recordings in each direction (E109's thresholds: pattern part (calibrated) +0.003, raw-count loss +0.005,
  mean APL effect -0.01, response rate ρ with APL on and off -0.05, every run at 0.6 of the four searches inside the
  bands; the pending items of 2026-11-17 14:44 as E109). Classes by E109's `class_of`.
- Part H's gate: each arm on all eight wirings against the family at the lab's input: a gain on the eight's population
  and on at least 6 of 8 wirings against all recordings, A and B; the other ρ as above; the controls (the mean of the
  four permutations of H_prof) below half the arm's gain against all recordings; the other measures against all
  recordings. H_prof carries the class; H_meas and H_dm6 are classed beside against the same controls.
- Step 10: the mean log2 factor of each of the 72 entries over the four searches against E109's mean over its six:
  Pearson r over the 72, and the number of E109's eight recurring entries with the same sign; the null permutes the
  four searches' mean factors over the glomeruli within each odor (10000 draws, seed 0); "they agree" when both have
  p at most 0.05 under the null (the share of draws at least as large, the observed counted). Each entry also in the validation2 ORN flies' SDs (E109's `flies_spread`).
- Outcomes (`outcomes`): (e) when H_prof passes; for part R, (a) passes and agree, (b) passes and do not agree, (c) the
  gain holds with APL on only, (d) the gain against Remy's unseen half holds and no gain against Yang's flies, (f)
  nothing carries and the controls gain at least as much.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no connectome table beyond the model's hemibrain wiring and
E72's re-draws of it; nothing of FlyWire, FAFB or BANC; no table of Yang's own Claude's project; no raw-only recording.
Hallem's receptor data (drosolf, as E38 reads it) is not a connectome table. Nothing is published to claude.ai.

Provenance: E109's `commits_clean` (the other checkouts clean, this repository clean in every file a process loads);
every row records `H.code_commits()`; every run's provenance has `H.model_provenance(orn_deltas=...)`.

Outputs: <shared>/results/E116_input_second_flies/; runs under outputs/E116_input_second_flies/ (H, S, R, V), then
data/model_outputs/E116_input_second_flies/.

Usage (from the repository root; one 16-thread lane on viki, under systemd-run --user):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E116_input_second_flies.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E116_input_second_flies.py --part H
    nice -n 10 ./venv/bin/python -u experiments/E116_input_second_flies.py --part search --search Ra Rb Rc Rd
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E116_input_second_flies.py --part read
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E116_input_second_flies.py --part remy
    ./venv/bin/python -u experiments/E116_input_second_flies.py --part score
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E83_searches as SR
from experiments import E109_input_search_rho as E109

EXPERIMENT = 'E116_input_second_flies'
SCRIPT = 'experiments/E116_input_second_flies.py'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
SEARCH_ODORS = ('paa', 'B-myr')
PAA_BMYR = ('paa', 'B-myr')
REMY_SPLITS = {0: ('a', 'b'), 1: ('c', 'd')}
REMY_HALVES = ('a', 'b', 'c', 'd')
REMY_UNSEEN = {'a': 'b', 'b': 'a', 'c': 'd', 'd': 'c'}
SEARCHES = {
    'Ra': dict(half='a', split=0, seed=1161),
    'Rb': dict(half='b', split=0, seed=1162),
    'Rc': dict(half='c', split=1, seed=1163),
    'Rd': dict(half='d', split=1, seed=1164),
}
FIRST_SPLIT, SECOND_SPLIT = ('Ra', 'Rb'), ('Rc', 'Rd')
POP, SIGMA0, BOUNDS = E109.POP, E109.SIGMA0, E109.BOUNDS
GENS = 12                                  # per odor: 24 generations in all
GATE_SIZE, FULL_SIZE, INVERT_SIZE = E109.GATE_SIZE, E109.FULL_SIZE, E109.INVERT_SIZE
PERM_SEEDS = E109.PERM_SEEDS
GAIN, OTHER_SLACK, CONTROL_RATIO = E109.GAIN, E109.OTHER_SLACK, E109.CONTROL_RATIO
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = E109.EQUAL_PATTERN, E109.EQUAL_RAW, E109.EQUAL_MEAN
RATE_RHO_SLACK, BLOCK_SLACK, STRUCTURE_SLACK = E109.RATE_RHO_SLACK, E109.BLOCK_SLACK, E109.STRUCTURE_SLACK
R_GAIN_WIRINGS, H_GAIN_WIRINGS = 3, 6
YANG_TARGETS = ('all', 'A', 'B')
H_ARMS = ('H_prof', 'H_meas', 'H_dm6')
H_CLASS_ARM = 'H_prof'
H_PERM_SEEDS = (1, 2, 3, 4)
H_PERMS = tuple(f'H_perm{s}' for s in H_PERM_SEEDS)
N_NULL = 10000
AGREE_P = 0.05
E109_SIX = E109.MAIN
E109_RECURRING = {'paa': ('DM1', 'DL5', 'VM7d', 'DP1l', 'DP1m'), 'B-myr': ('DM1', 'VC3', 'VM5d')}
COMBINED_WITH = 'B2'                       # E109's search whose factors the best of part R is combined with
CHECKS = (('pattern part (calibrated, all recordings)', 'd_pattern_cal', +1, EQUAL_PATTERN, 'fails'),
          ('raw-count loss (all recordings)', 'd_loss_raw', +1, EQUAL_RAW, 'fails'),
          ('mean APL effect (calibrated)', 'd_mean_effect_cal', -1, EQUAL_MEAN, 'fails'),
          ('response rate ρ, APL on', 'd_rate_rho_on', -1, RATE_RHO_SLACK, 'fails'),
          ('response rate ρ, APL off', 'd_rate_rho_off', -1, RATE_RHO_SLACK, 'fails'),
          ('the block', 'd_block', -1, BLOCK_SLACK, 'pending'),
          ('natural-source alignment, APL on', 'd_vcf_on', -1, STRUCTURE_SLACK, 'pending'),
          ('natural-source alignment, APL off', 'd_vcf_off', -1, STRUCTURE_SLACK, 'pending'),
          ('structure vs wild-type somata', 'd_struct_wt_on', -1, STRUCTURE_SLACK, 'pending'))
PRELOAD = ('experiments.E38_hallem_input', 'experiments.E106_rho_gap_parts', 'experiments.E107_validation2_panel',
           'experiments.E85_orn_fly_population', 'experiments.E108_paa_extrapolation', 'experiments.E64_input_oracle')

_CACHE: Dict[object, object] = {}


# ---------------------------------------------------------------------------
# folders, files, provenance

def out_dir() -> Path:
    """results/E116_input_second_flies/ (E116_OUT_DIR in the environment instead, for smoke tests)."""
    d = Path(os.environ['E116_OUT_DIR']) if os.environ.get('E116_OUT_DIR') else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    """outputs/E116_input_second_flies/<sub> (E116_LOCAL_DIR in the environment instead)."""
    root = Path(os.environ['E116_LOCAL_DIR']) if os.environ.get('E116_LOCAL_DIR') else H.outputs_dir(EXPERIMENT)
    d = root / sub if sub else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _write_json(path: Path, obj) -> None:
    E109._write_json(path, obj)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def commits_clean() -> dict:
    """E109's `commits_clean` (E108's check: the other checkouts clean, this repository clean in every file this process
    loads) after importing this experiment's modules (PRELOAD)."""
    import importlib
    if 'commits' not in _CACHE:
        for m in PRELOAD:
            importlib.import_module(m)
        _CACHE['commits'] = E109.commits_clean()
    return dict(_CACHE['commits'])


def _provenance_add(path: Path, part: str, tag: str, prov: dict) -> None:
    d = _read_json(path) or dict(part=part, script=SCRIPT, points={})
    d['points'][tag] = prov
    _write_json(path, d)


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment: this host's outputs, else the model-output archive; or None."""
    for d in (local_dir(part), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


# ---------------------------------------------------------------------------
# Remy's flies: halves, targets, the objective

def split_flies(flies: Sequence[str], held: pd.DataFrame, seed: int, pairs: Sequence[Tuple[str, str]],
                min_per_half: int = 1, tries: int = 1000) -> Tuple[List[str], List[str]]:
    """E77's rule on Remy's flies: the first of up to `tries` permutations (`default_rng(seed)`, the flies sorted) whose
    first floor(n/2) flies and the rest each hold every pair of `pairs` with `min_per_half` or more flies (`held`: flies
    x pairs, True where a fly holds a pair). Returns (first half, second half), each in sorted order."""
    flies = sorted(map(str, flies))
    hv = held.loc[flies, list(pairs)].astype(bool)
    rng = np.random.default_rng(seed)
    n = len(flies)
    for _ in range(tries):
        perm = [flies[i] for i in rng.permutation(n)]
        a = sorted(perm[:n // 2])
        b = sorted(perm[n // 2:])
        if (hv.loc[a].sum().min() >= min_per_half) and (hv.loc[b].sum().min() >= min_per_half):
            return a, b
    raise RuntimeError(f'no split of {n} flies keeps {min_per_half} fly per half on every pair (seed {seed})')


def _remy_flies() -> Tuple[pd.DataFrame, List[Tuple[str, str]], int]:
    """(her per-fly cosine table, flies x pairs; the 136 pairs as her table orients them; the lock's min_per_half)."""
    if 'flies' not in _CACHE:
        with E109.remy_copy() as E107:
            R = E107.per_fly('cosine')
            lock = E107.require_lock()
            sp = E107.set_pairs(E107.carrying_set())
            if E107.carrying_set() != '136':
                raise RuntimeError("Remy's 136 pairs are not identified (E107's sets)")
        cols = []
        have = set(R.columns)
        for a, b in sp:
            p = (a, b) if (a, b) in have else (b, a)
            if p not in have:
                raise RuntimeError(f'pair {a}, {b} not in her per-fly table')
            cols.append(p)
        _CACHE['flies'] = (R, cols, int(lock['definitions']['halves']['min_per_half']))
    return _CACHE['flies']


def remy_halves() -> Dict[str, List[str]]:
    """{'a', 'b', 'c', 'd': her flies} (`split_flies` on her per-fly cosine table, the 136 pairs; split 1 seed 0, split 2
    seed 1)."""
    if 'halves' not in _CACHE:
        R, cols, mph = _remy_flies()
        out = {}
        for seed, (x, y) in REMY_SPLITS.items():
            out[x], out[y] = split_flies(list(R.index), R.notna(), seed, cols, min_per_half=mph)
        _CACHE['halves'] = out
    return dict(_CACHE['halves'])


def remy_targets() -> Dict[str, pd.Series]:
    """{'a', 'b', 'c', 'd': the half's mean over flies of each pair (the flies holding it), 'all': her mean table (cosine,
    E107's `remy_v2('table_cos')`)} as canonical pair Series."""
    if 'targets' not in _CACHE:
        R, _, _ = _remy_flies()
        out = {h: H.canonical_pairs(R.loc[fl].mean()) for h, fl in remy_halves().items()}
        with E109.remy_copy() as E107:
            out['all'] = H.canonical_pairs(E107.remy_v2('table_cos'))
        _CACHE['targets'] = out
    return dict(_CACHE['targets'])


def remy_pairs(which: str = 'all') -> List[Tuple[str, str]]:
    """Her 136 pairs ('all') or the 31 of paa and B-myr among them ('pb'), E107's orientation."""
    from experiments import E107_validation2_panel as E107
    return E107.subset_of(E107.set_pairs(E107.carrying_set()), which)


def pair_dist(counts: pd.DataFrame) -> pd.Series:
    """Cosine distance on spike counts, canonical pairs (as E107 reads the family)."""
    return H.canonical_pairs(H.pair_distances(counts, 'cosine'))


def population_on(counts: Mapping[str, pd.DataFrame]) -> pd.Series:
    """The mean over wirings of the cosine distance on spike counts (APL on): E106's `mean_pairs` of the wirings."""
    df = pd.concat([pair_dist(c) for c in counts.values()], axis=1)
    return df.dropna().mean(axis=1)


def remy_rho(pop: pd.Series, target: pd.Series, which: str = 'all') -> float:
    """ρ on odor-odor distances (E107's `rho`) of a population against one of Remy's matrices over her 136 pairs (or the
    31 of paa and B-myr)."""
    from experiments import E107_validation2_panel as E107
    return float(E107.rho(pop, target, remy_pairs(which))[0])


def pb_distance(pop: pd.Series) -> float:
    """The distance of paa from B-myr in a pair Series (either orientation)."""
    for k in (PAA_BMYR, PAA_BMYR[::-1]):
        if k in pop.index:
            return float(pop.loc[k])
    return np.nan


class RemyObjective:
    """ρ on odor-odor distances of a population (APL on; `population_on`) over her 136 pairs against one half's mean
    (`remy_targets`)."""

    def __init__(self, half: str):
        self.half = half
        self.data = remy_targets()[half]
        self.pairs = remy_pairs('all')

    def __call__(self, on: Mapping[str, pd.DataFrame]) -> float:
        from experiments import E107_validation2_panel as E107
        return float(E107.rho(population_on(on), self.data, self.pairs)[0])


def score_of(rho_on: float, fracs: Mapping[str, Tuple[float, float]], base_fracs: Mapping[str, Tuple[float, float]]) -> dict:
    """A point's score: ρ (APL on, the searched half) minus E83's rate penalty on the searched odor's response fractions
    (per wiring, APL on and APL off, against the wiring's fit's)."""
    pr = SR.rate_penalty(fracs, base_fracs)
    return dict(rho=float(rho_on), rate_penalty=float(pr), score=float(rho_on - pr))


# ---------------------------------------------------------------------------
# the search (part R)

def v2_fitted(wirings: Sequence[str]) -> Dict[str, Dict[str, tuple]]:
    """E107's saved as-fitted validation2 runs ({wiring: {'on', 'off': (responses, spike counts)}}, columns abbreviated)."""
    from experiments import E107_validation2_panel as E107
    return E107.family_runs('validation2', 'fitted', list(wirings))


def _fracs(inc: Mapping[str, Mapping[str, tuple]], odor: str) -> Dict[str, Tuple[float, float]]:
    oc = SR._odor_col(next(iter(inc.values()))['on'][1], odor)
    return {w: (float(inc[w]['on'][0][oc].astype(float).gt(0.5).mean()),
                float(inc[w]['off'][0][oc].astype(float).gt(0.5).mean())) for w in inc}


def candidate_scores(inc: Mapping[str, Mapping[str, tuple]], odor: str, res: Mapping[str, Tuple[np.ndarray, ...]],
                     obj: RemyObjective, base_fracs: Mapping[str, Tuple[float, float]]) -> List[dict]:
    """Each candidate's score: the odor's column of every wiring's incumbent replaced by the candidate's (column i of
    `res[w]`, InputPool's layout: responses on, counts on, responses off, counts off)."""
    oc = SR._odor_col(next(iter(inc.values()))['on'][1], odor)
    n_cand = next(iter(res.values()))[1].shape[1]
    rows = []
    for i in range(n_cand):
        on = {w: SR._with_column(inc[w]['on'][1], oc, res[w][1][:, i]) for w in inc}
        fr = {w: (float((res[w][0][:, i] > 0.5).mean()), float((res[w][2][:, i] > 0.5).mean())) for w in inc}
        rows.append(dict(score_of(obj(on), fr, base_fracs), fracs=fr))
    return rows


def _sdir(name: str) -> Path:
    return local_dir(f'S/state_{name}')


def factors_path(name: str) -> Path:
    return out_dir() / f'input_factors_{name}.csv'


def write_factors(name: str, done: Mapping[str, Sequence[float]]) -> pd.DataFrame:
    """input_factors_<name>.csv (E109's layout: search, half, odor, glomerulus, lab_delta, factor, log2_factor,
    searched_delta) and inputs/searched_<name>.csv (the lab's input with the factors, E38's layout)."""
    lab = H.apl16_orn_deltas()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    now = lab.astype(float).copy()
    rows = []
    for odor in SEARCH_ODORS:
        if odor not in done:
            continue
        f = np.asarray(done[odor], dtype=float)
        j = ab.index(odor)
        now.iloc[:, j] = lab.iloc[:, j].to_numpy(dtype=float) * f
        for gi, g in enumerate(lab.index):
            rows.append(dict(search=name, half=SEARCHES[name]['half'], odor=odor, glomerulus=g,
                             lab_delta=float(lab.iloc[gi, j]), factor=float(f[gi]), log2_factor=float(np.log2(f[gi])),
                             searched_delta=float(lab.iloc[gi, j] * f[gi])))
    fac = pd.DataFrame(rows)
    fac.to_csv(factors_path(name), index=False)
    E83._write_input(out_dir(), f'searched_{name}', now)
    return fac


def search_finished(name: str) -> bool:
    p = factors_path(name)
    return p.exists() and set(pd.read_csv(p)['odor']) == set(SEARCH_ODORS)


def search_factors(name: str) -> pd.DataFrame:
    """A finished search's factors, glomeruli x odors (the lab input's layout, 1 where not searched)."""
    if not search_finished(name):
        raise RuntimeError(f'search {name} is not finished ({factors_path(name)})')
    return E83.read_factors(factors_path(name))


def part0_passed() -> bool:
    return bool(_read_json(out_dir() / 'part0.json').get('all_exact'))


def search(name: str, resume: bool = True, threads: int = E109.WORKER_THREADS, gens: int = GENS,
           pop: int = POP) -> pd.DataFrame:
    """Search `name` (SEARCHES; the module docstring): the pseudo-odor path checked at factor 1, then per odor the (1 +
    `pop`) evolution strategy for `gens` generations on S's four wirings, maximizing `score_of` against the searched half
    of Remy's flies; checkpointed every generation, resumed by default. Returns the history."""
    cfg = SEARCHES[name]
    half, unseen = cfg['half'], REMY_UNSEEN[cfg['half']]
    if not part0_passed():
        raise RuntimeError('run --part 0 first (part0.json: the identities must hold before any search)')
    cc = commits_clean()
    out = out_dir()
    sdir = _sdir(name)
    ck = sdir / 'checkpoint.json'
    hist_path, cand_path = out / f'search_history_{name}.csv', out / f'search_candidates_{name}.csv'
    lab = H.apl16_orn_deltas()
    n_gl = len(lab.index)
    saved = v2_fitted(S_WIRINGS)
    obj, other, allr = RemyObjective(half), RemyObjective(unseen), RemyObjective('all')
    lo, hi = np.log(BOUNDS[0]), np.log(BOUNDS[1])
    t_start = time.time()
    pool = SR.InputPool(S_WIRINGS, threads)
    try:
        for odor in SEARCH_ODORS:
            res = pool.run(odor, np.ones((1, n_gl)))
            oc = SR._odor_col(saved['real']['on'][1], odor)
            for w in S_WIRINGS:
                for c, kr, kc in (('on', 0, 1), ('off', 2, 3)):
                    if not (np.array_equal(res[w][kc][:, 0].astype(float), saved[w][c][1][oc].to_numpy(dtype=float))
                            and np.array_equal(res[w][kr][:, 0].astype(float), saved[w][c][0][oc].to_numpy(dtype=float))):
                        raise RuntimeError(f'{w}: the pseudo-odor run of {odor} at factor 1 differs from the {c} validation2 run')
        print(f"{name}: the pseudo-odor runs of paa and B-myr at factor 1 reproduce E107's four validation2 runs", flush=True)
        if resume and ck.exists():
            st = json.loads(ck.read_text())
            done = {o: list(v) for o, v in st['done'].items()}
            rng = SR._rng_from(st['rng'])
            oi, g = int(st['odor_index']), int(st['generation'])
            frames = st['frames']
            inc = {w: {c: (SR._load_frame(sdir, f'{w}__{frames}', f'{c}_responses'),
                           SR._load_frame(sdir, f'{w}__{frames}', f'{c}_spike_counts')) for c in ('on', 'off')}
                   for w in S_WIRINGS}
            base_fracs = {o: {w: tuple(v) for w, v in d.items()} for o, d in st['base_fracs'].items()}
            parent, parent_score, sigma = None, None, SIGMA0
            if g + 1 < gens:
                start_odor, start_gen = oi, g + 1
                parent, parent_score, sigma = np.array(st['parent']), float(st['parent_score']), float(st['sigma'])
                now = score_of(obj({w: inc[w]['on'][1] for w in S_WIRINGS}), _fracs(inc, SEARCH_ODORS[oi]),
                               base_fracs[SEARCH_ODORS[oi]])['score']
                if abs(now - parent_score) > 1e-12:
                    raise RuntimeError(f'{name}: the checkpoint parent scores {now} now, {parent_score} when saved')
            else:
                start_odor, start_gen = oi + 1, 0
            pos = oi * gens + g
            order = {o: k for k, o in enumerate(SEARCH_ODORS)}
            SR._trim_csv(hist_path, lambda d: (d.odor == 'start') | (d.odor.map(order) * gens + d.generation <= pos))
            SR._trim_csv(cand_path, lambda d: d.odor.map(order) * gens + d.generation <= pos)
            print(f'{name}: resuming at odor {start_odor} generation {start_gen}', flush=True)
        else:
            for p in (hist_path, cand_path, factors_path(name)):
                if p.exists():
                    p.unlink()
            SR._reset_state(sdir)
            inc = {w: {c: tuple(x.copy() for x in saved[w][c]) for c in ('on', 'off')} for w in S_WIRINGS}
            base_fracs = {o: _fracs(inc, o) for o in SEARCH_ODORS}
            done, rng = {}, np.random.default_rng(cfg['seed'])
            start_odor, start_gen, parent, parent_score, sigma = 0, 0, None, None, SIGMA0
            frames = 'start'
            for w in S_WIRINGS:
                SR._save_frames(sdir, f'{w}__{frames}', {'on_responses': inc[w]['on'][0], 'on_spike_counts': inc[w]['on'][1],
                                                         'off_responses': inc[w]['off'][0], 'off_spike_counts': inc[w]['off'][1]})
            on0 = {w: inc[w]['on'][1] for w in S_WIRINGS}
            SR._append_csv(hist_path, [dict(search=name, half=half, odor='start', generation=-1, rho_half=obj(on0),
                                            rho_unseen=other(on0), rho_all=allr(on0), seconds=time.time() - t_start)])
            print(f'{name}: start, S against half {half} {obj(on0):.4f}, against {unseen} {other(on0):.4f}', flush=True)
        old = _read_json(out / f'provenance_search_{name}.json')
        invocations = list(old.get('invocations', [])) + [dict(generated=_now(), resumed=bool(resume and old), **cc)]
        prov = dict(generated=_now(), script=SCRIPT, search=name, **cfg, unseen=unseen, **cc, invocations=invocations,
                    model=H.model_provenance(orn_deltas=lab, _wPNKC=E83.wiring_parts('real')['W']),
                    base_fracs=base_fracs, halves={h: remy_halves()[h] for h in (half, unseen)},
                    incumbents="E107's as-fitted validation2 runs (validation2_fitted_<wiring>)",
                    wirings={w: dict(folder=E83.saved_fit(w)['folder'], tag=E83.saved_fit(w)['tag'],
                                     params=E83.saved_fit(w)['p']) for w in S_WIRINGS},
                    pop=pop, generations_per_odor=gens, bounds=BOUNDS, sigma0=SIGMA0, sigma_min=SR.SIGMA_MIN,
                    sigma_max=SR.SIGMA_MAX, success_up=SR.SUCCESS_UP, failure_down=SR.FAILURE_DOWN,
                    rate_tol=SR.I_RATE_TOL, penalty=SR.I_PENALTY, workers=len(S_WIRINGS), threads_per_worker=threads,
                    objective=("ρ on odor-odor distances of S's population with APL on (cosine on spike counts, the "
                               "validation2 panel as fitted) over Remy's 136 pairs against the searched half's mean over "
                               "flies, minus E83's rate penalty on the searched odor's response fractions"),
                    reads_yang=False, approx_conc_note=H.APPROX_CONC_NOTE)
        _write_json(out / f'provenance_search_{name}.json', prov)
        for oi in range(start_odor, len(SEARCH_ODORS)):
            odor = SEARCH_ODORS[oi]
            oc = SR._odor_col(inc['real']['on'][1], odor)
            if oi == start_odor and parent is not None:
                g_first = start_gen
            else:
                parent = np.ones(n_gl)
                parent_score = score_of(obj({w: inc[w]['on'][1] for w in S_WIRINGS}), _fracs(inc, odor),
                                        base_fracs[odor])['score']
                sigma, g_first = SIGMA0, 0
            for g in range(g_first, gens):
                t0 = time.time()
                offspring = np.exp(SR.es_offspring(rng, np.log(parent), sigma, pop, lo, hi))
                res = pool.run(odor, offspring)
                rows = candidate_scores(inc, odor, res, obj, base_fracs[odor])
                scores = np.array([r['score'] for r in rows])
                b = int(np.argmax(scores))
                success = bool(scores[b] > parent_score)
                if success:
                    parent, parent_score = offspring[b].copy(), float(scores[b])
                    for w in S_WIRINGS:
                        r_on, s_on = inc[w]['on']
                        r_off, s_off = inc[w]['off']
                        inc[w] = {'on': (SR._with_column(r_on, oc, res[w][0][:, b]), SR._with_column(s_on, oc, res[w][1][:, b])),
                                  'off': (SR._with_column(r_off, oc, res[w][2][:, b]), SR._with_column(s_off, oc, res[w][3][:, b]))}
                sigma = SR.es_sigma(sigma, success)
                on = {w: inc[w]['on'][1] for w in S_WIRINGS}
                rh, ru, ra = obj(on), other(on), allr(on)
                fr = _fracs(inc, odor)
                hrow = dict(search=name, half=half, odor=odor, generation=g, best_score=float(scores[b]),
                            parent_score=parent_score, success=success, sigma=sigma, best_rho=rows[b]['rho'],
                            best_rate_penalty=rows[b]['rate_penalty'], rho_half=rh, rho_unseen=ru, rho_all=ra,
                            **{f'frac_on_{w}': fr[w][0] for w in S_WIRINGS}, **{f'frac_off_{w}': fr[w][1] for w in S_WIRINGS},
                            **{f'base_frac_on_{w}': base_fracs[odor][w][0] for w in S_WIRINGS},
                            **{f'base_frac_off_{w}': base_fracs[odor][w][1] for w in S_WIRINGS},
                            log2_factor_rms=float(np.sqrt(np.mean(np.log2(parent) ** 2))),
                            n_at_bound=int((np.abs(np.abs(np.log(parent)) - hi) < 1e-12).sum()), seconds=time.time() - t0)
                SR._append_csv(cand_path, [dict(search=name, half=half, odor=odor, generation=g, candidate=i, score=r['score'],
                                                rho=r['rho'], rate_penalty=r['rate_penalty']) for i, r in enumerate(rows)])
                if success:
                    frames = f'o{oi}g{g}'
                    for w in S_WIRINGS:
                        SR._save_frames(sdir, f'{w}__{frames}', {'on_responses': inc[w]['on'][0], 'on_spike_counts': inc[w]['on'][1],
                                                                 'off_responses': inc[w]['off'][0], 'off_spike_counts': inc[w]['off'][1]})
                if g == gens - 1:
                    done[odor] = parent.tolist()
                SR._write_json(ck, dict(search=name, odor_index=oi, generation=g, frames=frames, parent=parent.tolist(),
                                        parent_score=parent_score, sigma=sigma, rng=SR._rng_state(rng), done=done,
                                        base_fracs={o: {w: list(v) for w, v in d.items()} for o, d in base_fracs.items()}))
                SR._drop_other_frames(sdir, frames)
                SR._append_csv(hist_path, [hrow])
                print(f"{name} {odor} gen {g}: best {scores[b]:.4f} (parent {parent_score:.4f}); S against half {half} "
                      f"{rh:.4f}, against {unseen} {ru:.4f}, all {ra:.4f}; sigma {sigma:.2f}; {time.time() - t0:.0f}s", flush=True)
            parent = None
            write_factors(name, done)
    finally:
        pool.close()
    rows = []
    for w in S_WIRINGS:
        o = {c: inc[w][c] for c in ('on', 'off')}
        tag = f'{name}_{w}'
        save_run('S', tag, o)
        P = E83.saved_fit(w)
        rows.append(dict(search=name, half=half, wiring=w, tag=tag, panel='validation2', reading='fitted',
                         fit_folder=P['folder'], fit_tag=P['tag'], **P['p'],
                         rate_on=float(o['on'][0].to_numpy(dtype=float).mean()),
                         rate_off=float(o['off'][0].to_numpy(dtype=float).mean()), **cc))
    pd.DataFrame(rows).to_csv(out / f'search_runs_{name}.csv', index=False)
    return pd.read_csv(hist_path)


# ---------------------------------------------------------------------------
# inputs

def hallem_inputs(lab: Optional[pd.DataFrame] = None, hal: Optional[pd.Series] = None
                  ) -> Dict[str, Tuple[pd.DataFrame, dict]]:
    """Part H's inputs, {arm: (input, description)}, from the lab's input (`lab`, `H.apl16_orn_deltas()` by default) and
    Hallem's paa (`hal`: E38's `hallem_deltas`, the 16 glomeruli the model and Hallem's receptors share, 10^-2):
    H_meas (Hallem's on the 16, the lab's elsewhere), H_prof (the same times the lab's sum over the 16 / Hallem's sum),
    H_dm6 (the lab's with DM6 alone Hallem's), H_perm<s> (H_prof's 16 values permuted over the 16, `default_rng(s)`)."""
    if lab is None:
        lab = H.apl16_orn_deltas()
    if hal is None:
        from experiments import E38_hallem_input as E38
        hal = E38.hallem_deltas(lab.index)['paa']
    hal = hal.astype(float)
    col = next(c for c in lab.columns if H.odor_abbrev(c) == 'paa')
    shared = list(hal.index)
    lab_sum = float(lab.loc[shared, col].astype(float).sum())
    scale = lab_sum / float(hal.sum())
    prof = hal * scale

    def with_paa(values: pd.Series) -> pd.DataFrame:
        orn = lab.astype(float).copy()
        orn.loc[list(values.index), col] = values.to_numpy(dtype=float)
        return orn
    base = dict(shared_glomeruli=shared, lab_sum=lab_sum, hallem_sum=float(hal.sum()), profile_scale=scale,
                hallem_conc_log10=-2.0)
    out = {'H_meas': (with_paa(hal), dict(base, how="Hallem's paa on the 16 shared glomeruli as measured")),
           'H_prof': (with_paa(prof), dict(base, how="Hallem's paa on the 16, rescaled to the lab's sum over them")),
           'H_dm6': (with_paa(hal.loc[['DM6']]), dict(base, how="the lab's paa with DM6 alone Hallem's"))}
    for s in H_PERM_SEEDS:
        rng = np.random.default_rng(s)
        v = prof.to_numpy()[rng.permutation(len(prof))]
        out[f'H_perm{s}'] = (with_paa(pd.Series(v, index=prof.index)),
                             dict(base, how=f"H_prof's 16 values permuted over the 16 glomeruli (seed {s})", seed=s))
    return out


def held(label: str, orn: pd.DataFrame, name: str) -> dict:
    """A held reading (E85's `held_reading` with this experiment's breadth cache)."""
    from experiments import E85_orn_fly_population as E85
    return E85.held_reading(label, orn, E83.input_tag_for(orn, f'E116_{name}'), cache_dir=local_dir())


def v2_with(F: pd.DataFrame, size: float) -> pd.DataFrame:
    """The validation2 ORN input (E107's `panel_input`) with paa's and B-myr's columns times F ** size."""
    from experiments import E107_validation2_panel as E107
    v2 = E107.panel_input('validation2').astype(float).copy()
    for c in v2.columns:
        o = H.odor_abbrev(c)
        if o in PAA_BMYR:
            v2[c] = v2[c].to_numpy(dtype=float) * np.power(F[o].reindex(v2.index).to_numpy(dtype=float), float(size))
    return v2


# ---------------------------------------------------------------------------
# part 0

def part0(threads: int = E109.WORKER_THREADS) -> dict:
    """Part 0 (the module docstring); part0.json, halves.json, remy_reference.csv. Raises if an identity fails."""
    from experiments import E106_rho_gap_parts as E106
    from experiments import E107_validation2_panel as E107
    cc = commits_clean()
    out = out_dir()
    res = dict(script=SCRIPT, generated=_now(), **cc, columns={}, v2_path={}, held_path={}, search_path={}, objective={})
    lab = H.apl16_orn_deltas()
    v2 = E107.panel_input('validation2')
    for o in PAA_BMYR:
        cl = next(c for c in lab.columns if H.odor_abbrev(c) == o)
        cv = next(c for c in v2.columns if H.odor_abbrev(c) == o)
        res['columns'][o] = bool(cl == cv and np.array_equal(lab[cl].to_numpy(float), v2[cv].to_numpy(float)))
    # Remy's halves
    hv = remy_halves()
    R, cols, mph = _remy_flies()
    _write_json(out / 'halves.json', dict(halves=hv, min_per_half=mph, n_flies=int(len(R)), splits=REMY_SPLITS,
                                          held_per_half={h: int(R.loc[fl, cols].notna().sum().min()) for h, fl in hv.items()}))
    print(f'E116 0 halves: ' + json.dumps({h: len(v) for h, v in hv.items()}), flush=True)
    # E107's as-fitted validation2 runs on all eight
    runs = v2_fitted(ALL_WIRINGS)
    for w in ALL_WIRINGS:
        o = E109.fitted_on(w, v2)
        ns, nr = E83._same({'on': o['on'], 'off': o['off']}, runs[w])
        res['v2_path'][w] = dict(count_diffs=ns, response_diffs=nr)
        print(f"E116 0 validation2 path {w}: {res['v2_path'][w]}", flush=True)
        del o
        gc.collect()
    # the held reading of the lab's input on all eight
    for w in ALL_WIRINGS:
        hr = held(w, lab, 'I0')
        ns, nr = E83._same({'on': hr['on'], 'off': hr['off']}, E83.saved_fit(w)['saved'])
        res['held_path'][w] = dict(count_diffs=ns, response_diffs=nr, thr_factor=hr['thr_factor'], apl_factor=hr['apl_factor'])
        print(f"E116 0 held path {w}: {res['held_path'][w]}", flush=True)
        del hr
        gc.collect()
    # the search's path: the pseudo-odor runs at factor 1 against E107's validation2 columns
    pool = SR.InputPool(S_WIRINGS, threads)
    try:
        for odor in PAA_BMYR:
            r = pool.run(odor, np.ones((1, len(lab.index))))
            oc = SR._odor_col(runs['real']['on'][1], odor)
            diffs = {}
            for w in S_WIRINGS:
                n = 0
                for c, kr, kc in (('on', 0, 1), ('off', 2, 3)):
                    n += int((r[w][kc][:, 0].astype(float) != runs[w][c][1][oc].to_numpy(dtype=float)).sum())
                    n += int((r[w][kr][:, 0].astype(float) != runs[w][c][0][oc].to_numpy(dtype=float)).sum())
                diffs[w] = n
            res['search_path'][odor] = diffs
            print(f'E116 0 search path {odor}: differences {diffs}', flush=True)
    finally:
        pool.close()
    # the objective at the lab's input against E107's ρ of its saved runs (computed by E107's own readers)
    per = {w: {c: H.canonical_pairs(H.pair_distances(runs[w][c][1], 'cosine')) for c in ('on', 'off')} for w in S_WIRINGS}
    pop_e107 = E106.mean_pairs(list(per.values()))['on']
    sp = remy_pairs('all')
    for h in REMY_HALVES:
        v = RemyObjective(h)({w: runs[w]['on'][1] for w in S_WIRINGS})
        ref = float(E107.rho(pop_e107, E107.fly_mean(R.loc[hv[h]], min_flies=1), sp)[0])
        res['objective'][h] = dict(objective=v, e107=ref, ok=bool(abs(v - ref) <= 1e-12))
    print('E116 0 objective: ' + json.dumps(res['objective']), flush=True)
    # the family's reference against Remy's matrices (E107's runs, as fitted): all eight, S, T and each wiring
    rows = []
    tg = remy_targets()
    for unit, ws in (('all8', ALL_WIRINGS), ('S', S_WIRINGS), ('T', T_WIRINGS)) + tuple((w, (w,)) for w in ALL_WIRINGS):
        pop = population_on({w: runs[w]['on'][1] for w in ws})
        for t in ('all',) + REMY_HALVES:
            rows.append(dict(change='I0', unit=unit, target=t, rho=remy_rho(pop, tg[t]), rho_pb=remy_rho(pop, tg[t], 'pb'),
                             paa_bmyr=pb_distance(pop)))
    pd.DataFrame(rows).to_csv(out / 'remy_reference.csv', index=False)
    ok = (all(res['columns'].values())
          and all(d['count_diffs'] == 0 and d['response_diffs'] == 0 for d in res['v2_path'].values())
          and all(d['count_diffs'] == 0 and d['response_diffs'] == 0 and d['thr_factor'] == 1.0 and d['apl_factor'] == 1.0
                  for d in res['held_path'].values())
          and all(n == 0 for d in res['search_path'].values() for n in d.values())
          and all(d['ok'] for d in res['objective'].values()))
    res['all_exact'] = bool(ok)
    _write_json(out / 'part0.json', res)
    print(json.dumps(dict(all_exact=ok), indent=1), flush=True)
    if not ok:
        raise RuntimeError('E116 part 0: an identity does not hold (part0.json); nothing else runs')
    return res


# ---------------------------------------------------------------------------
# the held readings (parts H and read) and Remy's panel (part remy)

def _held_points(rp: Path, pp: Path, part: str, plan: List[dict], inputs: Mapping[str, pd.DataFrame],
                 max_points: Optional[int] = None) -> pd.DataFrame:
    """Runs the held readings of `plan` ([dict(change, wiring, tag, **meta)]) not yet in `rp`, wiring by wiring, each
    wiring's simulator first checked against its saved fit."""
    from experiments.E77_input_screen import frame_sha256
    cc = commits_clean()
    out = out_dir()
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    todo = [p for p in plan if p['tag'] not in done]
    if max_points is not None:
        todo = todo[:max_points]
    print(f'E116 {part}: {len(todo)} points to run', flush=True)
    for w in [w for w in ALL_WIRINGS if any(p['wiring'] == w for p in todo)]:
        sim, info = E83.build_simulator(w)
        E83._check_saved(w, sim, info['params'], E83.saved_fit(w)['saved'])
        del sim, info
        gc.collect()
        for pt in [p for p in todo if p['wiring'] == w]:
            orn = inputs[pt['change']]
            path = E83._write_input(out, pt['change'], orn)
            res = held(w, orn, pt['change'])
            o = {'on': res['on'], 'off': res['off']}
            save_run(part, pt['tag'], o)
            meta = {k: v for k, v in pt.items() if k not in ('change', 'wiring', 'tag')}
            row = E83._run_row(part, pt['change'], w, pt['tag'], res, o, **meta, input_tag=res['input_tag'],
                               orn_input=str(path.relative_to(H.shared_dir())) if H.shared_dir() in path.parents else str(path),
                               input_frame_sha256=frame_sha256(orn),
                               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}))
            row.update(E109.run_extras(o))
            row.update({k: v for k, v in cc.items() if k == 'queue_dirty_paths'})
            E83._append_row(rp, row)
            _provenance_add(pp, part, pt['tag'], H.model_provenance(orn_deltas=orn))
            print(f"E116 {part} {pt['tag']}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/"
                  f"{res['target_on']} on, {res['n_off']}/{res['target_off']} off); flag '{row.get('realism_flag') or ''}'; "
                  f"{res['seconds']:.0f}s", flush=True)
            del res, o
            gc.collect()
    return E83._read_rows(rp)


def part_h(max_points: Optional[int] = None) -> pd.DataFrame:
    """Part H: the held readings of Hallem's paa (`hallem_inputs`) on all eight wirings; the arms first, then the
    permutations. Rows in runs_H.csv; the inputs' descriptions in hallem_inputs.json."""
    if not part0_passed():
        raise RuntimeError('run --part 0 first')
    hi = hallem_inputs()
    _write_json(out_dir() / 'hallem_inputs.json', {k: v[1] for k, v in hi.items()})
    order = list(H_ARMS) + list(H_PERMS)
    plan = [dict(change=ch, wiring=w, tag=f'{ch}_{w}', kind='arm' if ch in H_ARMS else 'perm') for w in ALL_WIRINGS
            for ch in order]
    return _held_points(out_dir() / 'runs_H.csv', out_dir() / 'provenance_H.json', 'H', plan,
                        {k: v[0] for k, v in hi.items()}, max_points)


def best_search() -> Optional[str]:
    """The search of part R with the highest last objective on its own searched half (its history's parent score; no
    reading of Yang's), once all four are finished."""
    if not all(search_finished(n) for n in SEARCHES):
        return None
    sc = {n: float(pd.read_csv(out_dir() / f'search_history_{n}.csv').query("odor != 'start'")['parent_score'].iloc[-1])
          for n in SEARCHES}
    return max(sc, key=sc.get)


def change_meta(names: Optional[Sequence[str]] = None) -> Dict[str, dict]:
    """{change: meta} of the held readings on T: each finished search at GATE_SIZE and FULL_SIZE, its permutations at
    GATE_SIZE and its inversion; the combined input once all four are finished."""
    out = {}
    for n in (names or SEARCHES):
        if not search_finished(n):
            continue
        m = dict(search=n, half=SEARCHES[n]['half'], split=SEARCHES[n]['split'])
        out[f'{n}_s{GATE_SIZE:g}'] = dict(m, kind='search', size=GATE_SIZE, seed=None)
        out[f'{n}_s{FULL_SIZE:g}'] = dict(m, kind='search', size=FULL_SIZE, seed=None)
        for s in PERM_SEEDS:
            out[f'{n}_perm{s}'] = dict(m, kind='perm', size=GATE_SIZE, seed=s)
        out[f'{n}_inv'] = dict(m, kind='inverted', size=INVERT_SIZE, seed=None)
    b = best_search()
    if b is not None and names is None:
        out[f'comb_{b}_{COMBINED_WITH}'] = dict(search=b, half=SEARCHES[b]['half'], split=SEARCHES[b]['split'],
                                                kind='combined', size=GATE_SIZE, seed=None)
    return out


def change_input(change: str, m: Mapping) -> pd.DataFrame:
    """A reading's input on the lab's panel: the lab's input times the search's factors (permuted for 'perm'; times
    E109's B2 factors for 'combined') at the size."""
    F = search_factors(m['search'])
    if m['kind'] == 'perm':
        F = E83.permuted_factors(F, int(m['seed']))
    elif m['kind'] == 'combined':
        F = F * E109.search_factors(COMBINED_WITH).reindex(index=F.index, columns=F.columns).fillna(1.0)
    return E83.input_at(F, m['size'])


def read_plan(names: Optional[Sequence[str]] = None) -> List[dict]:
    meta = change_meta(names)
    rank = {'search': 0, 'perm': 1, 'inverted': 1, 'combined': 3}
    order = sorted(meta, key=lambda ch: (2 if meta[ch]['kind'] == 'search' and meta[ch]['size'] == FULL_SIZE
                                         else rank[meta[ch]['kind']], list(meta).index(ch)))
    return [dict(change=ch, wiring=w, tag=f'{ch}_{w}', **meta[ch]) for w in T_WIRINGS for ch in order]


def part_read(names: Optional[Sequence[str]] = None, max_points: Optional[int] = None) -> pd.DataFrame:
    """The held readings on T (resumable): rows in runs_R.csv."""
    if not part0_passed():
        raise RuntimeError('run --part 0 first')
    plan = read_plan(names)
    meta = change_meta(names)
    inputs = {ch: change_input(ch, meta[ch]) for ch in meta}
    return _held_points(out_dir() / 'runs_R.csv', out_dir() / 'provenance_R.json', 'R', plan, inputs, max_points)


def part_remy(names: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 8 (resumable): each finished search's factors at GATE_SIZE and FULL_SIZE on the validation2 input's paa and
    B-myr (`v2_with`), the four wirings of T as fitted (E109's `fitted_on`); rows in runs_V.csv."""
    from experiments.E77_input_screen import frame_sha256
    if not part0_passed():
        raise RuntimeError('run --part 0 first')
    cc = commits_clean()
    out = out_dir()
    rp, pp = out / 'runs_V.csv', out / 'provenance_V.json'
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    for n in [n for n in (names or SEARCHES) if search_finished(n)]:
        for size in (GATE_SIZE, FULL_SIZE):
            orn = v2_with(search_factors(n), size)
            ch = f'{n}_s{size:g}'
            for w in T_WIRINGS:
                tag = f'V_{ch}_{w}'
                if tag in done:
                    continue
                t0 = time.time()
                o = E109.fitted_on(w, orn)
                save_run('V', tag, {'on': o['on'], 'off': o['off']})
                row = dict(part='V', change=ch, search=n, size=size, wiring=w, tag=tag, panel='validation2',
                           reading='fitted', input_frame_sha256=frame_sha256(orn), fit_folder=E83.saved_fit(w)['folder'],
                           fit_tag=E83.saved_fit(w)['tag'], **{k: o['params'][k] for k in E83.PARAM_KEYS},
                           rate_on=float(o['on'][0].to_numpy(dtype=float).mean()),
                           rate_off=float(o['off'][0].to_numpy(dtype=float).mean()), seconds=time.time() - t0, **cc)
                E83._append_row(rp, row)
                _provenance_add(pp, 'V', tag, H.model_provenance(orn_deltas=orn))
                print(f'E116 V {tag}: {time.time() - t0:.0f}s', flush=True)
    return E83._read_rows(rp)


# ---------------------------------------------------------------------------
# scoring (no simulations)

def runs_table() -> pd.DataFrame:
    """Every reading run on Yang's panel: the saved fits as 'I0' (all eight), runs_H.csv's and runs_R.csv's, with
    'flagged'."""
    rows = [dict(change='I0', wiring=w, tag=f'I0_{w}', kind='reference', **E83.realism(E83.saved_fit(w)['saved']),
                 **E109.run_extras(E83.saved_fit(w)['saved'])) for w in ALL_WIRINGS]
    parts = [pd.DataFrame(rows)] + [E83._read_rows(out_dir() / f) for f in ('runs_H.csv', 'runs_R.csv')]
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates('tag', keep='last').reset_index(drop=True)
    df['flagged'] = [E109._flagged(x) for x in df['realism_flag']]
    return df


def run_outputs(change: str, wiring: str) -> Dict[str, tuple]:
    if change == 'I0':
        return E83.saved_fit(wiring)['saved']
    for part in ('H', 'R'):
        o = load_run(part, f'{change}_{wiring}')
        if o is not None:
            return o
    raise RuntimeError(f'no saved run {change}_{wiring}')


def readings(runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings.csv: per change, unit (each wiring; 'all8' for a change on all eight, 'T' for one on T or on all eight)
    and readout ('raw', 'cal'), E108's `unit_scores`, the mean realism quantities over the unit's runs, the runs flagged
    and the distance of paa from B-myr."""
    from experiments import E108_paa_extrapolation as E108
    runs = runs_table() if runs is None else runs
    rows = []
    for ch, g in runs.groupby('change', sort=False):
        ws = set(g['wiring'])
        if ws == set(ALL_WIRINGS):
            wirings, pops = ALL_WIRINGS, (('all8', ALL_WIRINGS), ('T', T_WIRINGS))
        elif ws == set(T_WIRINGS):
            wirings, pops = T_WIRINGS, (('T', T_WIRINGS),)
        else:
            continue
        meta = {r.wiring: r for r in g.itertuples()}
        for ro in ('raw', 'cal'):
            ds = {w: E83.pair_distances(run_outputs(ch, w), ro) for w in wirings}
            units = [(w, (w,)) for w in wirings] + list(pops)
            for unit, uw in units:
                d = ds[unit] if unit in ds else E83.population([ds[w] for w in uw])
                sub = g[g['wiring'].isin(uw)]
                pb = E83._cs(d['on']).index.get_loc(tuple(sorted(PAA_BMYR)))
                rows.append(dict(change=ch, unit=unit, kind='wiring' if unit in ds else 'population', readout=ro,
                                 n_flagged=int(sum(meta[w].flagged for w in uw)),
                                 **{f'run_{k}': float(sub[k].astype(float).mean()) for k in E109.RUN_MEANS if k in sub},
                                 paa_bmyr_on=float(E83._cs(d['on']).iloc[pb]), paa_bmyr_off=float(E83._cs(d['off']).iloc[pb]),
                                 **E108.unit_scores(d)))
    rd = pd.DataFrame(rows)
    rd.to_csv(out_dir() / 'readings.csv', index=False)
    return rd


def remy_readings() -> pd.DataFrame:
    """remy_readings.csv: the family (E107's as-fitted runs, 'I0') and each search's input at 0.6 and 1 (runs_V.csv) on
    T's population and each wiring of T, APL on, against each of Remy's matrices (her halves and her mean table): ρ over
    the 136 and the 31 of paa and B-myr, the distance of paa from B-myr."""
    tg = remy_targets()
    rv = E83._read_rows(out_dir() / 'runs_V.csv')
    changes = {'I0': v2_fitted(T_WIRINGS)}
    for ch, g in (rv.groupby('change', sort=False) if len(rv) else []):
        if set(g['wiring']) == set(T_WIRINGS):
            changes[ch] = {w: load_run('V', f'V_{ch}_{w}') for w in T_WIRINGS}
    rows = []
    for ch, runs in changes.items():
        for unit, ws in [('T', T_WIRINGS)] + [(w, (w,)) for w in T_WIRINGS]:
            pop = population_on({w: runs[w]['on'][1] for w in ws})
            for t in ('all',) + REMY_HALVES:
                rows.append(dict(change=ch, unit=unit, target=t, rho=remy_rho(pop, tg[t]), rho_pb=remy_rho(pop, tg[t], 'pb'),
                                 paa_bmyr=pb_distance(pop)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'remy_readings.csv', index=False)
    return df


def _v(rd: pd.DataFrame, change: str, unit: str, col: str, readout: str = 'raw') -> float:
    return E109._v(rd, change, unit, col, readout)


def _d(rd: pd.DataFrame, change: str, unit: str, col: str, readout: str = 'raw', ref: str = 'I0') -> float:
    return _v(rd, change, unit, col, readout) - _v(rd, ref, unit, col, readout)


def change_values(rd: pd.DataFrame, runs: pd.DataFrame, change: str, unit: str, wirings: Sequence[str],
                  n_gain_wirings: int) -> Optional[dict]:
    """One change's values against the family at the lab's input ('I0') on `unit` (T or all8): each ρ's change against
    Yang's all recordings, A and B, on the unit and on each wiring (the count gaining), the APL-effect statistics'
    changes against all recordings, the pending items', the runs flagged and the response rate ρ's change (mean over the
    unit's runs). None if a reading is missing."""
    if not ((rd['change'] == change) & (rd['unit'] == unit)).any():
        return None
    out = dict(change=change, unit=unit)
    for c in ('on', 'off'):
        for t in YANG_TARGETS:
            out[f'd_rho_{c}_{t}'] = _d(rd, change, unit, f'rho_{c}_{t}')
            per = {w: _d(rd, change, w, f'rho_{c}_{t}') for w in wirings}
            out[f'd_rho_{c}_{t}_wirings'] = per
            out[f'n_wirings_gain_{c}_{t}'] = int(sum(v >= GAIN[c] for v in per.values()))
        out[f'gain_wirings_ok_{c}'] = bool(all(out[f'n_wirings_gain_{c}_{t}'] >= n_gain_wirings for t in YANG_TARGETS))
    for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson_effect', 'rho_effect', 'loss'):
        out[f'd_{k}_cal'] = _d(rd, change, unit, f'all_{k}', 'cal')
    out['d_loss_raw'] = _d(rd, change, unit, 'all_loss', 'raw')
    out['d_mean_effect_cal'] = _d(rd, change, unit, 'mean_effect', 'cal')
    for k in ('block', 'vcf_on', 'vcf_off', 'struct_wt_on', 'struct_wt_off', 'paa_bmyr_on', 'paa_bmyr_off'):
        out[f'd_{k}'] = _d(rd, change, unit, k)
    out['d_vcf_on_minus_off'] = out['d_vcf_on'] - out['d_vcf_off']
    r = runs[(runs['change'] == change) & runs['wiring'].isin(wirings)]
    r0 = runs[(runs['change'] == 'I0') & runs['wiring'].isin(wirings)]
    out['n_flagged'] = int(r['flagged'].sum())
    out['n_runs'] = int(len(r))
    for c in ('on', 'off'):
        out[f'd_rate_rho_{c}'] = float(r[f'rate_rho_{c}'].astype(float).mean() - r0[f'rate_rho_{c}'].astype(float).mean())
    return out


def _gain(vals: Sequence[dict], c: str) -> Tuple[bool, bool]:
    """(the gain on ρ `c` holds in every value of `vals` against all three of Yang's targets, on the unit and on enough
    wirings; the other ρ is not lower by more than OTHER_SLACK against any)."""
    o = 'off' if c == 'on' else 'on'
    ok = all(v[f'd_rho_{c}_{t}'] >= GAIN[c] for v in vals for t in YANG_TARGETS) and all(v[f'gain_wirings_ok_{c}'] for v in vals)
    other = all(v[f'd_rho_{o}_{t}'] >= -OTHER_SLACK for v in vals for t in YANG_TARGETS)
    return bool(ok), bool(other)


def _measures(vals: Sequence[dict], out: dict) -> None:
    """The other measures (the worst over `vals`) into out['fails'] / out['pending_fails']."""
    fails, pend = {}, {}
    for label, key, sign, slack, kind in CHECKS:
        vs = [v[key] for v in vals]
        w = float(max(vs) if sign > 0 else min(vs))
        out[f'worst_{key}'] = w
        if (sign > 0 and w > slack) or (sign < 0 and w < -slack):
            (fails if kind == 'fails' else pend)[label] = w
    nf = int(sum(v['n_flagged'] for v in vals))
    out['n_flagged'] = nf
    if nf:
        fails['realism (runs outside the bands)'] = nf
    out['fails'], out['pending_fails'] = fails, pend


def _remy_d(remy: pd.DataFrame, change: str, half: str, unit: str = 'T') -> float:
    r = lambda ch: remy[(remy['change'] == ch) & (remy['unit'] == unit) & (remy['target'] == half)]['rho']
    a, b = r(change), r('I0')
    return float(a.iloc[0] - b.iloc[0]) if len(a) and len(b) else np.nan


def gate_r(rd: pd.DataFrame, runs: pd.DataFrame, remy: pd.DataFrame, size: float = GATE_SIZE) -> dict:
    """Part R's gate (the module docstring) on readings (`rd`: change, unit, readout and E108's `unit_scores` columns, T
    and its wirings), runs (change, wiring, flagged, rate_rho_on / _off) and Remy's readings (`remy`: change, unit,
    target, rho). Returns the gate's dict with 'cls' (E109's `class_of`), 'gain_yang_<c>', 'gain_remy'."""
    out = dict(size=size, complete=False)
    vals = {n: change_values(rd, runs, f'{n}_s{size:g}', 'T', T_WIRINGS, R_GAIN_WIRINGS) for n in SEARCHES}
    if any(v is None for v in vals.values()):
        out['note'] = 'readings missing: ' + ' '.join(n for n, v in vals.items() if v is None)
        return out
    for n, v in vals.items():
        v['d_remy_unseen'] = _remy_d(remy, f'{n}_s{size:g}', REMY_UNSEEN[SEARCHES[n]['half']])
        v['d_remy_all'] = _remy_d(remy, f'{n}_s{size:g}', 'all')
    if any(not np.isfinite(vals[n]['d_remy_unseen']) for n in SEARCHES):
        out['note'] = "Remy's readings missing"
        return out
    out['complete'] = True
    out['directions'] = {n: vals[n] for n in FIRST_SPLIT}
    out['second'] = {n: vals[n] for n in SECOND_SPLIT}
    dirs = [vals[n] for n in FIRST_SPLIT]
    remy_first = all(v['d_remy_unseen'] >= GAIN['on'] for v in dirs)
    remy_second = float(np.mean([vals[n]['d_remy_unseen'] for n in SECOND_SPLIT]))
    out['gain_remy'] = bool(remy_first and remy_second > 0)
    out['remy_second_mean'] = remy_second
    for c in ('on', 'off'):
        ok, other = _gain(dirs, c)
        sec = float(np.mean([vals[n][f'd_rho_{c}_all'] for n in SECOND_SPLIT]))
        out[f'second_d_rho_{c}'] = sec
        out[f'gain_yang_{c}'] = bool(ok and other and sec > 0)
        out[f'gain_{c}_parts'] = dict(directions=ok, other=other, second_split=bool(sec > 0), remy=out['gain_remy'])
        out[f'gain_{c}'] = bool(out[f'gain_yang_{c}'] and out['gain_remy'])
        out[f'mean_gain_{c}'] = float(np.mean([vals[n][f'd_rho_{c}_all'] for n in SEARCHES]))
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    ctl = {}
    for c in ('on', 'off'):
        g0 = out[f'mean_gain_{c}']
        perm = [_d(rd, f'{n}_perm{s}', 'T', f'rho_{c}_all') for n in SEARCHES for s in PERM_SEEDS]
        inv = [_d(rd, f'{n}_inv', 'T', f'rho_{c}_all') for n in SEARCHES]
        v = dict(permuted=float(np.mean(perm)) if np.isfinite(perm).all() else np.nan,
                 inverted=float(np.mean(inv)) if np.isfinite(inv).all() else np.nan)
        ctl[c] = dict(candidate=g0, **v, ok=bool(g0 > 0 and all(np.isfinite(x) and x < CONTROL_RATIO * g0 for x in v.values())))
    out['controls'] = ctl
    _measures(dirs, out)
    out['n_flagged_all_four'] = int(sum(v['n_flagged'] for v in vals.values()))
    out['cls'] = E109.class_of(out)
    out['class_without_pending'] = E109.class_of(out, with_pending=False)
    return out


def gate_h(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """Part H's gate for one arm (the module docstring): all eight wirings against the family at the lab's input."""
    out = dict(arm=arm, complete=False)
    v = change_values(rd, runs, arm, 'all8', ALL_WIRINGS, H_GAIN_WIRINGS)
    perms = {p: change_values(rd, runs, p, 'all8', ALL_WIRINGS, H_GAIN_WIRINGS) for p in H_PERMS}
    if v is None or any(x is None for x in perms.values()):
        out['note'] = 'readings missing'
        return out
    out['complete'] = True
    out['directions'] = {arm: v}
    for c in ('on', 'off'):
        ok, other = _gain([v], c)
        out[f'gain_{c}_parts'] = dict(gain=ok, other=other)
        out[f'gain_{c}'] = bool(ok and other)
        out[f'mean_gain_{c}'] = float(v[f'd_rho_{c}_all'])
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    ctl = {}
    for c in ('on', 'off'):
        g0 = out[f'mean_gain_{c}']
        pm = float(np.mean([x[f'd_rho_{c}_all'] for x in perms.values()]))
        ctl[c] = dict(candidate=g0, permuted=pm, ok=bool(g0 > 0 and pm < CONTROL_RATIO * g0))
    out['controls'] = ctl
    _measures([v], out)
    out['cls'] = E109.class_of(out)
    out['class_without_pending'] = E109.class_of(out, with_pending=False)
    return out


def entries(F: Mapping[str, pd.DataFrame], F109: Mapping[str, pd.DataFrame], n_null: int = N_NULL,
            seed: int = 0, recurring: Mapping[str, Sequence[str]] = E109_RECURRING) -> Tuple[pd.DataFrame, dict]:
    """Step 10: per entry (odor x glomerulus of SEARCH_ODORS) the log2 factor of each search of `F`, their mean, the mean
    as read at 0.6, E109's mean over `F109`; and the agreement: Pearson r of the two means over the entries, the number
    of `recurring` entries with the same sign, each with its null (the means of `F` permuted over the glomeruli within
    each odor, `n_null` draws, `default_rng(seed)`; p the share of draws at least as large, with the observed counted)
    and 'agree' (both p at most AGREE_P)."""
    odors = list(SEARCH_ODORS)
    gl = list(next(iter(F.values())).index)
    L = np.mean([np.log2(F[n][odors].to_numpy(dtype=float)) for n in F], axis=0)          # glomeruli x odors
    L9 = np.mean([np.log2(F109[n].reindex(index=gl)[odors].to_numpy(dtype=float)) for n in F109], axis=0)
    rec = [(gl.index(g), odors.index(o)) for o, gs in recurring.items() for g in gs if g in gl]

    def stats_of(M):
        r = float(np.corrcoef(M.ravel(), L9.ravel())[0, 1])
        s = int(sum(np.sign(M[i, j]) == np.sign(L9[i, j]) and M[i, j] != 0 for i, j in rec))
        return r, s
    r0, s0 = stats_of(L)
    rng = np.random.default_rng(seed)
    rn, sn = np.empty(n_null), np.empty(n_null)
    for k in range(n_null):
        M = np.column_stack([L[rng.permutation(len(gl)), j] for j in range(len(odors))])
        rn[k], sn[k] = stats_of(M)
    p_r = float((1 + (rn >= r0).sum()) / (1 + n_null))
    p_s = float((1 + (sn >= s0).sum()) / (1 + n_null))
    rows = []
    for j, o in enumerate(odors):
        for i, g in enumerate(gl):
            row = dict(odor=o, glomerulus=g, **{f'log2_{n}': float(np.log2(F[n].loc[g, o])) for n in F},
                       mean_log2=float(L[i, j]), mean_log2_read=float(L[i, j] * GATE_SIZE), e109_mean_log2=float(L9[i, j]),
                       e109_recurring=bool(g in recurring.get(o, ())), same_sign=bool(np.sign(L[i, j]) == np.sign(L9[i, j])))
            rows.append(row)
    agree = dict(pearson_r=r0, pearson_null_q95=float(np.quantile(rn, 0.95)), pearson_p=p_r, n_recurring=len(rec),
                 sign_agree=s0, sign_null_q95=float(np.quantile(sn, 0.95)), sign_p=p_s, n_null=n_null, seed=seed,
                 agree=bool(p_r <= AGREE_P and p_s <= AGREE_P))
    return pd.DataFrame(rows), agree


def outcomes(gr: Mapping, gh: Mapping, agree: Optional[bool]) -> List[str]:
    """The outcome letters (the module docstring) from part R's gate `gr`, H_prof's `gh` and step 10's agreement."""
    out = []
    if gh.get('complete') and gh.get('cls') == 'passes':
        out.append('e')
    if not gr.get('complete'):
        return out
    if gr['cls'] == 'passes':
        out.append('a' if agree else 'b')
    if gr.get('gain_on') and not gr.get('gain_off'):
        out.append('c')
    if gr.get('gain_remy') and not (gr.get('gain_yang_on') or gr.get('gain_yang_off')):
        out.append('d')
    ctl = gr.get('controls') or {}
    as_well = all(max(ctl.get(c, {}).get('permuted', -np.inf), ctl.get(c, {}).get('inverted', -np.inf))
                  >= ctl.get(c, {}).get('candidate', np.inf) for c in ('on', 'off'))
    if gr['cls'] == 'does not carry' and gh.get('cls', 'does not carry') == 'does not carry' and as_well:
        out.append('f')
    return out


def searches_table() -> pd.DataFrame:
    """searches.csv: per search, its start and last parent on S against the searched half, the unseen half and her mean
    table; per odor the rms of its log2 factors and the factors at a bound."""
    rows = []
    for n, cfg in SEARCHES.items():
        p = out_dir() / f'search_history_{n}.csv'
        if not p.exists():
            continue
        h = pd.read_csv(p)
        s0, last = h[h.odor == 'start'], h[h.odor != 'start'].tail(1)
        row = dict(search=n, half=cfg['half'], unseen=REMY_UNSEEN[cfg['half']], finished=search_finished(n),
                   generations_done=int((h.odor != 'start').sum()))
        for k in ('rho_half', 'rho_unseen', 'rho_all'):
            row[f'start_{k}'] = float(s0[k].iloc[0]) if len(s0) else np.nan
            row[f'last_{k}'] = float(last[k].iloc[0]) if len(last) else np.nan
        row['last_parent_score'] = float(last['parent_score'].iloc[0]) if len(last) else np.nan
        if search_finished(n):
            f = pd.read_csv(factors_path(n))
            for o, g in f.groupby('odor', sort=False):
                row[f'{o}_log2_rms'] = float(np.sqrt(np.mean(g.log2_factor ** 2)))
                row[f'{o}_n_at_bound'] = int((np.abs(np.abs(g.log2_factor) - 1.0) < 1e-9).sum())
        rows.append(row)
    st = pd.DataFrame(rows)
    st.to_csv(out_dir() / 'searches.csv', index=False)
    return st


def e112_ended() -> Tuple[bool, str]:
    """Whether E112 has ended (its index row in experiment_priorities.md says 'completed'), and the row's state."""
    p = H.shared_dir() / 'experiment_priorities.md'
    for line in p.read_text().splitlines():
        if line.startswith('| E112 |'):
            st = line.split('|')[3].strip()
            return st == 'completed', st
    return False, 'not found'


def e115_class_table(runs: pd.DataFrame) -> pd.DataFrame:
    """Step 11, beside: the 66 megamat pairs by E115's classes (results/E115_two_yardsticks/a_classes.csv, her cosine
    reading where it has one): per change (the family and the four searches at 0.6, T's population, raw counts) and
    condition, the summed contribution to 1 - ρ on odor-odor distances over the 66 against all of Yang's recordings
    (E106's `rank_residuals`) per class."""
    from experiments import E106_rho_gap_parts as E106
    from experiments import E108_paa_extrapolation as E108
    p = H.results_dir('E115_two_yardsticks') / 'a_classes.csv'
    if not p.exists():
        return pd.DataFrame()
    cl = pd.read_csv(p)
    metric = 'cos' if 'cos' in set(cl['remy']) else sorted(set(cl['remy']))[0]
    cl = cl[cl['remy'] == metric]
    dt = E108.dist_targets()['all']
    rows = []
    for ch in ['I0'] + [f'{n}_s{GATE_SIZE:g}' for n in SEARCHES]:
        if set(runs[runs['change'] == ch]['wiring']) < set(T_WIRINGS):
            continue
        d = E83.population([E83.pair_distances(run_outputs(ch, w), 'raw') for w in T_WIRINGS])
        for c in ('on', 'off'):
            pairs = [tuple(x) for x in cl[cl['condition'] == c][['odor_a', 'odor_b']].to_numpy()]
            m = E83._cs(d[c])
            idx = E106.select_pairs(m.index, pairs)
            rr = E106.rank_residuals(m.loc[idx], dt[c])
            cls = {tuple(sorted(k)): v for k, v in zip(pairs, cl[cl['condition'] == c]['cls'])}
            rr['cls'] = [cls.get(tuple(sorted(k)), 'neither') for k in rr.index]
            for k, g in rr.groupby('cls'):
                rows.append(dict(change=ch, condition=c, cls=k, n_pairs=len(g), contribution=float(g['contribution'].sum()),
                                 remy_metric=metric))
    df = pd.DataFrame(rows)
    if len(df):
        ref = df[df.change == 'I0'].set_index(['condition', 'cls'])['contribution']
        df['d_contribution'] = [r.contribution - ref.get((r.condition, r.cls), np.nan) for r in df.itertuples()]
    return df


def combined_table(rd: pd.DataFrame) -> pd.DataFrame:
    """Step 11, beside: the two ρ's against all of Yang's recordings on T (raw counts) for the family, the best search of
    part R at 0.6, E109's B2 at 0.6 (its archived held runs) and the two combined (`comb_<best>_B2`), changes from the
    family."""
    from experiments import E108_paa_extrapolation as E108
    b = best_search()
    if b is None:
        return pd.DataFrame()
    rows = []
    base = {c: _v(rd, 'I0', 'T', f'rho_{c}_all') for c in ('on', 'off')}
    b2 = {w: E109.load_run('R', f'{COMBINED_WITH}_s{GATE_SIZE:g}_{w}') for w in T_WIRINGS}
    sc_b2 = (E108.unit_scores(E83.population([E83.pair_distances(o, 'raw') for o in b2.values()]))
             if all(o is not None for o in b2.values()) else None)
    for label, ch in (('the family at the lab input', 'I0'), (f'{b} at 0.6', f'{b}_s{GATE_SIZE:g}'),
                      (f'{COMBINED_WITH} (E109) at 0.6', None), (f'{b} x {COMBINED_WITH} at 0.6', f'comb_{b}_{COMBINED_WITH}')):
        if ch is None:
            if sc_b2 is None:
                continue
            v = {c: sc_b2[f'rho_{c}_all'] for c in ('on', 'off')}
        else:
            v = {c: _v(rd, ch, 'T', f'rho_{c}_all') for c in ('on', 'off')}
        rows.append(dict(configuration=label, rho_on_all=v['on'], rho_off_all=v['off'], d_rho_on=v['on'] - base['on'],
                         d_rho_off=v['off'] - base['off']))
    return pd.DataFrame(rows)


def _f(x, spec: str = '.3f') -> str:
    return E109._f(x, spec)


def _md(df: pd.DataFrame, spec: str = '.4f') -> str:
    return E109._md_table(df, spec)


def comparison_table(rd: pd.DataFrame, runs: pd.DataFrame, remy: pd.DataFrame, gr: dict, gh: Dict[str, dict],
                     ent: pd.DataFrame, agree: dict, letters: List[str], beside: Dict[str, object]) -> str:
    """comparison_table.md: the tables the note reads first."""
    L = ['# E116: the input of paa and B-myr searched against a second set of flies and taken from Hallem\'s receptor data',
         '', 'Changes are from the family at the lab\'s input on the same wirings (raw counts unless said). "all", A, B: '
         "Yang's recordings (all, and E77's halves). Remy's halves a, b (split 1) and c, d (split 2).", '']
    L += ['## Part H: Hallem\'s paa, no search (all eight wirings)', '',
          '| arm | Δρ on, all / A / B | Δρ off, all / A / B | wirings gaining on / off (all) | pattern part (cal) | '
          'raw loss | mean APL effect | response rate ρ on / off | runs flagged | paa from B-myr on / off | class |',
          '|---|---|---|---|---:|---:|---:|---|---:|---|---|']
    for arm in H_ARMS + H_PERMS:
        v = change_values(rd, runs, arm, 'all8', ALL_WIRINGS, H_GAIN_WIRINGS)
        if v is None:
            continue
        g = gh.get(arm, {})
        L.append(f"| {arm} | {' / '.join(_f(v[f'd_rho_on_{t}'], '+.4f') for t in YANG_TARGETS)} | "
                 f"{' / '.join(_f(v[f'd_rho_off_{t}'], '+.4f') for t in YANG_TARGETS)} | "
                 f"{v['n_wirings_gain_on_all']} / {v['n_wirings_gain_off_all']} | {_f(v['d_pattern_cal'], '+.4f')} | "
                 f"{_f(v['d_loss_raw'], '+.4f')} | {_f(v['d_mean_effect_cal'], '+.4f')} | "
                 f"{_f(v['d_rate_rho_on'], '+.3f')} / {_f(v['d_rate_rho_off'], '+.3f')} | {v['n_flagged']} | "
                 f"{_f(v['d_paa_bmyr_on'], '+.3f')} / {_f(v['d_paa_bmyr_off'], '+.3f')} | {g.get('cls', '')} |")
    L += ['', f"Family on all eight against all recordings: ρ {_f(_v(rd, 'I0', 'all8', 'rho_on_all'))} / "
          f"{_f(_v(rd, 'I0', 'all8', 'rho_off_all'))}; paa from B-myr {_f(_v(rd, 'I0', 'all8', 'paa_bmyr_on'))} / "
          f"{_f(_v(rd, 'I0', 'all8', 'paa_bmyr_off'))}. {H_CLASS_ARM} carries part H's class.", '']
    for arm, g in gh.items():
        if g.get('complete'):
            L.append(f"- {arm}: class '{g['cls']}' (without the pending items '{g['class_without_pending']}'); controls "
                     + '; '.join(f"{c}: candidate {_f(g['controls'][c]['candidate'], '+.4f')}, permuted "
                                 f"{_f(g['controls'][c]['permuted'], '+.4f')}, ok {g['controls'][c]['ok']}" for c in ('on', 'off'))
                     + f"; fails {g['fails'] or 'none'}; pending {g['pending_fails'] or 'none'}")
    st = searches_table()
    if len(st):
        L += ['', '## Part R: the searches on S against Remy\'s halves (as fitted, APL on, 136 pairs)', '',
              _md(st[[c for c in st.columns if c in ('search', 'half', 'unseen', 'generations_done', 'start_rho_half',
                                                         'last_rho_half', 'start_rho_unseen', 'last_rho_unseen',
                                                         'start_rho_all', 'last_rho_all', 'paa_n_at_bound',
                                                         'B-myr_n_at_bound')]], '.4f')]
    rows = []
    for n in SEARCHES:
        for size in (GATE_SIZE, FULL_SIZE):
            ch = f'{n}_s{size:g}'
            v = change_values(rd, runs, ch, 'T', T_WIRINGS, R_GAIN_WIRINGS)
            if v is None:
                continue
            rows.append(dict(change=ch, **{f'd_rho_{c}_{t}': v[f'd_rho_{c}_{t}'] for c in ('on', 'off') for t in YANG_TARGETS},
                             wirings_on=v['n_wirings_gain_on_all'], wirings_off=v['n_wirings_gain_off_all'],
                             d_remy_unseen=_remy_d(remy, ch, REMY_UNSEEN[SEARCHES[n]['half']]),
                             d_remy_all=_remy_d(remy, ch, 'all'), d_paa_bmyr_on=v['d_paa_bmyr_on'],
                             d_paa_bmyr_off=v['d_paa_bmyr_off'], d_rate_rho_on=v['d_rate_rho_on'],
                             d_rate_rho_off=v['d_rate_rho_off'], n_flagged=v['n_flagged']))
    if rows:
        L += ['', "## The held readings on T (Yang's panel, both response counts re-held; Remy's panel as fitted)", '',
              f"Family on T: ρ against all recordings {_f(_v(rd, 'I0', 'T', 'rho_on_all'))} / "
              f"{_f(_v(rd, 'I0', 'T', 'rho_off_all'))}; against A {_f(_v(rd, 'I0', 'T', 'rho_on_A'))} / "
              f"{_f(_v(rd, 'I0', 'T', 'rho_off_A'))}; against B {_f(_v(rd, 'I0', 'T', 'rho_on_B'))} / "
              f"{_f(_v(rd, 'I0', 'T', 'rho_off_B'))}; paa from B-myr {_f(_v(rd, 'I0', 'T', 'paa_bmyr_on'))} / "
              f"{_f(_v(rd, 'I0', 'T', 'paa_bmyr_off'))}.", '', _md(pd.DataFrame(rows), '+.4f')]
    if gr.get('complete'):
        L += ['', '## Part R\'s gate (0.6)', '', '| statistic | ' + ' | '.join(FIRST_SPLIT) + ' | second split (mean of '
              + ', '.join(SECOND_SPLIT) + ') | threshold |', '|---|' + '---:|' * (len(FIRST_SPLIT) + 1) + '---|']
        D, S2 = gr['directions'], gr['second']
        sm = lambda k: float(np.mean([S2[n][k] for n in SECOND_SPLIT]))
        for lab, k, thr in (('Δρ APL on, all recordings', 'd_rho_on_all', '≥ 0.010'), ('Δρ APL on, A', 'd_rho_on_A', '≥ 0.010'),
                            ('Δρ APL on, B', 'd_rho_on_B', '≥ 0.010'), ('Δρ APL off, all recordings', 'd_rho_off_all', '≥ 0.015'),
                            ('Δρ APL off, A', 'd_rho_off_A', '≥ 0.015'), ('Δρ APL off, B', 'd_rho_off_B', '≥ 0.015'),
                            ("Δρ against Remy's unseen half (APL on, 136 pairs)", 'd_remy_unseen', '≥ 0.010 (second split > 0)'),
                            ('pattern part (cal)', 'd_pattern_cal', '≤ +0.003'), ('amplitude part', 'd_amplitude_cal', 'beside'),
                            ('shape part', 'd_shape_cal', 'beside'), ('spread ratio', 'd_sd_ratio_cal', 'beside'),
                            ('Pearson r on the per-pair APL effect', 'd_pearson_effect_cal', 'beside'),
                            ('ρ on the per-pair APL effect', 'd_rho_effect_cal', 'beside'),
                            ('raw-count loss', 'd_loss_raw', '≤ +0.005'), ('mean APL effect (cal)', 'd_mean_effect_cal', '≥ -0.01'),
                            ('response rate ρ, APL on', 'd_rate_rho_on', '≥ -0.05'), ('response rate ρ, APL off', 'd_rate_rho_off', '≥ -0.05'),
                            ('the block (pending)', 'd_block', '≥ -0.010'), ('natural-source alignment on (pending)', 'd_vcf_on', '≥ -0.05'),
                            ('natural-source alignment off (pending)', 'd_vcf_off', '≥ -0.05'),
                            ('structure vs wild-type somata (pending)', 'd_struct_wt_on', '≥ -0.05')):
            L.append(f"| {lab} | " + ' | '.join(_f(D[n][k], '+.4f') for n in FIRST_SPLIT) + f" | {_f(sm(k), '+.4f')} | {thr} |")
        L.append('| wirings gaining, on / off (all recordings) | ' + ' | '.join(
            f"{D[n]['n_wirings_gain_on_all']} / {D[n]['n_wirings_gain_off_all']}" for n in FIRST_SPLIT) + ' | | ≥ 3 of 4 |')
        L.append('| runs outside the bands | ' + ' | '.join(str(D[n]['n_flagged']) for n in FIRST_SPLIT) + ' | | 0 |')
        L += ['', '| ρ on odor-odor distances | candidate (mean of four, all recordings) | permuted (16) | inverted (4) | '
              'controls hold |', '|---|---:|---:|---:|---|']
        for c in ('on', 'off'):
            k = gr['controls'][c]
            L.append(f"| APL {c} | {_f(k['candidate'], '+.4f')} | {_f(k['permuted'], '+.4f')} | {_f(k['inverted'], '+.4f')} | "
                     f"{'yes' if k['ok'] else 'no'} |")
        L += ['', f"Class: '{gr['cls']}' (without the pending items '{gr['class_without_pending']}'); gain on / off: "
              f"{gr['gain_on']} / {gr['gain_off']} (Yang: {gr['gain_yang_on']} / {gr['gain_yang_off']}; Remy's unseen half: "
              f"{gr['gain_remy']}); fails {gr['fails'] or 'none'}; pending {gr['pending_fails'] or 'none'}."]
    elif gr:
        L += ['', f"Part R's gate: not complete ({gr.get('note', '')})."]
    if len(ent):
        L += ['', "## Step 10: the entries against E109's", '',
              f"Pearson r over the 72 entries {_f(agree['pearson_r'])} (null 95th percentile {_f(agree['pearson_null_q95'])}, "
              f"p {_f(agree['pearson_p'], '.4f')}); E109's eight recurring entries with the same sign {agree['sign_agree']} of "
              f"{agree['n_recurring']} (null 95th percentile {_f(agree['sign_null_q95'], '.1f')}, p {_f(agree['sign_p'], '.4f')}); "
              f"agree: {agree['agree']}.", '']
        sub = ent[ent['e109_recurring'] | (ent['mean_log2'].abs() >= 0.5)].copy()
        cols = ['odor', 'glomerulus', 'lab_delta', 'mean_log2', 'mean_log2_read', 'e109_mean_log2', 'e109_recurring',
                'flies_sd_log2', 'read_in_flies_sd']
        L.append(_md(sub[[c for c in cols if c in sub.columns]], '.3f'))
    for k, df in beside.items():
        if isinstance(df, pd.DataFrame) and len(df):
            L += ['', f'## Beside: {k}', '', _md(df, '.4f')]
        elif isinstance(df, str):
            L += ['', f'## Beside: {k}', '', df]
    L += ['', f"## Outcomes: {', '.join('(' + x + ')' for x in letters) if letters else 'none of (a) to (f)'}", '']
    return '\n'.join(L) + '\n'


def realism_table(runs: pd.DataFrame) -> pd.DataFrame:
    """realism.csv: per change (the family on all eight and on T, every reading), the mean and range over its runs of
    the response rates, their ratio, response rate ρ, the odors-per-KC distances and the largest counts, and the runs
    flagged."""
    rows = []
    for ch, g in runs.groupby('change', sort=False):
        units = [('all8', g), ('T', g[g['wiring'].isin(T_WIRINGS)])] if ch == 'I0' else [('runs', g)]
        for unit, s in units:
            row = dict(change=ch, unit=unit, n_runs=len(s), n_flagged=int(s['flagged'].sum()))
            for k in ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv',
                      'max_count_on', 'max_count_off', 'thr_factor', 'apl_factor'):
                if k in s:
                    x = s[k].astype(float)
                    row.update({k: float(x.mean()), f'{k}_min': float(x.min()), f'{k}_max': float(x.max())})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'realism.csv', index=False)
    return df


def score() -> dict:
    """Part score (no simulations; the module docstring)."""
    out = out_dir()
    runs = runs_table()
    rd = readings(runs)
    remy = remy_readings()
    gh = {arm: gate_h(rd, runs, arm) for arm in H_ARMS}
    gr = gate_r(rd, runs, remy)
    ent, agree = pd.DataFrame(), {}
    if all(search_finished(n) for n in SEARCHES):
        F = {n: search_factors(n) for n in SEARCHES}
        F109 = {n: E109.search_factors(n) for n in E109_SIX}
        ent, agree = entries(F, F109)
        fs = E109.flies_spread()
        ent = ent.merge(fs[['odor', 'glomerulus', 'flies_sd_log2']], on=['odor', 'glomerulus'], how='left')
        ent['read_in_flies_sd'] = ent['mean_log2_read'] / ent['flies_sd_log2']
        ent['mean_in_flies_sd'] = ent['mean_log2'] / ent['flies_sd_log2']
        ent.to_csv(out / 'entries.csv', index=False)
        _write_json(out / 'entries_agreement.json', agree)
    ended, state = e112_ended()
    beside = {"E115's classes of the megamat pairs (summed contributions to 1 - ρ, T, all recordings)": e115_class_table(runs),
              "the best search of part R combined with E109's B2 (T, all recordings)": combined_table(rd),
              "E112's reading": ('E112 had ended: its reading is to be applied by its own tables.' if ended else
                                 f"E112 had not ended at scoring (its index row: '{state}'), so the two ρ's are not given "
                                 'under its reading.')}
    for k, df in beside.items():
        if isinstance(df, pd.DataFrame) and len(df):
            df.to_csv(out / ('beside_e115_classes.csv' if 'E115' in k else 'beside_combined.csv'), index=False)
    letters = outcomes(gr, gh.get(H_CLASS_ARM, {}), agree.get('agree') if agree else None)
    rows = [E109._flat(g, which='H', arm=a, carrying=bool(a == H_CLASS_ARM)) for a, g in gh.items()]
    rows.append(E109._flat(gr, which='R', carrying=True))
    pd.DataFrame(rows).to_csv(out / 'gate.csv', index=False)
    verdict = dict(part_h={a: {k: g.get(k) for k in ('complete', 'cls', 'class_without_pending', 'gain_on', 'gain_off',
                                                       'mean_gain_on', 'mean_gain_off', 'controls', 'fails', 'pending_fails')}
                           for a, g in gh.items()},
                   part_r={k: gr.get(k) for k in ('complete', 'cls', 'class_without_pending', 'gain_on', 'gain_off',
                                                   'gain_yang_on', 'gain_yang_off', 'gain_remy', 'mean_gain_on',
                                                   'mean_gain_off', 'controls', 'fails', 'pending_fails', 'note')},
                   entries=agree, outcomes=letters, e112=dict(ended=ended, state=state))
    _write_json(out / 'gate.json', dict(verdict=verdict, gates=dict(H=gh, R=gr), **H.code_commits()))
    realism_table(runs)
    (out / 'comparison_table.md').write_text(comparison_table(rd, runs, remy, gr, gh, ent, agree, letters, beside))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), script=SCRIPT, **H.code_commits(),
                                                    approx_conc_note=H.APPROX_CONC_NOTE))
    print(json.dumps(verdict, indent=1, default=str), flush=True)
    return verdict


def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive (parts H, R, S, V) and model_pair_distances_long.csv (raw spike counts,
    Yang's panel only: parts H and R) for kcq/apl_effect.py."""
    out = out_dir()
    rows, parts = [], []
    for part in ('H', 'R'):
        for r in E83._read_rows(out / f'runs_{part}.csv').to_dict('records'):
            rows.append(dict(tag=r['tag'], variant=r['tag'], part=part, change=r['change'], kind=r.get('kind'),
                             search=r.get('search'), size=r.get('size'), seed=r.get('seed'), wiring=r['wiring'], base='bouton',
                             k=E83.BREADTH_K, breadth_k=E83.BREADTH_K, ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'],
                             orn_input=r.get('orn_input'), **{k: r[k] for k in E83.PARAM_KEYS}, kc_kc_inh=r['wKCKC_scale'],
                             thr_factor=r['thr_factor'], apl_factor=r['apl_factor'], rate_on=r['rate_on'],
                             rate_off=r['rate_off'], realism_flag=r.get('realism_flag'), subdir=part, panel='APL16',
                             kind_text=(f"C2's family's {r['wiring']} fit with input {r['change']}; breadth recomputed for the "
                                        "input; the somatic threshold and the tied APL scale re-held at the fit's APL-off and "
                                        "APL-on counts")))
    for part, f in (('S', None), ('V', 'runs_V.csv')):
        tabs = ([pd.read_csv(out / f'search_runs_{n}.csv') for n in SEARCHES if (out / f'search_runs_{n}.csv').exists()]
                if f is None else [E83._read_rows(out / f)])
        for t in tabs:
            for r in t.to_dict('records'):
                rows.append(dict(tag=r['tag'], variant=r['tag'], part=part, change=r.get('change', f"{r['search']}_search"),
                                 kind='search_final' if part == 'S' else 'remy_panel', search=r['search'],
                                 size=r.get('size', 1.0), wiring=r['wiring'], base='bouton', k=E83.BREADTH_K,
                                 breadth_k=E83.BREADTH_K, ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'],
                                 **{k: r[k] for k in E83.PARAM_KEYS}, kc_kc_inh=r['wKCKC_scale'], thr_factor=1.0, apl_factor=1.0,
                                 rate_on=r['rate_on'], rate_off=r['rate_off'], subdir=part, panel='validation2',
                                 kind_text=(f"C2's family's {r['wiring']} fit on the validation2 panel as fitted (nothing "
                                            f"re-held), search {r['search']}'s paa and B-myr factors")))
    for r in rows:
        if r['panel'] != 'APL16':
            continue
        o = load_run(r['subdir'], r['tag'])
        if o is None:
            continue
        parts.append(pd.concat({(r['tag'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['tag'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'H', 'search', 'read', 'remy', 'score', 'variants'))
    ap.add_argument('--search', nargs='+', default=None, choices=list(SEARCHES))
    ap.add_argument('--threads', type=int, default=E109.WORKER_THREADS, help="part search / 0: each worker's OpenMP threads")
    ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--max-points', type=int, default=None)
    ap.add_argument('--generations', type=int, default=GENS, help='part search: generations per odor (smoke tests)')
    ap.add_argument('--pop', type=int, default=POP, help='part search: offspring per generation (smoke tests)')
    a = ap.parse_args()
    if a.part == '0':
        part0(a.threads)
    elif a.part == 'H':
        part_h(a.max_points)
    elif a.part == 'search':
        for n in (a.search or list(SEARCHES)):
            search(n, resume=not a.no_resume, threads=a.threads, gens=a.generations, pop=a.pop)
    elif a.part == 'read':
        part_read(a.search, a.max_points)
    elif a.part == 'remy':
        part_remy(a.search)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()


if __name__ == '__main__':
    main()
