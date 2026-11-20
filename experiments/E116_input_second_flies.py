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
  four searches' mean factors over the glomeruli within each odor (10000 draws, seed 0); "they agree" when both are
  above the null's 95th percentile. Each entry also in the validation2 ORN flies' SDs (E109's `flies_spread`).
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
REMY_SPLITS = {0: ('a', 'b'), 1: ('c', 'd')}
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
YANG_TARGETS = ('all', 'A', 'B')
H_ARMS = ('H_prof', 'H_meas', 'H_dm6')
H_PERMS = tuple(f'H_perm{s}' for s in (1, 2, 3, 4))
H_GAIN_WIRINGS = 6
N_NULL = 10000
E109_RECURRING = {'paa': ('DM1', 'DL5', 'VM7d', 'DP1l', 'DP1m'), 'B-myr': ('DM1', 'VC3', 'VM5d')}
COMBINED_WITH = 'B2'                       # E109's search whose factors the best of part R is combined with


def out_dir() -> Path:
    """results/E116_input_second_flies/ (E116_OUT_DIR in the environment instead, for smoke tests)."""
    raise NotImplementedError


def local_dir(sub: str = '') -> Path:
    """outputs/E116_input_second_flies/<sub> (E116_LOCAL_DIR in the environment instead)."""
    raise NotImplementedError


def split_flies(flies: Sequence[str], held: pd.DataFrame, seed: int, pairs: Sequence[Tuple[str, str]],
                min_per_half: int = 1, tries: int = 1000) -> Tuple[List[str], List[str]]:
    """E77's rule on Remy's flies: the first of up to `tries` permutations (`default_rng(seed)`) whose first floor(n/2)
    flies and the rest each hold every pair of `pairs` with `min_per_half` or more flies (`held`: flies x pairs, True
    where a fly holds a pair)."""
    raise NotImplementedError


def remy_halves() -> Dict[str, List[str]]:
    """{'a', 'b', 'c', 'd': her flies} (`split_flies` on her per-fly cosine table, the 136 pairs)."""
    raise NotImplementedError


def remy_targets() -> Dict[str, pd.Series]:
    """{'a', 'b', 'c', 'd': the half's mean over flies, 'all': her mean table (cosine)} as canonical pair Series."""
    raise NotImplementedError


def population_on(counts: Mapping[str, pd.DataFrame]) -> pd.Series:
    """The mean over wirings of the cosine distance on spike counts (APL on), E106's `mean_pairs`."""
    raise NotImplementedError


class RemyObjective:
    """ρ on odor-odor distances of a population (APL on) over her 136 pairs against one half's mean."""

    def __init__(self, half: str):
        raise NotImplementedError

    def __call__(self, on: Mapping[str, pd.DataFrame]) -> float:
        raise NotImplementedError


def score_of(rho_on: float, fracs, base_fracs) -> dict:
    """A point's score: ρ minus E83's rate penalty on the searched odor's fractions."""
    raise NotImplementedError


def search(name: str, resume: bool = True, threads: int = E109.WORKER_THREADS, gens: int = GENS,
           pop: int = POP) -> pd.DataFrame:
    """Search `name` (SEARCHES; the module docstring)."""
    raise NotImplementedError


def hallem_inputs(lab: Optional[pd.DataFrame] = None) -> Dict[str, Tuple[pd.DataFrame, dict]]:
    """Part H's inputs: {arm: (input, description)} for H_ARMS and H_PERMS."""
    raise NotImplementedError


def part0(threads: int = E109.WORKER_THREADS) -> dict:
    raise NotImplementedError


def part_h(max_points: Optional[int] = None) -> pd.DataFrame:
    raise NotImplementedError


def part_read(names: Optional[Sequence[str]] = None, max_points: Optional[int] = None) -> pd.DataFrame:
    raise NotImplementedError


def part_remy(names: Optional[Sequence[str]] = None) -> pd.DataFrame:
    raise NotImplementedError


def gate_r(rd: pd.DataFrame, runs: pd.DataFrame, remy: pd.DataFrame) -> dict:
    """Part R's gate (the module docstring) on readings (`rd`: change, unit, readout and E108's `unit_scores` columns),
    runs (change, wiring, flagged, rate_rho_on / _off) and Remy's readings (`remy`: change, unit, half, rho)."""
    raise NotImplementedError


def gate_h(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """Part H's gate for one arm."""
    raise NotImplementedError


def entries(F: Mapping[str, pd.DataFrame], F109: Mapping[str, pd.DataFrame], n_null: int = N_NULL,
            seed: int = 0) -> Tuple[pd.DataFrame, dict]:
    """Step 10: the entries' table and the agreement with E109's (correlation, sign agreement, their nulls)."""
    raise NotImplementedError


def outcomes(gr: Mapping, gh: Mapping, agree: Optional[bool]) -> List[str]:
    """The outcome letters (the module docstring)."""
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
