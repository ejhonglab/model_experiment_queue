#!/usr/bin/env python3
"""E93: the measured input for the glomeruli that the model's input leaves at zero.

The entry is item E93 of experiment_priorities.md (proposed by viki 2026-11-13 20:37; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-13 20:12). The model's ORN input (`H.apl16_orn_deltas()`) has the 36
glomeruli both of the lab's panels have; the hemibrain wiring has 54 and `fit_mb_model` zero-fills the other 18. For ten
of them responses to the panel's odors have been measured: VA4 and VL1 in the lab's megamat ORN recordings
(`megamat_orn_deltas()`, the 12 megamat odors), DL4 and VM4 in its validation recordings (`validation2_orn_deltas()`,
paa, -aPine, euc, B-myr), and DA4m, DA3, DA4l, VA1v, DL3, VA1d (15 odors, all but euc) and DL4 (the 12 megamat odors) in
Hallem's receptor data (drosolf, as E38's `hallem_deltas`; 10^-2), each odor's values multiplied by the slope through the
origin of the lab's input on Hallem's over the 16 glomeruli both have for that odor. Unmeasured entries stay 0 (26 of the
160).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `saved_fit`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's (`recording_halves()`,
seed 0), C / D the second split's (seed 1; E89's `targets89`). Every input is fixed before any run and none is fitted to
the KC recordings.

The inputs (INPUTS; glomeruli x the 16 odors, the lab's 36 rows unchanged in all of them; inputs/<name>.csv with sha256):
- I0, the lab's input (the reference: the saved runs); I0z, I0 with the ten rows written out as zeros (part 0).
- F, the fuller input (the primary): I0 plus the ten rows above.
- Bounds, on S only: L (the lab's four rows, no Hallem), U (F with Hallem's values as recorded at 10^-2), V (VA4's row
  alone).
- Controls, on T: PL<seed> (F with the ten rows' labels permuted among the ten glomeruli, no row left in place) and
  PO<seed> (F with each added row's values permuted over the odors, megamat and validation odors separately). Seed 1
  each; seed 2 each too if F moves T's population (pattern part under the readout against all recordings) by
  CONTROL_DRAW2 or more either way.
- DD, only if part A2 says "in the list's direction": F with the PN dendrites' departures put on the input at the
  entries whose halves agree in sign, bounded to a factor of DD_BOUND (on all eight wirings).

Parts (CLI `--part`):
- inputs: the inputs, added_rows.csv, hallem_scaling.csv, inputs/ (no simulations).
- 0, the identity: on each wiring, I0z's simulator (the breadth recomputed for it, `E85.build`, which also checks that
  `J.bouton_base` returns the same wiring and APL weights for a 46-row input) reproduces the saved APL-on and APL-off runs
  with 0 count differences, and its breadth equals the lab input's exactly (identity.json). If not, the experiment stops
  there (outcome (h)).
- A (no simulations): A1 what is added (a1_reliability.csv: halves of the flies for VA4, VL1 (olfdb dataset 1), DL4,
  VM4 (dataset 4), VA4 in the ORNs against the PN dendrites; hallem_scaling.csv; input_shares.csv: each odor's added
  share of claw-weighted positive input; input_distances.csv: the input distances by group with and without the rows;
  who_responds.csv: in the saved runs of the eight wirings the response rate and the odors per KC by a KC's share of
  zero-filled claws, of claws on the ten, on VA4 and on DA4l or DA3). A2 the PN dendrites at the searched input's
  entries (E86's partE_entries.csv on the glomeruli the GH146, pebbled and model tables share): the departure, log2 of
  the dendrites' relative response over the ORNs' (E90's `al_table`), normalized within a glomerulus, within an odor and
  with the model's own PN stage (E90's saved PN rates) in place of the ORNs; intervals from DEP_BOOT resamples of both
  tables' flies; the sign agreement and ρ with the six searches' mean log2 factor against a null of the departures
  permuted within each odor (NULL_DRAWS); the halves' agreement; the rule (A2_RULE, fixed in the entry).
  dendrite_departures.csv, dendrite_reading.json.
- B: the inputs at fixed scales with both response rates re-held (E85's `held_reading`: the breadth recomputed for the
  input, the somatic threshold bisected to the fit's APL-off count, the tied APL scale to its APL-on count); F on all
  eight, L, U, V on S, the controls on T (and DD if A2 says so). Resumable; runs_B.csv (absolute parameters, the factors,
  realism under both rules, the largest spike counts, the named odors' response fractions, the share of KCs responding to
  no odor, `H.code_commits()` of the process, which must be clean); runs under outputs/E93_missing_glomeruli_input/B/.
- score (C): readings.csv (per input and unit, every target of E89's `targets89`, the 120 and the 66 megamat pairs, raw
  counts, the calibrated readout and E73's per-KC reading), spread.csv (the pattern part's amplitude and shape parts
  beside every pattern part, `H.apl_effect_spread_shape` with the target through `H.canonical_pairs`), apl_free.csv (E87's
  target without APL), coresponse.csv (E89's co-response change), line.csv (the groups' APL-on distances against E92's
  part A line), who_responds_runs.csv, gate.csv (C's gate on T), equal.csv, bounds.csv, verdict.json,
  comparison_table.md, the figure.
- D: the eight wirings at F's held points again with the stage matrices saved (E90's `read_stages`; each run's counts
  checked against part B's); `--part Dscore`: E90's claw readings and verdict by group on all claws, on the claws of the
  46 glomeruli with input and on those of the 36 (stages.csv, the verdict's 'D').
- E: the family fitted again with F (`E10_joint_fit.py --orn-input <F> -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss
  apl --apl-readout both --tag _cap3`, E_MAXFEV evaluations, from each wiring's own fit), on S (stage S), then on T (stage
  T) only if S's fitted population is E_S_GAIN or more below the reference against all recordings; refused unless the
  gate passes or F is an equal (`--force` overrides); `--part Escore`: fits.csv's rows scored, the population of eight by
  E91's part D rule.
- F: only on a pass by (C) or (E) (not written).
- variants: variants.csv and model_pair_distances_long.csv for the archive:
  `./venv/bin/python archive_tools/archive_model_outputs.py E93_missing_glomeruli_input`, then
  `./venv/bin/python archive_tools/migrate_model_output_params.py --apply`.
- table: comparison_table.md and the figure from what exists.

Outputs: results/E93_missing_glomeruli_input/; spike counts and responses under outputs/E93_missing_glomeruli_input/
(B/ for part B, D/ for the stage runs, E's fits at the folder's root as E10_joint_fit.py writes them).

Usage (from the repository root; one 16-thread process at a time on viki):
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part 0
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part B
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part score
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E93_missing_glomeruli_input'

# ---------------------------------------------------------------------------
# the inputs (the entry's Design; fixed before any run)

ADDED = ('VA4', 'VL1', 'DL4', 'VM4', 'DA4m', 'DA3', 'DA4l', 'VA1v', 'DL3', 'VA1d')
LAB_MEGAMAT = ('VA4', 'VL1')            # megamat_orn_deltas(), the 12 megamat odors
LAB_VALIDATION = ('DL4', 'VM4')         # validation2_orn_deltas(), paa, -aPine, euc, B-myr
HALLEM_ONLY = ('DA4m', 'DA3', 'DA4l', 'VA1v', 'DL3', 'VA1d')   # Hallem's 15 odors
HALLEM_MEGAMAT = ('DL4',)               # Hallem for the 12 megamat odors (the validation odors are the lab's)
HALLEM_NAMES = {'2h': '2-heptanone', 'IaA': 'isopentyl acetate', 'pa': 'pentyl acetate', '2-but': '2-butanone',
                'eb': 'ethyl butyrate', 'ep': 'ethyl propionate', '6al': 'hexanal', 't2h': 'E2-hexenal', '1-8ol': '1-octanol',
                '1-5ol': '1-pentanol', '1-6ol': '1-hexanol', 'ms': 'methyl salicylate', 'paa': 'phenylacetaldehyde',
                '-aPine': 'a-pinene', 'B-myr': 'b-myrcene'}
HALLEM_CONC = -2.0
INPUTS = ('I0', 'I0z', 'F', 'L', 'U', 'V')
CONTROL_KINDS = ('PL', 'PO')            # the ten rows' labels permuted; each row's values permuted over the odors
CONTROL_SEEDS = (1, 2)
NAMED_ODORS = ('2h', 'IaA', 'pa', '1-5ol', '1-6ol')
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# part A

SHARE_BINS = (-0.01, 0.0, 0.2, 0.4, 1.0)
SHARE_LABELS = ('none', 'up to 0.2', '0.2 to 0.4', 'above 0.4')
OLFDB_MEGAMAT, OLFDB_VALIDATION = 1, 4  # the lab's ORN datasets (olfdb), megamat and validation2
N_HALF_SPLITS = 200
DEP_FLOOR = 0.05                        # the departure's floor on the relative responses
DEP_NORMS = ('glom', 'odor', 'model_pn')
DEP_BOOT = 1000
NULL_DRAWS = 2000
A2_SUBSETS = ('meets', 'meets_read')    # E86's rule's entries; at 0.6
A2_MIN_NORMS = 2                        # the rule: in at least two of the three normalizations
DD_BOUND = 2.0                          # the extra arm's factors bounded to [1 / DD_BOUND, DD_BOUND]
SEED = 0

# ---------------------------------------------------------------------------
# part C, the gate (fixed in the entry)

REF_ALL8 = 0.0482                       # the reference's eight wirings (E83), pattern part, readout, all recordings
GATE_GAIN = 0.003                       # (i) T's population against all recordings and the mean over the four halves
GATE_CELLS = 12                         # (ii) of the 16 cells (T's four wirings by the four halves)
GATE_RAW_SLACK, GATE_MEAN_SLACK = 0.005, 0.01   # (iii)
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075  # (iv)
CONTROL_RATIO = 0.5                     # (v)
CONTROL_DRAW2 = 0.0015                  # a second draw of each control if F moves T this much either way
BOUND_GAIN = 0.0015                     # a bound is preferred to F only this far below on S against one half
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
HALVES = ('A', 'B', 'C', 'D')
HALF_OTHER = {'A': 'B', 'B': 'A', 'C': 'D', 'D': 'C'}
FLIES_MEAN = 0.178

# ---------------------------------------------------------------------------
# part E

E_MAXFEV = 150
E_STAGES = {'S': S_WIRINGS, 'T': T_WIRINGS}
E_S_GAIN = 0.0015
E_PASS_GAIN = 0.003
E_WIRINGS_LOWER = 6
E_RAW_REF, E_RAW_SLACK, E_MEAN_SLACK = 0.073, 0.005, 0.01


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# the inputs

def lab_input() -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()`: 36 glomeruli x (panel, odor))."""
    raise NotImplementedError


def lab_rows() -> pd.DataFrame:
    """VA4 and VL1 (megamat odors) and DL4 and VM4 (validation odors) from the lab's tables, in the lab input's columns;
    NaN where the panel did not measure the glomerulus."""
    raise NotImplementedError


def hallem_table() -> pd.DataFrame:
    """Hallem's spike deltas at 10^-2 (drosolf, `add_sfr=False`), the uniglomerular receptors' glomeruli x the 15 odors
    of HALLEM_NAMES (abbreviated)."""
    raise NotImplementedError


def slope_through_origin(x, y) -> float:
    """The least-squares slope of y on x through the origin (sum x y / sum x x)."""
    raise NotImplementedError


def hallem_scaling(lab: Optional[pd.DataFrame] = None, hal: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per odor of Hallem's 15 (and 'pooled'): the slope through the origin of the lab's input on Hallem's over the
    glomeruli both have, Pearson r, the number of glomeruli and their largest values (hallem_scaling.csv)."""
    raise NotImplementedError


def added_rows(kind: str = 'F') -> pd.DataFrame:
    """The ten added rows (ADDED x the lab input's columns): 'F' (the lab's rows and Hallem's scaled per odor), 'U'
    (Hallem's as recorded), 'L' (the lab's four rows only), 'V' (VA4's row only); unmeasured entries 0."""
    raise NotImplementedError


def measured_mask() -> pd.DataFrame:
    """Which of F's 160 added entries are measured (bool, ADDED x the lab input's columns; 134 True)."""
    raise NotImplementedError


def full_input(rows: pd.DataFrame, lab: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """The lab's 36 rows followed by `rows` (the lab input's columns)."""
    raise NotImplementedError


def permute_labels(rows: pd.DataFrame, seed: int) -> pd.DataFrame:
    """The rows' labels permuted among themselves with no row left in place (a derangement drawn with `seed`)."""
    raise NotImplementedError


def permute_odors(rows: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Each row's values permuted over the odors (generator `seed`, rows in order), the megamat and the validation odors
    separately."""
    raise NotImplementedError


def input_of(name: str) -> pd.DataFrame:
    """An input by name: INPUTS, 'PL<seed>' / 'PO<seed>' (the controls) or 'DD' (part A2's arm)."""
    raise NotImplementedError


def write_inputs(names: Optional[Sequence[str]] = None) -> dict:
    """inputs/<name>.csv (glomeruli x abbreviated odors) with their sha256 (inputs/sha256.json), added_rows.csv (every
    added entry of F, U, L, V with its source and whether measured), hallem_scaling.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part 0

def identity(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """Part 0 (the module docstring); identity.json; raises on a difference."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A1

def claw_counts(label: str) -> pd.DataFrame:
    """KC x glomerulus claw counts of a wiring (E90's units, archived; else from E83's `wiring_parts`)."""
    raise NotImplementedError


def kc_classes(M: pd.DataFrame) -> pd.DataFrame:
    """Per KC: its claws, its share of claws on the zero-filled glomeruli and on the ten, whether it has a claw on VA4
    and on DA4l or DA3, and the share bins (SHARE_BINS)."""
    raise NotImplementedError


def who_responds(o: Dict[str, tuple], M: pd.DataFrame, **meta) -> List[dict]:
    """Rows of the response rate and the odors per KC by class (`kc_classes`), APL on and off, with the named odors'
    response rates for the VA4 and the DA4l / DA3 classes."""
    raise NotImplementedError


def input_shares(orn: pd.DataFrame, base: pd.DataFrame, claws: pd.Series) -> pd.Series:
    """Per odor the added share of claw-weighted positive input: (sum over glomeruli of claws x max(orn, 0) minus the
    same for `base`) over the former."""
    raise NotImplementedError


def weighted_distances(orn: pd.DataFrame, claws: pd.Series) -> pd.Series:
    """1 - the claw-weighted Pearson correlation over the glomeruli with input, per pair of odors (canonical order)."""
    raise NotImplementedError


def half_agreement(values: np.ndarray, n_splits: int = N_HALF_SPLITS, seed: int = SEED) -> float:
    """The mean over random splits of the flies in two halves of Pearson r between the halves' mean profiles (flies x
    odors, NaN where not measured)."""
    raise NotImplementedError


def part_a1() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A2

def departure(P: pd.DataFrame, O: pd.DataFrame, norm: str, floor: float = DEP_FLOOR) -> pd.DataFrame:
    """log2 of the dendrites' relative response over the ORNs' (odors x glomeruli, both clipped at 0): relative within a
    glomerulus over the odors ('glom'), within an odor over the glomeruli ('odor'), or over the whole table ('global',
    also the model's PN stage's normalization: `O` is then the model's PN rate changes), with `floor` added to both."""
    raise NotImplementedError


def sign_agreement(dep: np.ndarray, factor: np.ndarray) -> float:
    """The share of entries whose departure and searched log2 factor have the same sign (zeros count as disagreement)."""
    raise NotImplementedError


def null_agreement(dep: pd.Series, odors: Sequence[str], factor: np.ndarray, mask: np.ndarray, n: int = NULL_DRAWS,
                   seed: int = SEED) -> np.ndarray:
    """The sign agreement over `mask`'s entries with the departures permuted within each odor, `n` draws."""
    raise NotImplementedError


def a2_rule(results: Sequence[dict]) -> str:
    """The entry's rule over the normalizations' results (dicts with 'agree', 'null_p5', 'null_p95'): "in the list's
    direction" if the agreement is above the null's 95th percentile in at least A2_MIN_NORMS, "opposite" if below its
    5th in at least A2_MIN_NORMS, "neither" otherwise."""
    raise NotImplementedError


def part_a2() -> dict:
    raise NotImplementedError


def part_a() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B

def held93(label: str, orn: pd.DataFrame, name: str) -> dict:
    """E85's `held_reading` with this experiment's cache (the breadth's APL-off run under outputs/E93.../)."""
    raise NotImplementedError


def b_plan() -> List[Tuple[str, str]]:
    """[(input, wiring)]: F on all eight; L, U, V on S; the controls' first draws on T; their second draws if F moves
    T's population by CONTROL_DRAW2 or more (once F's T runs exist); DD on all eight if part A2 says so."""
    raise NotImplementedError


def run_extras(o: Dict[str, tuple]) -> dict:
    """The named odors' response fractions and the share of KCs responding to no odor, APL on and off."""
    raise NotImplementedError


def part_b(inputs: Optional[Sequence[str]] = None, wirings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# scoring

def targets93() -> Dict[str, pd.Series]:
    """E89's `targets89`: halves A to D, all recordings, the six typical, the noise-free effect, baseline-noise units."""
    raise NotImplementedError


def distances(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances: 'raw', 'cal' (E83's), 'zcell' (E73's per-KC reading, graded)."""
    raise NotImplementedError


def unit_scores(d: Dict[str, pd.Series], tg: Dict[str, pd.Series]) -> dict:
    """E83's `effect_scores` plus, per target, the pattern part, the loss and ρ over the 66 megamat pairs."""
    raise NotImplementedError


def spread_rows(d: Dict[str, pd.Series], tg: Dict[str, pd.Series], pairs: str = 'all', **meta) -> List[dict]:
    """Per target, the pattern part with its amplitude and shape parts, the spread ratio, Pearson r and ρ (on the 120
    pairs or the 66 megamat pairs; the target through `H.canonical_pairs`; the parts checked to add up)."""
    raise NotImplementedError


def readings(force: bool = False) -> pd.DataFrame:
    raise NotImplementedError


def gate93(rd: pd.DataFrame, runs: pd.DataFrame, name: str = 'F', ref: str = 'I0') -> dict:
    """The gate (C) on T for input `name` against the reference: (i) to (vi), with S's population and the number of the
    eight wirings lower beside."""
    raise NotImplementedError


def equal93(rd: pd.DataFrame, runs: pd.DataFrame, name: str = 'F', ref: str = 'I0') -> dict:
    """The equal rule for input `name`."""
    raise NotImplementedError


def bounds93(rd: pd.DataFrame, runs: pd.DataFrame) -> List[dict]:
    """The rule for fitted changes on the bounds L, U, V."""
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part D

def part_d(wirings: Sequence[str] = ALL_WIRINGS) -> None:
    raise NotImplementedError


def score_d() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part E

def e_allowed(stage: str) -> Tuple[bool, str]:
    raise NotImplementedError


def part_e(stage: str = 'S', force: bool = False, dry_run: bool = False) -> List[str]:
    raise NotImplementedError


def score_e() -> dict:
    raise NotImplementedError


def part_f() -> None:
    """Part F (only on a pass by (C) or (E)): the firing-ceiling family with the response rate with APL on re-held, S at
    the flies' 2 s pulse with both response rates re-held, and S through the trial noise as E72's part E ran it; each
    holds if the population keeps at least half of its gain. Not written: it runs only on a pass."""
    raise NotImplementedError('E93 part F runs only on a pass by part C or part E; it is not written')


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def comparison_table() -> str:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('inputs', '0', 'A', 'A1', 'A2', 'B', 'score', 'C', 'D', 'Dscore', 'E',
                                                      'Escore', 'F', 'variants', 'table'))
    ap.add_argument('--wirings', default='all8', help="0 / B / D: 'S', 'T', 'all8' or a comma-separated list")
    ap.add_argument('--inputs', default=None, help='B: comma-separated inputs (default: the plan)')
    ap.add_argument('--stage', default='S', choices=list(E_STAGES), help='E: which fits')
    ap.add_argument('--force', action='store_true', help="E: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='E: print the commands only')
    raise SystemExit('E93: stub (the interface); the parts are being written')


if __name__ == '__main__':
    main()
