#!/usr/bin/env python3
"""E133: what in a calyx the best family's standing depends on: the male CNS's left calyx against its own re-drawn
wirings, and one difference from the hemibrain moved at a time.

The entry is item E133 of experiment_priorities.md (moved in by viki 2026-11-24 05:51 under the idle-host rule; planned
by a Fable subagent, claude-fable-5-2). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight
wirings (E83's `saved_fit`). "Fixed scales" is E130's step 1: the hemibrain fit's (C2's) factors on the calyx's own base
(E130's bases.json), KC-KC at C2's scale on the calyx's own matrix (the male calyces' from Pratyush's KCclaw-to-KCclaw
table, E130's matrices; the hemibrain's olfsysm's), both response rates re-held (E77_input_screen.hold_rates) at C2's
response counts as fractions of the calyx's KCs x odors. The calyces: the hemibrain under multibouton_claws 'argmax'
and the male CNS's left and right (E10_joint_fit.py's `--connectome`).

STUB (2026-11-24): the parts and their signatures; bodies follow.

Parts (CLI `--part`):
- step0 (no simulation): (i) from E130's saved runs of the three calyces at fixed scales ('S1_<calyx>') and of their
  counted fits, and of the family's eight wirings: the calibrated pattern part over all 120 pairs and over the 39 pairs
  of one of 2-but, paa and B-myr with one of the other 13 odors and the other 81 (`pair_groups`), against all
  recordings and each of the four halves A to D (E108's `effect_targets`); for the 120 also the amplitude and shape
  parts, the spread ratio and Pearson r (E108's `unit_scores`). (ii) The three odors' responders with APL on and off,
  by KC type and claw count, and their overlap. (iii) From the model's tables, for each of the three odors' glomeruli
  (those carrying at least a tenth of the odor's summed positive ORN input) on each calyx: the claws by KC type, the KCs
  with two or more such claws, the boutons, and APL's synapses on those claws and boutons. step0.json, step0_*.csv.
- check (simulations): the identity. For each calyx, with no transplant and no re-draw, the runner (`calyx_tables`,
  `build_simulator`, the step-1 hold from E130's start) reproduces E130's fixed-scale run S1_<calyx>: the same held
  factors and 0 differing spike counts and responses, APL on and off; and `breadth_run` on the calyx's tables
  reproduces the breadth E130 used (the path step 2 takes). check_<host>.json.
- step1 (simulations, resumable, lanes `--calyces`, `--seeds`, `--reverse`; points claimed by mkdir of
  results/E133_calyx_transplants/claims/<tag>/): the 21 re-draws at fixed scales, seeds 1 to 7 of the within-type claw
  permutation (pyolfsysm's `shuffle_wPNKC`, via E10_joint_fit.bouton_base's `wiring='claw_permutation_within_type'`)
  of each calyx, E130's step 1 on each (the breadth recomputed on the permuted wiring, the rest as the calyx's).
  Rows in runs_1_<host>_<lane>.csv, provenance in provenance_1_<host>_<lane>.json, the runs locally
  (outputs/E133_calyx_transplants/step1/) and in the archive (data/model_outputs/E133_calyx_transplants/step1/), and
  each point's raw-count pair distances in results/.../pairs/<tag>.csv.
- score1 (no simulation): L, H and step 1's class; the right calyx's seven beside (no class); readings1.csv,
  runs1.csv, verdict1.json, step1_table.md.
- variants (no simulation): variants.csv and model_pair_distances_long.csv for kcq/apl_effect.py and the archive.
- step2 (simulations; on step 1's "kept" or "mixed"): the six arms on each host's four wirings, through
  experiments/E133_builders.py (imported only here).
- commands: the lanes' command lines.

Readings the entry leaves open, fixed here before any run (the reading closest to E130's code):
- "The 39 pairs that hold 2-but, paa or B-myr" are the pairs of one of the three with one of the other 13 odors
  (3 x 13 = 39); the 3 pairs among the three odors are in "the other 81" with the 78 pairs of neither. The planner's
  table has 39 and 81 and no other split adds to those counts; the 42 / 78 split (any of the three) is reported beside.
- A group's pattern error: the per-pair APL effect's pattern error over all 120 pairs (model minus data, each with its
  mean over the 120 removed, E83's `pattern_error`), root mean square over the group's pairs. So the 120 pairs' pattern
  part squared is the pair-count-weighted mean of the groups' squares.
- The "no class" guard of a step-1 point: its largest spike count per KC and odor above twice its own unpermuted
  calyx's (E130's S1_<calyx>) in either condition; or outside a band (E110's `bands_of`, E130's five); or at a corner (a
  bisection at its bracket's bound, or the held threshold outside 1/4 to 4 times the base's or APL->KC outside 1/30 to
  30 times the start's, E124's bounds read as E130's `bound_flags`).
- Step 1's class is read on 'all' (all recordings), the halves A to D beside. A class needs every one of the 14 runs of
  L and H to have one (E130's `class_population`); otherwise no class, with the would-be class and the runs named.
  "The hemibrain's median" is the median of the hemibrain's seven re-draws.

Holds (CLAUDE.md, "Connectome data"): the hemibrain and Pratyush's male CNS only (olfsysm's model tables and the shared
directory's data/Male-CNS_Tom/). No olfsysm change. Bouton identities are never inferred from position or tuning. No
recording outside olfdb is read. Nothing is published.

Concentrations: the model's input for paa, euc and -aPine is the ORN data within 1 log10 unit of the flies'.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E133_calyx_transplants.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E133_calyx_transplants.py --part check
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E133_calyx_transplants.py --part step1 --calyces malecns-left hemibrain-argmax malecns-right
    ./venv/bin/python -u experiments/E133_calyx_transplants.py --part score1
"""
from __future__ import annotations

import argparse
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E133_calyx_transplants'
SCRIPT = 'experiments/E133_calyx_transplants.py'
CALYCES = ('hemibrain-argmax', 'malecns-left', 'malecns-right')
LEFT, HEMI, RIGHT = 'malecns-left', 'hemibrain-argmax', 'malecns-right'
HOSTS = (HEMI, LEFT)
DONOR = {HEMI: LEFT, LEFT: HEMI}
SEEDS1 = tuple(range(1, 8))
HOST_SEEDS = (0, 1, 2, 3)
THREE = ('2-but', 'paa', 'B-myr')
GLOM_SHARE = 0.1
TARGETS_READ = ('all', 'A', 'B', 'C', 'D')
READOUTS = ('raw', 'cal')
# step 1's class (the entry's, fixed beforehand)
KEPT_GAP = 0.007
REDEALT_GAP = 0.004
KEPT_ABOVE = 5
MAX_COUNT_FACTOR = 2.0
# E124's bounds as E130 reads them
THR_BOUNDS = (0.25, 4.0)
APLKC_BOUNDS = (1 / 30, 30.0)
BOUND_TOL = 0.01

Point = namedtuple('Point', 'step calyx seed arm')     # arm None for step 1


def tag_of(pt: Point) -> str:
    """'S1_<calyx>_seed<n>' (step 1); 'S2_<arm>_<host>_seed<n>' (step 2)."""
    raise NotImplementedError


def pair_groups(pairs: pd.Index, odors: Sequence[str] = THREE) -> Dict[str, np.ndarray]:
    """Boolean masks over canonical pairs: 'p120', 'p39' (one of `odors` with one other odor), 'p81' (the rest),
    'among' (both in `odors`), 'none' (neither), 'p42' (any of `odors`)."""
    raise NotImplementedError


def group_patterns(d: Dict[str, pd.Series], targets: Optional[Mapping[str, pd.Series]] = None) -> dict:
    """Per effect target t: '<t>_pattern_<group>' for the groups of `pair_groups` (the module docstring's reading)."""
    raise NotImplementedError


def part_step0() -> dict:
    raise NotImplementedError


def calyx_tables(c: str, seed: int = 0) -> dict:
    """A calyx's tables at one of its wirings (seed 0 its own; n a within-type claw permutation): 'calyx', 'seed', 'orn',
    'base', 'wnorm', 'W', 'breadth', 'kckc' (the male calyx's matrix; None for the hemibrain: olfsysm's own), 'counts'."""
    raise NotImplementedError


def breadth_run(tables: Mapping) -> pd.Series:
    """The per-KC APL-off breadth of `tables` (E10_joint_fit.bouton_base's APL-off run at the base's parameters)."""
    raise NotImplementedError


def build_simulator(tables: Mapping, params: Mapping[str, float]) -> Tuple:
    """C2's family on `tables` as a fixed-parameter simulator (E130's `build_simulator` on given tables)."""
    raise NotImplementedError


def run_point(pt: Point, lane: str, table: Path, prov_path: Path, cc: dict, reclaim: bool = False) -> Optional[dict]:
    raise NotImplementedError


def part_check(calyces: Sequence[str] = CALYCES) -> dict:
    raise NotImplementedError


def part_step1(calyces: Sequence[str] = CALYCES, seeds: Sequence[int] = SEEDS1, lane: Optional[str] = None,
               reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    raise NotImplementedError


def class_step1(left: Sequence[float], hemi: Sequence[float], all_classed: bool = True) -> dict:
    """Step 1's class from the left's and the hemibrain's seven pattern parts (the entry's rule)."""
    raise NotImplementedError


def part_score1() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def part_step2(*args, **kwargs):
    raise NotImplementedError('E133 step 2 runs on experiments/E133_builders.py, not written yet')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'check', 'step1', 'score1', 'variants', 'step2', 'commands'))
    ap.add_argument('--calyces', nargs='+', default=list(CALYCES), choices=list(CALYCES))
    ap.add_argument('--seeds', nargs='+', type=int, default=list(SEEDS1))
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--lane', default=None)
    ap.add_argument('--reclaim', action='store_true')
    a = ap.parse_args()
    raise NotImplementedError(a.part)


if __name__ == '__main__':
    main()
