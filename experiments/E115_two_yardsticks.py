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
- The residual correlation is Pearson's r between the two sets of rank residuals over the 66 pairs (E61's residuals,
  E106's `rank_residuals`; ranks over the pairs both hold). The verdict reads Remy's somata as correlation distance (E90's
  18 flies; her published reading), the cosine reading beside.
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
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E115_two_yardsticks'
CONDITIONS = ('on', 'off')
SHARED_MIN = 0.30
LEVEL_TOL = 0.03
POOLED_MAX_FALL = 0.02
CLASS_NEAR, CLASS_FAR = 10, 15
COMPARABLE = (0.5, 2.0)
MIN_COMPARABLE = 10
N_BOOT, SEED = 1000, 0
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
CLASSES = ('both sets against the family', "Yang's flies alone", 'neither')


def out_dir() -> Path:
    raise NotImplementedError


def residual_correlation(model: pd.Series, yang: pd.Series, remy: pd.Series) -> float:
    """Pearson r over the pairs all three hold between the model's rank residuals against `yang` and `remy`'s rank
    residuals against `yang` (E106's `rank_residuals`)."""
    raise NotImplementedError


def two_set_rank(yang: pd.Series, remy: pd.Series) -> pd.Series:
    """Per pair (the pairs both hold), the mean of its rank among `yang`'s distances and its rank among `remy`'s."""
    raise NotImplementedError


def pair_classes(model: pd.Series, yang: pd.Series, remy: pd.Series, near: int = CLASS_NEAR, far: int = CLASS_FAR) -> pd.Series:
    """Per pair, one of CLASSES (the rule in the module docstring)."""
    raise NotImplementedError


def class_shares(model: pd.Series, yang: pd.Series, classes: pd.Series) -> pd.DataFrame:
    """E61's contributions to 1 - ρ of the model against `yang`, summed per class (rows 'class') and per odor within
    each class (rows 'odor')."""
    raise NotImplementedError


def reader_row(model: pd.Series, yang: pd.Series, remy: pd.Series) -> Dict[str, float]:
    """One reader's part-A numbers: ρ against Yang's, against Remy's, against the two sets' mean rank, and the residual
    correlation; Remy's flies against Yang's beside."""
    raise NotImplementedError


def resample_part_a(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """The joint resamples of Yang's recordings, Remy's flies and the eight wirings: one row per resample, reader,
    condition and Remy reading."""
    raise NotImplementedError


def input_size(change: pd.DataFrame, base: pd.DataFrame) -> float:
    """||change - base||_F / ||base||_F over the odors and glomeruli both hold."""
    raise NotImplementedError


def fractions(folder: str, tag: str) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the per-odor response fractions (two or more spikes) of a saved run, from the model-output archive's
    response matrices."""
    raise NotImplementedError


def rate_rho(frac: pd.Series, cond: str) -> float:
    """Response rate ρ: Spearman of per-odor response fractions with the flies' in baseline-noise units."""
    raise NotImplementedError


def pooled_rate_rho(fracs: Sequence[pd.Series], cond: str) -> float:
    """Response rate ρ of the mean of several runs' per-odor response fractions."""
    raise NotImplementedError


def percentile_of(x: float, values: Sequence[float]) -> float:
    """The percentage of `values` below `x` (ties counted half)."""
    raise NotImplementedError


def price_verdict(directions: Mapping[str, Mapping[str, float]], controls: pd.DataFrame) -> dict:
    """"The price is inside the controls" from each direction's fall, size and pooled fall, and the controls' changes and
    sizes."""
    raise NotImplementedError


def verdicts(a: pd.DataFrame, boot: pd.DataFrame, price: Mapping) -> dict:
    """The three verdicts and the outcome letters."""
    raise NotImplementedError


def part0() -> dict:
    raise NotImplementedError


def part_a() -> Dict[str, pd.DataFrame]:
    raise NotImplementedError


def part_b() -> Dict[str, pd.DataFrame]:
    raise NotImplementedError


def registry_entries() -> List[Tuple[str, float, str, Optional[Sequence[float]]]]:
    """The two yardsticks for kcq/data_targets.py, read from results/E115_two_yardsticks/registry.csv."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('0', 'A', 'B', 'verdict', 'all'), default='all')
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
