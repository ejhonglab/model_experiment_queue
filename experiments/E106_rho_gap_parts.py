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
   `kcq.apl_effect` (columns `spike_counts_cosine__dist_rho_<on|off>_<subset>`); this script reads them back;
2. the family and the ceiling family: the split between the conditions (each condition rebuilt from the other
   condition's distances and an APL effect, the model's or the flies'), E61's rank residuals by pair and by odor, the
   odor part and the pair part of each matrix with both hybrids, each odor left out, and each odor's pairs (and the 29
   pairs of paa and B-myr) put at the data's values;
3. the model's ORN input read the same way, as Pearson and as cosine distance;
4. the yardsticks, as registry entries (`registry_entries`, read by kcq/data_targets.py): the ceilings by subset of
   pairs, one recording against the mean of the others, Remy's wild-type KC somata against both matrices, the SD over
   wirings, the resampling SD and the split-half difference of the family's ρ's, the paired SD of a near-null change,
   and the agreement of two populations of E85's model flies whose inputs come from disjoint ORN flies;
5. the trade-off table: the two ρ's against the APL-effect loss, its parts and the bands' quantities over the archive;
   the inputs searched for the pattern part (E77, E83, E86) read for the two ρ's on their held-out wirings;
6. robustness: against the distances without the trial noise (E81) and against the six typical recordings with APL
   output blocked (TNT) (E76).

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no table of Yang's own Claude's connectome project, no
FlyWire / FAFB / BANC data and no raw-only recording is read. Archived variants whose wirings were built from such
tables (`kcq.apl_effect.HELD_VARIANTS`: every E97 variant, and E4's two FAFB wirings) get no new reading and are left
out of every table here.

Outputs: <shared>/results/E106_rho_gap_parts/ (archive.csv, archive_summary.csv, conditions.csv, residuals_pairs.csv,
residuals_odors.csv, odor_pair_parts.csv, odors.csv, input.csv, tradeoff.csv, searched_inputs.csv, robustness.csv,
yardsticks.csv, comparison_table.md, provenance.json).
Usage: ./venv/bin/python -u experiments/E106_rho_gap_parts.py
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E106_rho_gap_parts'
CONDITIONS = ('on', 'off')


# ---------------------------------------------------------------------------
# pairs, subsets and the splits' algebra (made-up matrices in the tests)

def subset_pairs(index: pd.Index, subset: str) -> pd.Index:
    """The pairs of `index` in `subset` (`kcq.apl_effect.PAIR_SUBSETS`)."""
    raise NotImplementedError


def rho(model: pd.Series, data: pd.Series, subset: str = 'all120') -> float:
    """ρ on odor-odor distances of `model` against `data` over the pairs of `subset`."""
    raise NotImplementedError


def rank_residuals(model: pd.Series, data: pd.Series) -> pd.DataFrame:
    """Per pair (E61): the ranks over the pairs, the rank residual d (model minus data) and the contribution
    6 d^2 / (n (n^2 - 1)) to 1 - ρ."""
    raise NotImplementedError


def condition_split(m_on: pd.Series, m_off: pd.Series, d_on: pd.Series, d_off: pd.Series) -> Dict[str, float]:
    """The split between the conditions: ρ with APL on of the flies' APL-off distances plus the model's per-pair APL
    effect, and of the model's APL-off distances plus the flies' effect (and the same for APL off); ρ between the two
    conditions' rank residuals."""
    raise NotImplementedError


def odor_pair_parts(pairs: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """A pair Series split into its odor part (the least-squares fit of value(a, b) = s_a + s_b over the pairs,
    `kcq.harness.odor_shares`) and its pair part (the rest): (odor part per pair, pair part per pair, s per odor)."""
    raise NotImplementedError


def odor_pair_hybrids(model: pd.Series, data: pd.Series) -> Dict[str, float]:
    """ρ against `data` of the model, of the flies' odor part plus the model's pair part, and of the model's odor part
    plus the flies' pair part; and of the model's odor part alone and of its pair part alone."""
    raise NotImplementedError


def leave_one_odor_out(model: pd.Series, data: pd.Series) -> pd.Series:
    """Per odor: ρ over the pairs without the odor."""
    raise NotImplementedError


def at_data_values(model: pd.Series, data: pd.Series, odors: Sequence[str]) -> pd.Series:
    """The model's distances with every pair involving one of `odors` put at the data's value."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# reading the archive

def load_pairs(folder: str, tag: str) -> Dict[str, pd.Series]:
    """{'on', 'off'}: a saved variant's pair distances (raw counts, cosine; results/<folder>/model_pair_distances_long.csv)."""
    raise NotImplementedError


def family_members(name: str = 'family') -> List[Tuple[str, str, str]]:
    """[(wiring label, results folder, tag)] of the family ('family', eight wirings) or the ceiling family ('ceiling',
    four)."""
    raise NotImplementedError


def population(members: Sequence[Tuple[str, str, str]]) -> Dict[str, pd.Series]:
    """{'on', 'off'}: the mean over the members of their pair distances per condition (E72's population)."""
    raise NotImplementedError


def data_matrices() -> Dict[str, pd.Series]:
    """{'on': the control flies' mean matrix, 'off': APL output blocked (TNT)'s} (`kcq.harness.yang_mean_rdm`)."""
    raise NotImplementedError


def input_matrices() -> Dict[str, pd.Series]:
    """The model's ORN input read as distances: {name: pair Series} for the registry's ORN reference (1 - Pearson of the
    olfdb fly-mean vectors), the cosine of those vectors, and 1 - Pearson and 1 - cosine of the model's own input
    (`kcq.harness.apl16_orn_deltas`)."""
    raise NotImplementedError


def archive(by_variant: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """by_variant.csv's rows with the two ρ's of step 1, the held variants out, one row per distinct variant (the same
    saved distances in several folders counted once)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the tables

def family_tables(name: str = 'family') -> Dict[str, pd.DataFrame]:
    """Step 2 for one family: conditions, residuals by pair and odor, odor and pair parts, odors (left out, at the
    data's values)."""
    raise NotImplementedError


def input_table() -> pd.DataFrame:
    """Step 3: the input's ρ's by subset and its rank residuals against the family's."""
    raise NotImplementedError


def tradeoff_table(arch: pd.DataFrame) -> pd.DataFrame:
    """Step 5: rank correlations over the archive (inside the bands and all) between the two ρ's and the APL-effect
    loss, its parts and the bands' quantities."""
    raise NotImplementedError


def searched_inputs_table() -> pd.DataFrame:
    """Step 5: the two ρ's of the inputs searched for the pattern part (E77, E83, E86) on their held-out wirings, paired
    with the lab's input on the same wirings."""
    raise NotImplementedError


def robustness_table() -> pd.DataFrame:
    """Step 6: the family's ρ's against the noise-free distances (E81) and the six typical TNT recordings (E76)."""
    raise NotImplementedError


def registry_entries() -> List[Tuple[str, float, str, Optional[Sequence[float]]]]:
    """Step 4's yardsticks for kcq/data_targets.py: [(name, value, definition, 90 % interval or None)]."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--out', type=Path, default=None)
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
