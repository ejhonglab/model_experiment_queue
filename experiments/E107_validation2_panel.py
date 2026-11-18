#!/usr/bin/env python3
"""E107: the four validation2 odors in the panel their ORN input was recorded in: the best family against Remy's
validation2 KC recordings. A diagnostic: it changes no model and moves no loss.

The entry is item E107 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08; moved into the list 2026-11-17 16:34). In Yang's panel the family
holds paa and B-myr 0.50 apart with APL on where the control flies have 0.93; this reads Remy's wild-type KC somata of
the validation2 panel (22 odors, every one recorded in the same ORN session as the model's input for paa, -aPine, euc
and B-myr) to say whether the gap is paa's concentration (a), the family's KC stage (b), the ORN recordings (c) or the
joining of two ORN sessions (d). "The family" is the breadth-2 APL-effect fit (C2, E67) on its eight wirings (E72's
part B: the hemibrain wiring and seven re-drawn), each at its own fit (E83's `saved_fit`, the absolute parameters of
E72's / E67's fits.csv), read as the mean matrix over the wirings (E72's population).

Parts (CLI `--part`):
- 0, the lock (no distance of Remy's opened): the files with their sha256, the map of Remy's odor names onto the ORN
  input's (olfdb's odors table), the odors of each set of pairs, the flies behind each pair by concentration, the least
  number of flies for a pair to enter a per-fly statistic, and every choice the entry leaves open; written to
  experiments/E107_part0_lock.json (LOCK_PATH), which is committed before any reading code runs (`require_lock`).
- data: the files read, copied to data/remy_validation2_kc/ in the shared directory with a README (sources, sizes,
  sha256); the 41 MB per-trial file is not copied.
- A (no simulations): Remy's validation2 KC matrix (her mean table and the per-fly matrices), the validation2 ORN input
  alone against it (ρ on odor-odor distances over each set of pairs and over the pairs of paa, B-myr, euc and -aPine),
  the distance of paa from B-myr, euc's and -aPine's pairs by concentration, the ceilings from halves of the flies, her
  megamat somata as the yardstick and the level check against Yang's control flies.
- B (simulations at fixed parameters): the family's eight wirings on the validation2 ORN input and on the megamat ORN
  input (`panel_input`), APL on and acutely off, as fitted (`run_fitted`) and with the response rates re-held
  (`run_held`); each as-fitted run is checked against the saved fit on the odors both panels share (0 count
  differences); resumable; runs under outputs/E107_validation2_panel/B/.
- C (no simulations): the runs read against Remy's matrices (ρ over each set of pairs, E61's rank residuals and shares
  by odor, the pairs of paa and B-myr), the outcome class (`classify`), the realism numbers with no band applied,
  verdict.json, comparison_table.md, the figures.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.

Outputs: <shared>/results/E107_validation2_panel/; spike counts and responses under outputs/E107_validation2_panel/.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part 0
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part data
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part A
    OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 nice -n 10 ./venv/bin/python -u experiments/E107_validation2_panel.py --part B
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part C
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E107_validation2_panel'
LOCK_PATH = H.HERE / 'experiments' / 'E107_part0_lock.json'
VERDICT_PAIR = ('paa', 'B-myr')
FOUR = ('paa', 'B-myr', 'euc', '-aPine')
CLOSE, APART = 0.70, 0.85          # the entry's bands for the distance of paa from B-myr
ORDER_TOL = 0.05                   # "within 0.05 of its megamat level"
LEVEL_TOL = 0.05                   # her megamat somata against Yang's control flies, difference of the means


def remy_dir() -> Path:
    """Remy's folder on tensor (`Remy/MB ms FA2023/data_for_autocheck/`), read only."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part 0, the lock (manifests and stimulus table only)

def file_table(paths: Sequence[str]) -> List[dict]:
    """[{'path', 'size', 'sha256'}] of Remy's files (paths relative to `remy_dir`)."""
    raise NotImplementedError


def olfdb_odors() -> pd.DataFrame:
    """olfdb's odors table (http://hal:8080/api/odors): name, abbrev, aliases, cid."""
    raise NotImplementedError


def odor_map(names: Sequence[str], olfdb: pd.DataFrame, input_names: Sequence[str]) -> Dict[str, str]:
    """{Remy's name: the ORN input's name}, both resolved to one olfdb odor by abbreviation or alias (case-sensitive
    first, then case-insensitive); raises when a name resolves to no odor or to several."""
    raise NotImplementedError


def stimulus_table(panel: str = 'validation2') -> pd.DataFrame:
    """Remy's stimulus table (stims_per_fly_long.csv) for the KC somata of `panel`: fly, stim, odor (mapped), conc,
    presented, in_rdm (after her quality control), qc_dropped."""
    raise NotImplementedError


def pair_sets(stims: pd.DataFrame) -> Dict[str, List[str]]:
    """The odors of each set of pairs ('231', '153', '136', and the four odors' pairs within each), by the lock's rules."""
    raise NotImplementedError


def flies_by_pair(stims: pd.DataFrame, odors: Sequence[str]) -> pd.DataFrame:
    """Per pair and per combination of the two odors' concentrations: the flies whose matrix holds the pair (both
    stimuli in_rdm), their number."""
    raise NotImplementedError


def lock_definitions() -> dict:
    """Every choice of the lock, as a JSON-able dict (the measures, the sets, the minimum flies, the readings of the
    family, the outcome rules)."""
    raise NotImplementedError


def write_lock(force: bool = False) -> dict:
    """Part 0: writes LOCK_PATH (the definitions, the files and their sha256, the odor map, the sets, the flies by pair)
    with its own sha256; refuses to overwrite unless `force`."""
    raise NotImplementedError


def require_lock() -> dict:
    """The committed lock, after checking that LOCK_PATH is committed and unchanged in git and that Remy's files still
    have the locked sha256; raises otherwise (the reading code calls this before opening any distance)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the data pull

def copy_data() -> Path:
    """Copies the files read to <shared>/data/remy_validation2_kc/ (not the per-trial file) with a README.md of source
    paths, sizes and sha256; returns the folder."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# reading Remy's matrices

def mean_table(measure: str = 'cosine', cells: str = 'rmap_nls_cluster') -> pd.Series:
    """Remy's mean table (figure-05/resampled/tidy_abbrev_dists.tsv), column `kc-<measure>-<cells>`, as a canonical
    pair Series over the 22 odors (the input's names; concentrations pooled as in her table)."""
    raise NotImplementedError


def per_fly(measure: str = 'cosine', cells: str = 'rmap_nls_cluster', by_conc: bool = False) -> pd.DataFrame:
    """Remy's per-fly matrices (xrds_stim_rdm_concat.nc with her quality-control drop list, the lock's) as flies x
    pairs (NaN where a fly lacks a pair); with `by_conc`, euc and -aPine keep their concentration in the odor name."""
    raise NotImplementedError


def fly_mean(rows: pd.DataFrame, min_flies: int) -> pd.Series:
    """The mean over flies of each pair held by at least `min_flies` flies."""
    raise NotImplementedError


def split_half_ceiling(rows: pd.DataFrame, pairs: Optional[Sequence[Tuple[str, str]]] = None, n_splits: int = 1000,
                       seed: int = 0, min_flies: int = 1) -> Dict[str, float]:
    """The ceiling of ρ on odor-odor distances from random halves of the flies (rows), each split read on the pairs both
    halves hold (at least `min_flies` flies each): the mean split-half ρ, its 5th and 95th percentiles, the
    Spearman-Brown value for all flies and its square root (the ceiling), as `hong2p.stats.split_half_reliability`."""
    raise NotImplementedError


def megamat_somata(measure: str = 'cosine') -> Tuple[pd.DataFrame, dict]:
    """Remy's megamat KC somata per fly (megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc): flies x pairs, and what was
    read (the measure and cell selection the file holds, by the lock's order of preference)."""
    raise NotImplementedError


def level_check() -> dict:
    """Her megamat somata against Yang's control flies on the megamat pairs both have, same measure: the two means and
    their difference, flagged above LEVEL_TOL."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the ORN input and the measures

def panel_input(panel: str) -> pd.DataFrame:
    """The model's ORN input for `panel` ('validation2': `validation2_orn_deltas`, 22 odors; 'megamat':
    `megamat_orn_deltas`, 17 odors) on the 36 glomeruli of the lab's input (`H.apl16_orn_deltas`), in its layout."""
    raise NotImplementedError


def input_pairs(panel: str, measure: str = 'pearson') -> pd.Series:
    """The ORN input alone as odor-odor distances: 1 - Pearson across the 36 glomeruli (the program's reading;
    'cosine' beside)."""
    raise NotImplementedError


def pairs_with(index: pd.Index, odors: Sequence[str]) -> pd.Index:
    """The pairs of `index` that involve at least one of `odors`."""
    raise NotImplementedError


def rho_over(model: pd.Series, data: pd.Series, pairs: Optional[Sequence[Tuple[str, str]]] = None) -> Tuple[float, int]:
    """ρ on odor-odor distances of `model` against `data` over `pairs` (all common pairs if None): (ρ, number of pairs)."""
    raise NotImplementedError


def odor_shares(model: pd.Series, data: pd.Series, odors: Sequence[str]) -> pd.DataFrame:
    """E61's shares by odor over the pairs of `odors`: per odor the mean absolute rank residual of its pairs, their
    share of 1 - ρ and the change in ρ when its pairs are left out (E61's `residuals`)."""
    raise NotImplementedError


def distance_band(d: float) -> str:
    """'close' (CLOSE or less), 'apart' (APART or more) or 'middle'."""
    raise NotImplementedError


def classify(d_cos: float, d_corr: float, orn: Dict[str, float], fam: Dict[str, float],
             level_flag: bool = False) -> dict:
    """The outcome class by the entry's rules (the lock's `outcome_rules`): from the distance of paa from B-myr in both
    measures and the ρ's of the ORN input alone and of the family (keys 'pb', over the pairs of paa and B-myr; 'all',
    over the carrying set; 'megamat', the megamat level); returns the class, the band used and each condition."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B, the runs

def run_fitted(label: str, panel: str) -> Dict[str, tuple]:
    """One wiring at its fit (E98's `ExtendedRun`, breadth held at the lab input's) on `panel`'s input, APL on and
    acutely off: {'on', 'off'}: (responses, spike counts)."""
    raise NotImplementedError


def run_held(label: str, panel: str, rate_on: float) -> dict:
    """One wiring on `panel`'s input with the response rates re-held (E77's `hold_rates` at the fit's APL-off fraction
    and `rate_on` with APL on, breadth held at the lab input's): hold_rates' dict."""
    raise NotImplementedError


def part_b(panels: Sequence[str] = ('validation2', 'megamat'), wirings: Optional[Sequence[str]] = None,
           readings: Sequence[str] = ('fitted', 'held')) -> None:
    """Runs (and saves) every wiring, panel and reading not yet saved."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part C, the tables

def family_pairs(panel: str, reading: str) -> Dict[str, object]:
    """The family's APL-on pair distances (cosine on spike counts) on `panel`: the population mean over the eight
    wirings and each wiring's."""
    raise NotImplementedError


def realism(panel: str, reading: str) -> pd.DataFrame:
    """The family's response rates with APL on and off, their ratio and its distribution of odors per KC on `panel`,
    per wiring and for the population, with Remy's per-odor response fractions where they exist; no band applied."""
    raise NotImplementedError


def verdict() -> dict:
    """verdict.json: the outcome class and the numbers it rests on, with every partner reading."""
    raise NotImplementedError


def comparison_table() -> str:
    """comparison_table.md: the ORN input alone, the family as fitted and re-held, Remy's flies (ceilings) side by
    side."""
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv for `archive_tools/archive_model_outputs.py`."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'data', 'A', 'B', 'C', 'variants'))
    ap.add_argument('--force', action='store_true')
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
