#!/usr/bin/env python3
"""E134: the recorded populations of every stage and panel, in one interface, with the lock.

A population is the table of trial-mean responses of all units of one stage (glomeruli, bouton regions, claw regions,
KCs) to the odors of one panel, one table per fly, in baseline-noise units (CLAUDE.md, "The main goal for new
experiments"). This module loads every population E134's entry names (item 2 of experiment_priorities.md, design
steps 0 and 1) and writes and checks the lock: the files with their sha256, the flies, the odor names, the exclusions,
each stage's response criterion, the breadth bins scaled to the panel's size, and for the ORN stage that the recordings
are of axon terminals.

The populations (`POPULATIONS`; keys (stage, panel, cohort)):
- megamat (17 odors at 10^-3): ORNs, PN dendrites, PN boutons, KC claws, KC somata; and the one fly recorded at claws
  and somata on the same trials (cohort 'same_fly');
- validation2 (22 odors): ORNs, KC somata;
- the kiwi / control mixtures: ORNs, KC somata;
- the APL16 panel: KC somata of the control flies (cohort 'control') and of APL output blocked (TNT) ('tnt').

Every reader (hal's kcq/population_fit.py, atlas's kcq/stage_readings.py, experiments/E134_flies.py) goes through
`load_population` / `load_populations`, which refuse to load once the lock exists and a file's sha256 differs from it.

Usage (from the repository root):
    ./venv/bin/python -u data_prep/stage_population_data.py --copy     # locate and copy what the shared dir lacks
    ./venv/bin/python -u data_prep/stage_population_data.py --lock     # write the lock (results/ and experiments/)
    ./venv/bin/python -u data_prep/stage_population_data.py --check    # every population loads and matches the lock
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

STAGES = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PANELS = ('megamat', 'validation2', 'mixtures', 'apl16')
# what one unit of each stage is; only glomeruli are named (the fourth view, "the matrix itself", reads them)
UNIT_KIND = {'orn': 'glomerulus', 'pn_dendrite': 'glomerulus', 'pn_bouton': 'bouton_region',
             'kc_claw': 'claw_region', 'kc_soma': 'kc'}
NAMED_UNITS = {s: UNIT_KIND[s] == 'glomerulus' for s in STAGES}
# (stage, panel, cohort): cohort 'wild_type' (Remy's and the lab's flies), 'same_fly' (the one megamat fly recorded at
# claws and somata on the same trials), 'control' / 'tnt' (Yang's APL16 flies)
POPULATIONS: Tuple[Tuple[str, str, str], ...] = (
    ('orn', 'megamat', 'wild_type'),
    ('pn_dendrite', 'megamat', 'wild_type'),
    ('pn_bouton', 'megamat', 'wild_type'),
    ('kc_claw', 'megamat', 'wild_type'),
    ('kc_soma', 'megamat', 'wild_type'),
    ('kc_claw', 'megamat', 'same_fly'),
    ('kc_soma', 'megamat', 'same_fly'),
    ('orn', 'validation2', 'wild_type'),
    ('kc_soma', 'validation2', 'wild_type'),
    ('orn', 'mixtures', 'wild_type'),
    ('kc_soma', 'mixtures', 'wild_type'),
    ('kc_soma', 'apl16', 'control'),
    ('kc_soma', 'apl16', 'tnt'),
)
# the bouton fly that every analysis leaves out (Remy's dashboards do; data/remy_per_roi/README.md)
EXCLUDED_FLIES = (('pn_bouton', 'megamat', '2022-10-21__fly02__megamat1_001__dsub3'),)
SHARED = Path('/local/tensor/Tom/claude_kc_remapping')
RESULTS = SHARED / 'results' / 'E134_population_baseline'
LOCK_PATH = RESULTS / 'lock.json'                       # the lock, in the shared directory ...
HERE = Path(__file__).resolve().parents[1]
LOCK_COPY = HERE / 'experiments' / 'E134_lock.json'     # ... and its committed copy (no model is read before it)


@dataclass
class Population:
    """One stage's recorded population on one panel: per fly, a units x odors table of trial means.

    Attributes:
        stage, panel, cohort: the key (one of `POPULATIONS`).
        odors: the panel's odors (abbreviations), in panel order; every table's columns are these, NaN where a fly
            was not given the odor.
        flies: the flies (recordings) kept, in order; excluded flies are not loaded.
        tables: fly -> DataFrame, index the units (glomerulus names where `named_units`, else the recording's ROI or
            cell numbers as strings), columns `odors`, values the trial-mean response in `units`.
        trials: fly -> DataFrame (index the units, columns a MultiIndex (odor, repeat)) of single-trial responses in
            the same units, or None when not loaded or not on hand (`load_population(trials=True)`).
        criterion: the response criterion in `units`: a unit responds to an odor if its trial mean is above it.
        units: what the values are, e.g. 'baseline-noise SD' (the window mean minus the pre-odor mean, over the SD of
            the frames in the pre-odor windows), or 'binary' where only responses are on hand.
        graded: False when only binary responses are on hand (the views on how strongly cannot be read).
        unit_kind: `UNIT_KIND[stage]`.
        named_units: True for glomeruli.
        recording_site: where the indicator was imaged ('axon terminals' for the ORNs, 'dendrites', 'boutons',
            'claws', 'somata').
        window: the response window in s after odor onset, and the pre-odor window, as (start, end) pairs.
        exclusions: what was left out (flies, units, odors, trials), one dict per item with a reason.
        files: the files read, as recorded in the lock (path, bytes, sha256).
        notes: anything a reader must know (concentration substitutions, how a fly's odors were pooled).
    """
    stage: str
    panel: str
    cohort: str
    odors: Tuple[str, ...]
    flies: Tuple[str, ...]
    tables: Dict[str, pd.DataFrame]
    criterion: float
    units: str
    graded: bool = True
    trials: Optional[Dict[str, pd.DataFrame]] = None
    unit_kind: str = ''
    named_units: bool = False
    recording_site: str = ''
    window: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    exclusions: List[dict] = field(default_factory=list)
    files: List[dict] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def key(self) -> Tuple[str, str, str]:
        return (self.stage, self.panel, self.cohort)

    def responses(self, fly: str) -> pd.DataFrame:
        """The fly's units x odors table of responses (1.0 above the criterion, 0.0 at or below, NaN where the
        odor was not given)."""
        raise NotImplementedError

    def breadth_bins(self) -> Tuple[Tuple[int, int], ...]:
        """`breadth_bins(len(self.odors))`."""
        return breadth_bins(len(self.odors))

    def n_units(self) -> pd.Series:
        """Units per fly."""
        return pd.Series({f: len(t) for f, t in self.tables.items()}, name='n_units')

    def fly_odors(self, fly: str) -> Tuple[str, ...]:
        """The odors the fly was given (columns not all NaN)."""
        t = self.tables[fly]
        return tuple(o for o in self.odors if t[o].notna().any())


def population_keys() -> Tuple[Tuple[str, str, str], ...]:
    """Every (stage, panel, cohort) this module loads (`POPULATIONS`)."""
    return POPULATIONS


def panel_odors(panel: str) -> Tuple[str, ...]:
    """The panel's odors, as abbreviations, in panel order (megamat 17; validation2 22; the mixtures panel's
    odors and mixtures as Remy names them; APL16 16, pfo dropped)."""
    raise NotImplementedError


def response_criterion(stage: str, panel: str, cohort: str = 'wild_type') -> float:
    """The stage's response criterion in its units (baseline-noise SDs): 1.9684 for the KC somata of the APL16
    panel (the threshold that keeps the control flies at 0.068; notes/response_criterion_decision.md), 2.0 for the
    per-ROI boutons and claws (E119's verdict criterion), and the lock's value for the others."""
    raise NotImplementedError


def breadth_bins(n_odors: int) -> Tuple[Tuple[int, int], ...]:
    """The bins of the number of odors a unit responds to, scaled to the panel's size: 0, 1, 2-3, 4-7, 8-16 on 16
    odors (the realism band's bins), with the upper three edges at 2/16, 4/16 and 8/16 of the panel, rounded and
    kept distinct. Returns closed (lo, hi) pairs covering 0 to n_odors."""
    raise NotImplementedError


def breadth_labels(n_odors: int) -> Tuple[str, ...]:
    """Labels of `breadth_bins(n_odors)` ('0', '1', '2-3', ...)."""
    raise NotImplementedError


def load_population(stage: str, panel: str, cohort: str = 'wild_type', *, trials: bool = False,
                    check_lock: bool = True) -> Population:
    """One recorded population (a key of `POPULATIONS`), excluded flies left out.

    Args:
        stage, panel, cohort: the key.
        trials: also load the single-trial tables (`Population.trials`) where they exist.
        check_lock: once the lock exists, refuse (ValueError) if a file read differs from its sha256 in the lock.
    """
    raise NotImplementedError


def load_populations(keys: Optional[Sequence[Tuple[str, str, str]]] = None, **kwargs) -> Dict[Tuple[str, str, str],
                                                                                                Population]:
    """`load_population` for every key (all of `POPULATIONS` by default); a population whose files are missing is
    left out, and listed by `missing_populations`."""
    raise NotImplementedError


def missing_populations() -> Dict[Tuple[str, str, str], str]:
    """Keys whose files are not on hand, each with what would be needed (outcome (f) of the entry)."""
    raise NotImplementedError


def lock_contents() -> dict:
    """The lock as a dict: per population the files (path, bytes, sha256), flies, odors, exclusions, criterion, units,
    window, breadth bins, recording site, and the module's commit."""
    raise NotImplementedError


def write_lock(path: Path = LOCK_PATH, copy: Path = LOCK_COPY) -> dict:
    """Write `lock_contents()` to the shared directory's lock and its committed copy in experiments/."""
    raise NotImplementedError


def read_lock(path: Path = LOCK_COPY) -> Optional[dict]:
    """The committed lock, or None before it exists."""
    raise NotImplementedError


def require_lock() -> dict:
    """The committed lock; raises RuntimeError if it is missing or uncommitted (no model is read before it)."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    """--copy, --lock, --check (see the module docstring)."""
    raise NotImplementedError


if __name__ == '__main__':
    main()
