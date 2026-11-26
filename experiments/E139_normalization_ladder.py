#!/usr/bin/env python3
"""E139, the weaker antennal-lobe normalization without the calyx refit, along its strength: the megamat stages, the
held-out panels and the APL effect at four rungs of the PN normalization's additive term, and the two calyx factors of
E136's A1 one at a time (experiment_priorities.md, E139; planned by a Fable subagent, claude-fable-5-2, spawned from
viki 2026-11-25 18:39).

E136's A1 holds (the megamat stage error on the T wirings against the held-out half falls from the unrefitted family's
1.049 to 0.934) at a point that joins three changes: the additive term at 32 times the family's (1017 against 31.4),
APL->bouton at x0.155 and KC-KC at x3.24. This reads them apart at fixed scales; nothing is fitted first.

Definitions (fixed before any run):
- a rung is a value of pn_inhadd with pn_inhsc following by E77's rule (`inhsc_for`), as in E102's ladder and A1's
  coordinate (E136's `pn_norm_kws` at z_inhsc = 0): RUNGS;
- the settings (SETTINGS): F, the family's own calyx factors; B, APL->bouton at x0.155 with KC-KC at x1; K, KC-KC at
  x3.24 with APL->bouton at x1 (E136's A1 has both);
- a point is a setting at a rung, read on one wiring and one panel's input exactly as E136 reads an arm (`read_point`:
  E136's read at the point's coordinates, both rates held at the panel's operating point: megamat 0.0915, validation2
  0.050, the mixtures 0.054, each with the wiring's own APL-off ratio, the APL16 input at the wiring's own fit's two
  counts), on all eight wirings;
- a setting is read on megamat against the fit half, the held-out half and all flies on all eight wirings; a rung is
  named on S against the fit half; every class is E136's, on T against the held-out half, against the unrefitted
  family; "a price on the APL effect": the calibrated pattern part above the family's by more than PRICE (E102's rule)
  against all recordings or either of halves A and B; "serves both": holds with no price.

Design:
- step 0 (no simulation, `--part step0`): A1's saved evaluations and their rescore binned by rung (`a1_by_rung`); E102's
  saved ladder scored as E136 scores a read (`e102_ladder`); the PN stage alone along the ladder on E135's port
  (`pn_stage_ladder`); the recorded step's signature, each odor's residual off a proportional step against its summed
  ORN input, beside each rung's (`recorded_step_signature`);
- step 1 (`--part read --step 1`): F at the four rungs on the megamat, validation2 and mixture inputs, eight wirings
  each (96 reads); F's APL16 rows are E102's saved runs;
- step 2 (`--part read --step 2`): B at 314 and 1000 and K at 1000 on the megamat, validation2, mixture and APL16
  inputs, eight wirings each (96 reads);
- beside (`--part read --step beside`): B at 1000 and A1 on E98's extended input, read by E102's rules for the flies'
  APL along concentration and the KCs at low concentration (16 reads);
- the score (`--part score`): E136's, with the unrefitted family as the reference, the kiwi mixtures and both mixture
  panels' second-block readings (E137) beside every setting;
- step 3 (`--part fit`, only on a candidate): if a rung below 1000 "improves" at F or B with a price under half of
  A1's, A0's two-factor fit with the normalization fixed at that rung, on S against the fit half, then E136's reads.

Outputs: results/E139_normalization_ladder/ (reads_<host>_<lane>.csv, claims/, step0_*.csv, score_*.csv,
classes.json); runs in data/model_outputs/E139_normalization_ladder/read/ (E134's step-3 layout, written in place, as
E136's).
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E139_normalization_ladder'
SCRIPT = 'experiments/E139_normalization_ladder.py'
# the rungs of pn_inhadd (olfsysm's default 31.4088 is the family's)
RUNGS = (100.0, 314.0, 1000.0, 3141.0)
# the calyx factors of the settings, as multiples of the family's (E136's A1: APL->bouton x0.155, KC-KC x3.24)
SETTINGS: Dict[str, Dict[str, float]] = {'F': {}, 'B': {'aplpn': 0.155}, 'K': {'kckc': 3.24}}
# per step: (setting, rungs) pairs and the panels read (F's APL16 rows are E102's saved runs)
STEPS: Dict[str, tuple] = {
    '1': ((('F', RUNGS),), ('megamat', 'validation2', 'mixtures')),
    '2': ((('B', (314.0, 1000.0)), ('K', (1000.0,))), ('megamat', 'validation2', 'mixtures', 'apl16')),
}
# "a price on the APL effect": the calibrated pattern part above the family's by more than this (E102's rule)
PRICE = 0.003
# step 3's candidate: a rung below this "improves" at F or B with a price under this share of A1's
CANDIDATE_BELOW, CANDIDATE_PRICE_SHARE = 1000.0, 0.5


@dataclass(frozen=True)
class Point:
    """A setting at a rung, read on one wiring and one panel's input."""
    setting: str
    rung: float
    panel: str
    wiring: str


def out_dir() -> Path:
    """<shared>/results/E139_normalization_ladder/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E139_normalization_ladder/."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def point_tag(pt: Point) -> str:
    """'<setting>_inhadd<rung>_<panel>_<wiring>', e.g. 'F_inhadd314_megamat_P3'."""
    if pt.setting not in SETTINGS:
        raise ValueError(f'{pt}')
    return f'{pt.setting}_inhadd{pt.rung:g}_{pt.panel}_{pt.wiring}'


def point_z(setting: str, rung: float) -> Dict[str, float]:
    """The point's coordinates in E136's A1 (log factors of each wiring's own scales; `pn_norm_kws`): aplpn and kckc
    from SETTINGS, inhadd = ln(rung / olfsysm's default), inhsc = 0 (pn_inhsc follows by `inhsc_for`)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    z = {k: float(np.log(v)) for k, v in SETTINGS[setting].items()}
    z.update(aplpn=z.get('aplpn', 0.0), kckc=z.get('kckc', 0.0), inhadd=float(np.log(rung / PN_INHADD_DEFAULT)),
             inhsc=0.0)
    return z


def planned_points(steps: Sequence[str] = ('1', '2')) -> List[Point]:
    """The reads of steps 1 and 2 in the order a forward lane takes them: per step, per setting and rung, per panel,
    the eight wirings."""
    from experiments import E136_stage_population_fit as F
    pts = []
    for s in steps:
        groups, panels = STEPS[s]
        for setting, rungs in groups:
            for r in rungs:
                for p in panels:
                    pts += [Point(setting, float(r), p, w) for w in F.WIRINGS]
    return pts


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read, claimed under out_dir()/claims/ as E136's reads are: E136's read (`read_at`) of A1's model at
    `point_z(pt.setting, pt.rung)` on the point's wiring and panel, both rates held at the panel's operating point, saved
    in the archive's read/; its row in reads_<host>_<lane>.csv."""
    raise NotImplementedError


def part_read(lane: str, steps: Sequence[str], reverse: bool = False, reclaim: bool = False) -> None:
    """A lane of reads: every planned point of `steps`, skipping what is done or claimed elsewhere."""
    raise NotImplementedError


def a1_by_rung() -> pd.DataFrame:
    """Step 0: A1's saved evaluations (fit_history_A1.csv) and their rescore against both halves (rescore_A1.csv),
    binned by the additive term: each stage's error along it, with the error over the boutons, claws and somata alone
    (single evaluations whose other factors differ: a lead)."""
    raise NotImplementedError


def e102_ladder() -> pd.DataFrame:
    """Step 0: E102's saved ladder (its archive's B/inhadd<rung>/, the family's calyx on the APL16 input) scored as E136
    scores a read: the APL16 somata's strength and all-view errors, the APL-effect statistics, rung by rung and wiring by
    wiring; F's APL16 rows of the score."""
    raise NotImplementedError


def pn_stage_ladder() -> pd.DataFrame:
    """Step 0: the PN stage alone along the ladder (E135's numpy port): the dendrites' strength, matrix and sharing
    errors on held-out flies, and the spread of the odors' totals (the family 0.074, the flies 0.312)."""
    raise NotImplementedError


def recorded_step_signature() -> pd.DataFrame:
    """Step 0: each odor's residual off a proportional step from the recorded ORN terminals to the recorded PN
    dendrites against the odor's summed ORN input, beside what the PN stage gives at each rung, after checking that the
    family's and A1's stages can be told apart through the flies' noise."""
    raise NotImplementedError


def part_score() -> dict:
    """The score: every setting at every rung through E136's score (`score_megamat`, `score_panels`, `apl16_rows`,
    `classify`) with the unrefitted family as the reference, F's APL16 rows from `e102_ladder`; the price on the APL
    effect; the E137 readings beside; classes.json."""
    raise NotImplementedError


def part_fit(rung: float) -> dict:
    """Step 3, only on a candidate: A0's two-factor fit (E136's `fit_arm`) with the normalization fixed at `rung`, on S
    against the fit half."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'read', 'score', 'fit'))
    ap.add_argument('--step', nargs='+', default=['1', '2'], choices=('1', '2', 'beside'))
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--rung', type=float)
    a = ap.parse_args(argv)
    if a.part == 'step0':
        for f in (a1_by_rung, e102_ladder, pn_stage_ladder, recorded_step_signature):
            f().to_csv(out_dir() / f'step0_{f.__name__}.csv', index=False)
    elif a.part == 'read':
        part_read(a.lane, a.step, a.reverse, a.reclaim)
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'fit':
        if a.rung is None:
            raise SystemExit('--part fit needs --rung')
        print(json.dumps(part_fit(a.rung), indent=1, default=str))


if __name__ == '__main__':
    main()
