#!/usr/bin/env python3
"""E92, parts B to F: APL's inhibition of the KCs as a gain on each claw's drive.

The entry is item E92 of experiment_priorities.md (proposed by hal 2026-11-13 19:04; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-13 18:38). Part A (no simulations) is experiments.E92_steps_like_for_like;
this module holds the parts that simulate. olfsysm's `apl_inh_divisive` (`mp.kc.apl_inh_divisive`) divides each claw's
PN drive x by 1 + max(i, 0), i its APL inhibition (APL->claw weight times APL's activity), in place of the clipped
subtraction max(x - i, 0); `bouton_inh_divisive` (`mp.pn.bouton_inh_divisive`) does the same at the PN boutons. With APL
off both are the subtractive runs (x / (1 + 0) = x - 0), so every APL-off run is the saved fit's.

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`, through E91's `WiringModel91`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the
recordings E77's (`recording_halves()`, seed 0), C / D the second split's (seed 1; E89's `targets89`). Every point is on
the fits as they are (one arm, 'B'); its APL->KC weights are E89's at the family's measure F and the point's exponent
(`apl_kc_weights`; F at exponent 2 gives the fit's weights exactly), and the tied APL scale (APL->KC, KC->APL,
bouton->APL; E83's TIED_ON) is bisected to the wiring's APL-on count (E89's `weighting_reading`, E69's `bisect_factor`)
in APL_BRACKET, wider than E83's because the scale changes units under a divisor.

The reference (`--reference`, default 'C2'): C2's clipped family, E89's saved C1 points at the same weightings
(variants B_F_k2 = the unchanged fits, B_F_k1, B_F_k0.5, B_none) and on T the unchanged fits. If E91's refitted family
passes as a population of eight, the reference is that family ('E91_refits'): a stub, since olfsysm refuses the divisive
form together with `claw_inh_pool`, so that branch compares against the pooled form's points at the refits' parameters,
and E91's runs are not read before its note is posted.

Parts (CLI `--part`):
- B1: `apl_inh_divisive=True` at F with exponents 2, 1, 0.5 and no weighting on S (16 points); runs_B1.csv.
- B2, for the record (E91's outcome (e)): `bouton_inh_divisive=True` at exponent 2, and both divisive forms at exponents
  2 and 1 (12 points); runs_B2.csv.
- B3: only if part A's line (A4) names a variant off the line that no gate has read (at most 8 points on S); a stub.
- smoke: B1's point at exponent 2 on the real wiring alone, into smoke/ (runs_smoke.csv; not runs_B1.csv): the
  simulator's identity at the fit, the APL-off run with the divisive form equal to the saved fit's, the bisection inside
  the bracket, the parameter set in place equal to a simulator built with it, the flags reset after.
- co: the co-response changes under APL of the four nearest ester pairs, the ester and the alcohol groups beside every S
  setting (E89's `coresponse_c` readings 1 and 1b, and the geometric mean of E88's two estimates at L1 and L3 with its
  known error from E89's estimator check); coresponse_B.csv.
- C: the choice per half of both splits (`choose92`: the lowest pattern part of S's population under the calibrated
  readout among B1's settings and the references, inside the bands and not at a bracket's bound); a divisive choice that
  gains CHOICE_MIN_GAIN or more over the reference is run on T and read with E83's gate against the other half in four
  directions (`gate92`, E91's `gate91`); the control on S and T (`control_matrix92`: the reference's subtractive form at
  the chosen weighting with a static inhibition at the soma, olfsysm's `ext_inh_kc`, present only with APL on, equal for
  every odor to each KC's mean loss of drive under the gain at the chosen setting, the APL scale re-held); then `score`.
- score: readings_B.csv (every setting per unit and readout, with the amplitude part, the shape part, the spread ratio
  and Pearson r beside every pattern part), beside_B.csv (part A's readings beside every S setting: the two groups'
  distances with APL on in A2's readings and the place against A4's line, when part A has them), coresponse_B.csv
  (unless --no-co), choice.csv, gate.csv, equal.csv, verdict_B.json, comparison_table_B.md.
- D: the family fitted again, only on a choice by C or an equal (`d_plan`): `E10_joint_fit.py -k <exponent> --aplpn-max 3
  --rate-target ctrl_frac --loss apl --apl-readout both --tag _cap3 --apl-divisive`, 150 evaluations, from each wiring's
  own fit of C2's family with its tied APL scale at B1's re-held factor for that exponent (`--init-params`; `--d-start
  exact` starts at the fit's own scales instead), one process at a time; `--pop S` first, `--pop T` only if S's fitted
  population is D_S_GAIN or more below the reference against all recordings. `--part Dscore` scores them (d_scores.csv).
- E: only on a pass by C or D (a stub). F: the summed form on part A's readings after E91's note is posted (a stub).
- variants: variants.csv and model_pair_distances_long.csv for the archive: `./venv/bin/python
  archive_tools/archive_model_outputs.py E92_steps_like_for_like`, then `./venv/bin/python
  archive_tools/migrate_model_output_params.py --apply`.
- table: comparison_table_B.md.

Every run row carries the absolute parameters (the tied APL scale's factor and the scales as run), the form, the
weighting and the spread of its weights over the responding KCs, the realism (E83's `realism`: the four bands under the
rule of 2026-11-10 16:10, `realism_flag`, with the earlier rule beside, `realism_flag_zcrit`), the largest spike counts,
the divisor 1 + i at the responders' claws during the pulse (10th percentile, median, 90th percentile, largest), the
loss of drive over the responders, and `H.code_commits()` of the process that ran it (which must be clean).

Outputs: results/E92_steps_like_for_like/; spike counts and responses under outputs/E92_steps_like_for_like/B/ (each
run's loss of drive per KC and odor beside its counts, <tag>_loss.parquet); D's fits where E10_joint_fit.py writes them
for this experiment name (fits.csv in the results folder, counts at the outputs folder's root).

Usage (from the repository root; one 16-thread process at a time, two at most on hal):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E92_gain --part B1 --wirings S
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E87_inhibition_form as E87
from experiments import E89_apl_breadth_weighting as E89
from experiments import E91_summed_drive_inhibition as E91

EXPERIMENT = 'E92_steps_like_for_like'
ARM = 'B'

FORMS = {'sub': (False, False), 'div': (True, False), 'bdiv': (False, True), 'both': (True, True)}
B1_SETTINGS = ((('F', 2.0), 'div'), (('F', 1.0), 'div'), (('F', 0.5), 'div'), (('F', 0.0), 'div'))
B2_SETTINGS = ((('F', 2.0), 'bdiv'), (('F', 2.0), 'both'), (('F', 1.0), 'both'))
REFERENCE = E89.REFERENCE
RAW_CRIT = E87.RAW_CRIT
PULSE = E91.PULSE
APL_BRACKET = (0.01, 10.0)
APL_STEP = E83.APL_STEP


def variant_name(w: Tuple[str, float], form: str = 'sub', control: bool = False) -> str:
    """A setting's variant name: 'B_<weighting>' (E89's `wtag`) for the subtractive form, so that E89's references read
    as 'B_F_k2' etc.; '_<form>' for a divisive form ('div', 'bdiv', 'both'); '_ctrl' for the control of the setting
    named before it (which itself runs in the subtractive form)."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    """`variant_name`'s fields back: arm, measure, exponent, weighting, wtag, form (as run: 'sub' for a control),
    form_of (the form in the name), apl_inh_divisive, bouton_inh_divisive, control, part ('ref', 'B1', 'B2', 'C')."""
    raise NotImplementedError


def set_form(mp, form: str) -> None:
    """Set the divisive flags of `form` (FORMS) on a simulator's olfsysm parameters (they stay set for its later
    simulations; `set_form(mp, 'sub')` resets them); refuses a divisive KC form with a pooled share set."""
    raise NotImplementedError


def claw_loss(x: np.ndarray, inh: np.ndarray, kc_div: bool) -> Tuple[np.ndarray, np.ndarray]:
    """(claw state, loss of drive) from the claws' PN drive `x` and APL inhibition `inh`: divisive max(x / (1 + max(i,
    0)), 0), whose loss x - x / (1 + i) = state * i; subtractive max(x - i, 0), loss max(x, 0) - state."""
    raise NotImplementedError


def gain_stats(rv, mp, odors: Sequence[int], responders: np.ndarray, form: str) -> dict:
    """From the last APL-on run of a simulator (APL's activity saved): per KC and odor the window means over PULSE of the
    loss of drive ('L') and of the summed drive before the KCs' inhibition ('D'); the divisor 1 + i (and i) at the
    responders' claws and, with the bouton form, at the boutons feeding them, over the pulse's time steps."""
    raise NotImplementedError


def control_matrix92(loss: pd.DataFrame) -> pd.DataFrame:
    """The control's static inhibition at the soma: per KC its mean over the odors of the loss of drive under the gain
    (`loss`, KCs x odors), the same for every odor."""
    raise NotImplementedError


def part_b1(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B1 (16 points on S)."""
    raise NotImplementedError


def part_b2(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B2 (12 points on S)."""
    raise NotImplementedError


def part_b3(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B3: only if A4 names a variant off the line that no gate has read (a stub)."""
    raise NotImplementedError('E92 part B3 runs only if part A4 names a variant off the line; not written')


def part_smoke(out: Path) -> dict:
    """The smoke run (the module docstring)."""
    raise NotImplementedError


def readings92(out: Path, runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings_B.csv."""
    raise NotImplementedError


def choose92(rd: pd.DataFrame, half: str) -> dict:
    """The choice on S against half `half`."""
    raise NotImplementedError


def gate_direction92(rd: pd.DataFrame, runs: pd.DataFrame, X: str, chosen: Optional[str], ref: str,
                     controls: Sequence[str]) -> dict:
    """One direction of the gate on T, with the parts of the pattern part beside."""
    raise NotImplementedError


def equal_rule92(rd: pd.DataFrame) -> List[dict]:
    """The equal rule (E91's)."""
    raise NotImplementedError


def part_c(out: Path) -> dict:
    """C (the module docstring)."""
    raise NotImplementedError


def score(out: Path, n_workers: int = 4, coresponse: bool = True) -> dict:
    """The readings, the choice, the gate, the equal rule and the verdict."""
    raise NotImplementedError


def d_command(label: str, k: float, start: str = 'reheld') -> List[str]:
    """The E10_joint_fit.py command of one D fit (`--apl-divisive`)."""
    raise NotImplementedError


def part_d(out: Path, pop: str = 'S', force: bool = False, dry_run: bool = False) -> List[str]:
    """D's fits of one population."""
    raise NotImplementedError


def score_d(out: Path, write: bool = True) -> dict:
    """D's scores and the verdict as a population of eight (E91's part D rule)."""
    raise NotImplementedError


def part_e(out: Path) -> None:
    """Part E, only on a pass by C or D (a stub)."""
    raise NotImplementedError('E92 part E runs only on a pass by part C or part D; not written')


def part_f(out: Path) -> None:
    """Part F, the summed form on part A's readings after E91's note is posted (a stub)."""
    raise NotImplementedError("E92 part F reads E91's runs only after its note is posted; not written")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('B1', 'B2', 'B3', 'smoke', 'co', 'C', 'score', 'D', 'Dscore', 'E', 'F', 'variants', 'table'))
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
