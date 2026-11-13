#!/usr/bin/env python3
"""E87: the form of APL's inhibition, read from the flies' own responses.

Which inhibition, applied to the KC responses with APL output blocked (TNT), gives the control flies' distances; what
the same inhibition makes of the models' responses without APL; and a target for those responses that needs no APL
(experiment_priorities.md, E87; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-12 20:10).

A recording's responses are its KC x 16-odor trial means (the solvent left out). A response is a value above the
criterion (1 z; the registry's noise_threshold in baseline-noise units); values not above it are set to 0 on both sides
(the registry's thresholded effect). A form reduces KC k's response to odor o by c * A_o^q * B_k^p, A_o the odor's mean
(thresholded) response over the recording's KCs with APL output blocked relative to its mean over odors, B_k the KC's
breadth with APL output blocked relative to its mean over responding KCs (the number of odors it responds to, or its
summed response). Subtractive: the reduced value is kept if it stays above the criterion; divisive: the criterion is
raised by the same amount and surviving values are unchanged. c is set so that the mean response rate over the
recordings is the control flies'. The predicted control matrix is the mean over the recordings with APL output blocked
of the cosine distances of the reduced responses; the predicted effect is that minus the same recordings' own.

Parts (CLI `--part`):
- A: the forms on the flies' responses (no simulations): A1 the grid, the per-odor oracle and the form from KCs of
  similar tuning, as trial means, with the trial noise removed and binarized, in z and in baseline-noise units, with
  floors and permuted controls; A2 the form chosen on one half of the recordings and scored on the other.
- B: the models' APL-off responses through the same forms (no simulations); the two-by-two.
- C: the effect of a form within the recordings with APL output blocked as a target that needs no APL; the models and
  the archive against it.
- D1: the forms as a static inhibition in place of APL on the eight wirings (simulations); D2 only if C3 has a variant
  that stands out.
Outputs: results/E87_inhibition_form/ (comparison_table.md first); large outputs under outputs/E87_inhibition_form/.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E87_inhibition_form'

# the grid of A1 and the thresholds, fixed before any run (the entry's Design)
QS = (0.0, 0.5, 1.0, 1.5, 2.0)
PS = (0.0, 0.25, 0.5, 1.0, 2.0)
KINDS = ('subtractive', 'divisive')
BREADTHS = ('count', 'summed')
SIMILAR_POWERS = (1, 2, 4)
UNITS = ('z', 'noise')
READINGS = ('trial_means', 'noise_free', 'binary')
N_PERM = 20
N_RANDOM_SPLITS = 200
HOLD_GAIN = 0.005          # A2: a form beyond a uniform inhibition holds
BREADTH_GAIN = 0.005       # A2: the breadth dependence is needed in a unit
SIMILAR_GAIN = 0.003       # A2: the similar-tuning form counts
B_RESPONSES_NOT = 0.004    # B: the flies' form on the model's responses is not lower than the model's APL by this
B_FLIES_GAIN = 0.010       # B: ... while on the flies' responses it is lower than the model's APL by this
B_FORM_GAIN = 0.006        # B: the error is in the form of the model's APL (one wiring draw)
B_FORM_WIRINGS = 6
C_FLOOR_RATIO = 0.5        # C1: usable if its pattern floor is at most half the main target's ...
C_CEILING = 0.90           # ... and its ceiling on ρ is at least this
C_STANDOUT = 0.008         # C3: a variant stands out by this against the breadth-2 fit on the same wiring
C_MEAN_SHARE = 2.0 / 3.0
D_GAIN = 0.003             # D1: a static form does better / worse than the family's APL
D_WIRINGS = 6
D_RAW_SLACK = 0.005
D_UNCHANGED = 0.0482       # the eight wirings' pattern part under the readout against all recordings (E83)
D_RAW_UNCHANGED = 0.0729
D_PROXY_R = 0.8
D_PROXY_DIFF = 0.006

NEAR4 = (('2h', 'IaA'), ('2h', 'pa'), ('IaA', 'pa'), ('eb', 'ep'))


@dataclass(frozen=True)
class Form:
    """One form of inhibition: `kind` 'subtractive' or 'divisive', exponents `q` (on the odor's response over all KCs)
    and `p` (on the KC's breadth), `breadth` 'count' (odors responded to) or 'summed' (summed response)."""
    kind: str = 'subtractive'
    q: float = 0.0
    p: float = 0.0
    breadth: str = 'count'

    @property
    def name(self) -> str:
        raise NotImplementedError


UNIFORM = Form('subtractive', 0.0, 0.0, 'count')


def criterion(unit: str) -> float:
    """The trial-mean response criterion of `unit` ('z': 1; 'noise': the registry's noise_threshold)."""
    raise NotImplementedError


def trialmean_mats(unit: str) -> Dict[str, Dict[tuple, np.ndarray]]:
    """{'TNT3c' / 'TNT': {recording: KCs x 16 trial means (APL16 order)}} in `unit`."""
    raise NotImplementedError


def trial_arrays(unit: str) -> Dict[str, Dict[tuple, np.ndarray]]:
    """{'TNT3c' / 'TNT': {recording: (repeats, KCs, 16) single trials}} in `unit`."""
    raise NotImplementedError


def single_trial_criterion(unit: str) -> float:
    """The single-trial criterion at which the control flies' single-trial response rate is their trial-mean rate."""
    raise NotImplementedError


def strengths(m: np.ndarray, crit: float) -> np.ndarray:
    """A_o: each odor's mean thresholded response over the KCs, relative to its mean over odors."""
    raise NotImplementedError


def breadths(m: np.ndarray, crit: float, measure: str = 'count') -> np.ndarray:
    """B_k: each KC's breadth (odors above the criterion, or summed thresholded response), relative to its mean over
    responding KCs."""
    raise NotImplementedError


def inhibition(m: np.ndarray, form: Form, crit: float) -> np.ndarray:
    """A_o^q * B_k^p (KCs x odors), the form's inhibition before its constant c."""
    raise NotImplementedError


def apply_form(m: np.ndarray, c: float, inh: np.ndarray, kind: str, crit: float) -> np.ndarray:
    """The reduced responses (thresholded convention: values not above the criterion are 0)."""
    raise NotImplementedError


def fit_c(mats: Sequence[np.ndarray], inhs: Sequence[np.ndarray], kind: str, crit: float, rate: float) -> float:
    """The constant c at which the mean response rate over `mats` after the form is `rate` (bisection)."""
    raise NotImplementedError


def similar_inhibition(m: np.ndarray, crit: float, power: float) -> np.ndarray:
    """Each KC's inhibition by the KCs of similar tuning: for odor o, the weighted mean of the other KCs' thresholded
    responses to o, weights the positive part of their correlation with the KC over the other 15 odors, to `power`."""
    raise NotImplementedError


def per_odor_oracle(m: np.ndarray, fracs: np.ndarray, kind: str, crit: float) -> np.ndarray:
    """Each odor's response fraction set to `fracs` (the control flies'), by a per-odor criterion or subtraction."""
    raise NotImplementedError


def mean_distances(mats: Sequence[np.ndarray]) -> pd.Series:
    """The mean over recordings of the cosine distances (canonical pair order) of thresholded-convention matrices."""
    raise NotImplementedError


def scores(pred: pd.Series, data: pd.Series) -> dict:
    """loss, mean, pattern part, ρ, the share of the data's variance explained and the slope, the ester and alcohol
    groups' errors, the four nearest ester pairs and the six near pairs."""
    raise NotImplementedError


def grid() -> List[Form]:
    """A1's forms: q x p x kind x breadth (breadth only matters for p > 0)."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """A1 and A2: forms.csv, forms_pairs.csv, controls.csv, recordings.csv, held_out.csv, floors."""
    raise NotImplementedError


def part_b(out: Path) -> dict:
    """The models' APL-off responses through the forms: model_forms.csv, two_by_two.csv."""
    raise NotImplementedError


def part_c(out: Path) -> dict:
    """The target without APL (apl_free_target.csv), the models against it, the archive scan (archive_scan.csv)."""
    raise NotImplementedError


def part_d1(out: Path) -> dict:
    """The forms as a static inhibition in place of APL on the eight wirings (runs_D.csv, readings_D.csv)."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B', 'C', 'D1', 'D2', 'verdict'), required=True)
    args = ap.parse_args()
    out = H.results_dir(EXPERIMENT)
    {'A': part_a, 'B': part_b, 'C': part_c, 'D1': part_d1}[args.part](out)


if __name__ == '__main__':
    main()
