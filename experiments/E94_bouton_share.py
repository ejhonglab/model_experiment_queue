#!/usr/bin/env python3
"""E94: the PN boutons' share of APL's inhibition in C2's family, from none to four times the fitted dose, as a
subtraction and as a gain, read with the gate for fitted changes; and the family fitted again with that share at the
connectome's.

The entry is item E94 of experiment_priorities.md (proposed by viki 2026-11-14 01:20; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 00:51). STUB: the interface; the bodies follow.

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`, through E91's `WiringModel91`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the
recordings E77's, C / D the second split's (E89's `targets89`). A dose m multiplies the wiring's fitted APL->bouton scale
(`wAPLPN_scale`) and nothing else; the tied APL scale (APL->KC, KC->APL, bouton->APL; E83's TIED_ON) is then bisected
to the wiring's APL-on count in E92's bracket (0.01 to 10). The somatic threshold and the KC-KC scale stay at the fit's.
The APL-off run of every point is the saved fit's.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E91_summed_drive_inhibition as E91

EXPERIMENT = 'E94_bouton_share'
ARM = 'B'

# the settings (the entry's Design; fixed before any run)
FORMS = {'sub': False, 'gain': True}                # bouton_inh_divisive
B1_DOSES = (0.0, 0.1, 0.25, 0.5, 2.0, 4.0)          # the subtraction at exponent 2 (dose 1: the saved fits)
B2_DOSES = (2.0, 4.0)                               # the subtraction at exponent 1 (dose 1: E89's saved point)
B3_DOSES = (0.02, 0.05, 0.15)                       # the gain at exponent 2 (dose 1: E92's B2 point)
ARMS = {'B1': ('sub', 2.0, B1_DOSES), 'B2': ('sub', 1.0, B2_DOSES), 'B3': ('gain', 2.0, B3_DOSES)}
APL_BRACKET = (0.01, 10.0)                          # E92's
BOUND_E34 = 3.0                                     # the APL->bouton bound of E34, in multiples of E3's base
CONNECTOME_SHARE = (0.10, 0.11)                     # Yang's own Claude's tables (hemibrain, male CNS)
ANATOMICAL_SHARE = 0.105
ANATOMICAL_DOSES = (0.1, 0.25, 0.5)
ZERO_SHARE_NOTE = 0.2                               # a chosen setting with more bouton-by-odor entries at zero is named
DIVISOR_HIGH = 10.0                                 # E92's rule

# part C (E91's rules)
CHOICE_MIN_GAIN = 0.0015
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
EQUAL_KINDS = ('lower_share', 'weaker_weighting', 'gain')
HALVES = ('A', 'B', 'C', 'D')

# part D
D_MAXFEV = 150
UNCHANGED_S = {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
D_PASS_GAIN, D_WIRINGS_LOWER, D_RAW_REF, D_RAW_SLACK, D_MEAN_SLACK = 0.003, 6, 0.073, 0.005, 0.01
D2_APLPN_MAX = 12.0
FLIES_MEAN = 0.178


def variant_name(form: str, k: float, dose: float, control: bool = False) -> str:
    """'B_<form>_k<k>_d<dose>' (+ '_ctrl')."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    raise NotImplementedError


class WiringModel94(E91.WiringModel91):
    """E91's WiringModel91 with a dose on the APL->bouton scale (`dose`, 1 the fit's): `run_on(x)` runs APL on with the
    tied scales times x and wAPLPN_scale times the dose."""


def b_reading94(model, wk, dose: float, form: str, ext=None, label: str = '', bracket=APL_BRACKET) -> dict:
    raise NotImplementedError


def claw_loss_shares(W, b_off, b_on, claw_on, kc_of_claw=None) -> dict:
    """The share of the claws' loss of drive under APL that is made at the boutons."""
    raise NotImplementedError


def control_matrix94(d_ref, d_chosen, d_off) -> Tuple[pd.DataFrame, str, float]:
    raise NotImplementedError


def part_0(out: Path, wirings: Sequence[str] = E83.ALL_WIRINGS) -> dict:
    """The identities (dose 1 subtraction = the saved fits; the gain at dose 1 = E92's B2 point on the real wiring;
    dose 0 = a simulator built without the path)."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """A1 the path in the fits, A2 the connectome's share, A3 the bouton stage under other doses, A4 the archive."""
    raise NotImplementedError


def part_b(out: Path, arm: str, wirings: Sequence[str]) -> pd.DataFrame:
    raise NotImplementedError


def score(out: Path) -> dict:
    raise NotImplementedError


def choose94(rd: pd.DataFrame, half: str) -> dict:
    raise NotImplementedError


def equal_rule94(rd: pd.DataFrame) -> List[dict]:
    raise NotImplementedError


def anatomical_dose(runs: pd.DataFrame) -> dict:
    raise NotImplementedError


def part_c(out: Path) -> dict:
    raise NotImplementedError


def part_d(out: Path, pop: str = 'S', force: bool = False, dry_run: bool = False) -> List[str]:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'hostcheck', 'score', 'C', 'D', 'Dscore', 'D2', 'E',
                                                      'steps', 'variants', 'table'))
    ap.add_argument('--arm', choices=tuple(ARMS), default='B1')
    ap.add_argument('--wirings', default='S')
    ap.parse_args()
    raise NotImplementedError('E94: stub')


if __name__ == '__main__':
    main()
