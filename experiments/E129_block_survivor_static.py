#!/usr/bin/env python3
"""E129: the inhibition that outlasts the block reduced to a static level: in the blocked condition only, the family's
inhibition among KCs replaced by its time average, a constant per KC, a constant per KC type, or one constant for all
KCs, against a pooled neuron and a slow matrix.

The entry is item E129 of experiment_priorities.md (proposed by hal 2026-11-22 13:29 at Tom's request of 12:57, moved
into the list under the idle-host rule; planned by a Fable subagent, claude-fable-5-2, spawned from hal 13:05). "The
family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S
is real, P1, P2, M1 and T is P3, P4, M2, M3. "KC-KC at zero" is E110's part A (KC-KC 0 in both conditions, both response
rates held at fixed scales).

The protocol (design item 1; nothing fitted). On each wiring the APL-on run is the family's saved run, untouched, and
the common somatic threshold (fixed_thr) and every other parameter stay the fit's. Only the APL-off condition is
simulated again (`OffRunner`, one of E83's simulators of the wiring at its fit, APL->KC and APL->bouton 0): the KC-KC
matrix is replaced by the arm's inhibition, with one scale bisected (E69's `bisect_factor`, from factor 1) to hold the
family's APL-off response count (`target_off`, the response rate 0.150). Any change of the per-pair APL effect then comes
from the distances with APL off.

The arms (tags S1_<arm>_<wiring>, S2_<arm>_<wiring>, ...):
- T (step 1a), the time course removed: a static inhibition (olfsysm's `ext_inh`, KCs x odors, subtracted from each
  KC's input while the odor is on, where KC-KC enters and in the same units, `E95_static_apl.set_ext_inh_ko`) equal to
  each KC's and odor's time-mean of the family's own KC-KC input with APL off, and the KC-KC scale 0.
- G (1a), one level: the same constant for every KC and odor.
- P (1a), a pooled neuron: E122's uniform neuron U (every entry total / N^2, the diagonal included, the hemibrain
  matrix's total weight) in the matrix's place, filtered with 50 ms as the matrix is.
- S (1a), slow: the hemibrain matrix filtered with 1 s (olfsysm's `kc_kc_tau`, set on the APL-off simulator only).
- K (1b, gated), a constant per KC: its mean over the 16 odors of T's inhibition. Y (1b), a constant per KC type: the
  mean over the type's KCs and the odors. KP (1b, gated again), K's constants permuted among the KCs of each type.
- Z (1c, gated), G's run with the level removed: the APL-off run at the fit's threshold with no inhibition among KCs and
  no static inhibition, not bisected (the acute reading).
- C50, C0 (step 2, gated; APL-on runs): the hemibrain matrix at 0.5 and 0 of the fit's KC-KC scale with APL on, the tied
  APL scales (APL->KC, KC->APL, bouton->APL) bisected to the family's APL-on count, the threshold the fit's; each paired
  with the APL-off held run of G and of Y where that arm at least keeps the gain (kinds 'C50|G', 'C0|G', 'C50|Y',
  'C0|Y').
- MB (step 3): MB-C1's counted synapses in the matrix's place with APL off. Refused without --tom-approved-mbc1 (Tom's
  word on Pratyush's hemibrain MB-C1 tables, which the loader does not read); not implemented until then.
- References, read and not run: the family; KC-KC at zero (E110's part A, 'E110_A'); the step-0 preview 'Gpre'.

Readings the entry leaves open, fixed here before any run:
1. T's time-mean. olfsysm filters each KC's spikes (one vesicle per spike, nves 1, checked: `mp.kc.ves_p` 0) into an
   activity whose integral over the trial is the KC's spike count, so the family's KC-KC input to KC i for odor o has the
   trial integral s x sum_j HB[i, j] n[j, o] (s the fit's KC-KC scale, HB the hemibrain matrix as the simulator holds it,
   n the family's saved APL-off spike counts, counted over the whole simulated trial). T's inhibition is that integral
   divided by the pulse length (0.5 s, the simulator's stim window, checked): the family's total KC-KC input delivered
   evenly while the odor is on (`time_mean_input`). It is the time-mean over the pulse wherever spikes and their 50 ms
   tails fall inside it, and conserves the trial's total otherwise.
2. The static arms' base levels at factor 1 all have T's grand mean: G T's mean over KCs and odors, K each KC's mean
   over the odors, Y the mean over the KCs of each type and the odors, KP K's values permuted within each type
   (numpy's default_rng(KP_SEED), the permutation drawn per type in TYPES order) (`static_base`). Their factors are
   therefore comparable, and each arm's held level is its factor times its base.
3. The one scale. Static arms: the factor multiplies the base matrix (ext_inh_scale stays 1), with the KC-KC scale 0.
   P and S: the factor multiplies the APL-off simulator's KC-KC scale (s x f). Every arm: E69's `bisect_factor` from 1,
   bracket SCALE_BRACKET (1/50 to 50), first step SCALE_STEP (1.25), to `target_off` (two or more spikes). Step 2's APL
   factor: E83's APL_BRACKET (0.1 to 10) and APL_STEP (1.25).
4. A corner: the scale at a bound of its bracket, or the response rate not held, i.e. the count more than RATE_TOL
   (0.2 %) of `target_off` away from it (about eight responses, a rate within 0.0003). An arm with corners on more than
   two wirings gets no class.
5. P's matrix is E122's U; S's 1 s is set as `mp.kc.kc_kc_tau` on the APL-off simulator after it is built (the
   simulator does not reset it per call; checked before every run), not through E10_joint_fit.py.
6. The classes (E123's and E126's, unchanged; calibrated pattern part, S's and T's populations apart, against all
   recordings and each of the four halves A to D):
   - "equal" (E123.classify, one-sided): not higher than the family by more than WITHIN (0.003) on both populations and
     on 6 of 8 wirings, on every target, with every run in the four bands;
   - "keeps the gain": at or below the half-way point between the family and KC-KC at zero on both populations, every
     target (0.0554 on S and 0.0577 on T against all recordings, checked);
   - "neither" otherwise; "corner (no class)" by reading 4;
   - two arms "differ" only by E126's rule (one lower on 7 of 8 wirings against all recordings and on both populations
     against all five targets), else "not told apart".
   "At least keeps the gain" is the class 'equal' or 'keeps the gain'. "T and G fall in different classes" (step 1b's
   gate) and "K and Y fall in different classes" (KP's) compare the class strings, a corner included.
7. Step 1c's gate: G at least keeps the gain. Step 2's gate: G or Y at least keeps the gain; its 16 APL-on runs are the
   same whichever static arm they are paired with, so they are run once and paired with each arm that opens it.
8. The step-0 preview 'Gpre': the family's saved APL-on runs paired with E110 part A's APL-off runs (KC-KC 0, the
   threshold bisected to the same count): a uniform level while the odor is on is the threshold raised while the odor
   is on, up to Vm's lag at the step. Scored and classed as an arm (its classes are a preview, not G's); score reports
   G minus Gpre per population once G is run.
9. The split (step 0): T's base matrix (KCs x odors, Hz) per wiring, over all KCs and over the KCs that respond with APL
   off in the family: the sum of squares about the grand mean in shares per KC type, per KC within its type, per odor
   and the rest (adding to 1; 'share_kc' the first two together, the static per-KC part), the per type x odor part of
   the rest beside, the CV over the odors at the median KC, and the mean level per type in Hz and against fixed_thr
   (`split_kckc`).
10. Realism of the static arms: the held level (mean, per type, largest) in Hz and against each KC's threshold as the
    simulator sets it (fixed_thr + 2 x spontaneous input, read from olfsysm after the run) and against fixed_thr; the
    share of KC x odor entries whose level is above the KC's threshold ('level_above_thr_share'), flagged as hard to
    defend when positive ('level_above_thr').

Parts (CLI `--part`):
- step0 (no simulation): the preview and the split; step0_preview.csv, step0_preview.json, step0_split.csv,
  provenance_step0.json.
- identity (simulations, `--wirings`): on every wiring, (i) the hemibrain matrix at scale 1 through `OffRunner` gives the
  saved APL-off run exactly, (ii) with a static inhibition of zeros set as well, also exactly; on 'real', (iii) T's
  inhibition at factor 1 with KC-KC 0 changes the counts (the static inhibition is read). identity_<host>_<lane>.json;
  step 1a refuses to run until every wiring has passed (`identity_passed`).
- step1a, step1b, step1c, step2 (simulations, resumable, per lane `--wirings`; a point is claimed by mkdir of
  results/E129_block_survivor_static/claims/<tag>/ and skipped when done or claimed, `--reclaim` taking over a claim left
  by a process that died): rows in runs_<step>_<host>_<lane>.csv, provenance (`H.model_provenance` per point) in
  provenance_<step>_<host>_<lane>.json; runs saved locally (outputs/E129_block_survivor_static/<step>/) and in the
  model-output archive (data/model_outputs/E129_block_survivor_static/<step>/), static inhibitions beside. Steps 1b, 1c
  and 2 refuse unless verdict.json (from `score`) opens their gate.
- step3: refused without --tom-approved-mbc1.
- score (no simulation): readings.csv, runs.csv, verdict.json (classes, gates, comparisons, outcomes), comparison_table.md,
  provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- plan: the points a lane would run.

The outcomes (flags for the write-up, which reads them): (a) G equal; (b) Y or K equal and G not; (c) T equal and K not;
(d) T not at least keeping the gain while P or S does; (e) S against P by the "differ" rule; (f) corners, or no arm
told apart from another. No correlation chooses anything, no pair is weighted, and no threshold moves after a run is
read.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it, and the model loader's hemibrain KC->KC counts and KC types (step 3's MB-C1 tables wait on Tom's word); no olfsysm
change; nothing is published.

Usage (from the repository root; 16 threads, under nice; one lane per population):
    ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part identity
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step1a --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step1a --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E129_block_survivor_static.py --part score
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E129_block_survivor_static'
SCRIPT = 'experiments/E129_block_survivor_static.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
TYPES = ('g', "a'b'", 'ab')
STATIC_ARMS = ('T', 'G', 'K', 'Y', 'KP')
STEPS = {'step1a': ('G', 'T', 'P', 'S'), 'step1b': ('K', 'Y', 'KP'), 'step1c': ('Z',), 'step2': ('C50', 'C0')}
PULSE_S = 0.5
SLOW_TAU_S = 1.0
KP_SEED = 0
SCALE_BRACKET, SCALE_STEP = (1 / 50, 50.0), 1.25
RATE_TOL = 0.002
WITHIN = 0.003
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7


def out_dir() -> Path:
    """<shared>/results/E129_block_survivor_static/ (created)."""
    return H.results_dir(EXPERIMENT)


def time_mean_input(counts: np.ndarray, M: np.ndarray, scale: float, pulse_s: float = PULSE_S) -> np.ndarray:
    """Reading 1: scale x M @ counts / pulse_s (KCs x odors, Hz), M post x pre, counts KCs x odors."""
    raise NotImplementedError


def static_base(arm: str, Tm: np.ndarray, types: np.ndarray, seed: int = KP_SEED) -> np.ndarray:
    """Reading 2: the base matrix (KCs x odors) of a static arm from T's `Tm`: 'T', 'G', 'K', 'Y', 'KP'."""
    raise NotImplementedError


def split_kckc(A: np.ndarray, types: np.ndarray) -> dict:
    """Reading 9: the shares per type, per KC within type, per odor, the rest, the type x odor part of the rest."""
    raise NotImplementedError


class OffRunner:
    """The APL-off condition of one wiring at its fit for one arm (the module docstring): `run(f)` -> (responses, spike
    counts) at the arm's scale factor f."""

    def __init__(self, w: str, arm: str):
        raise NotImplementedError

    def run(self, f: float = 1.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        raise NotImplementedError


def part_step0() -> dict:
    """The preview and the split (no simulation)."""
    raise NotImplementedError


def part_identity(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> dict:
    """The identities (simulations)."""
    raise NotImplementedError


def identity_passed() -> bool:
    raise NotImplementedError


def part_runs(step: str, wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False
              ) -> pd.DataFrame:
    """Steps 1a, 1b, 1c and 2 (resumable, claimed per point)."""
    raise NotImplementedError


def part_step3(tom_approved_mbc1: bool = False) -> None:
    """Step 3, refused without Tom's word on the MB-C1 tables."""
    if not tom_approved_mbc1:
        raise SystemExit('E129 step 3 waits on Tom\'s word on Pratyush\'s hemibrain MB-C1 tables (--tom-approved-mbc1)')
    raise NotImplementedError


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    raise NotImplementedError


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    raise NotImplementedError


def gates(classes: Mapping[str, dict]) -> dict:
    raise NotImplementedError


def outcomes(classes: Mapping[str, dict], comps: Mapping[str, dict]) -> dict:
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
