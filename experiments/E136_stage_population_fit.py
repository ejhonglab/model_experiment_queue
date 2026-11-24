#!/usr/bin/env python3
"""E136, the best family fitted to all megamat stages at once: one loss over the PN, bouton, claw and soma populations,
five structural arms, a control fit on odor-permuted flies, and reads on held-out wirings, held-out flies and three
panels the fit never sees.

The entry is item E136 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04; its code claimed by disco 2026-11-24 13:29, after
atlas handed it on). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own
fit (E83's `members`, `saved_fit`, `build_simulator`); S is real, P1, P2, M1 and T is P3, P4, M2, M3. The measures, the
readings and the stage table are those of the models as populations at every recorded stage (E134): the model read as
each stage is recorded (atlas's kcq/stage_readings.py), the views, one-fly units and population errors (hal's
kcq/population_fit.py), the flies' populations and the lock (viki's data_prep/stage_population_data.py,
experiments/E134_lock.json), the family built and run with every stage saved (atlas's
experiments/E134_population_baseline.py: `build_family`, `run_stages`, `held_targets`). This module calls them and
copies none of them; a change those files need is asked of their owners.

The guard: nothing is fitted before E134's measure gate reads "usable" (results/E134_population_baseline/
measure_gate.json, written by viki from the gate planner's reading) and its stage table is in (score_classes.json,
score_stage_table.csv); `require_gate` refuses otherwise. A1 also waits for E135's class (results/
E135_antennal_lobe_dendrites/, `require_e135`); the other arms do not.

The loss (the entry's definition, fixed before any fit): the megamat stage error over the stages that E134 found
readable (`readable_stages`, from score_classes.json), on the S wirings' mean views (hal's `wirings_as_flies`, row
'common'), against a fixed half of each stage's flies (`fly_halves`, half 'fit'), each stage read in its default
reading (`R.DEFAULT_READINGS`) with its curve set anew at every evaluation on pooled values only and the half's own
noise and units; the PN dendrites on four views (`error_with_matrix`), the other stages on three (E134's
`megamat_stage_error`). Beside it, the operating point (`OperatingPoint`): the response rate with APL on at the
megamat flies' value (the lock's 'remy_megamat_criterion_rate', 0.0915, by two spikes), the rate with APL off at
OFF_ON_RATIO times that (2.2, the flies' ratio on the APL16 panel), as a realism band used as an operating point.

The arms (`ARMS`), one fit each, every one with at most MAX_FREE free parameters, none per odor or glomerulus:
- A0, the family as it is, refitted under the new loss;
- A1, with the antennal-lobe stage that E135 hands on (its fitted global parameters, with the kept glomerulus terms
  where they carry, all fixed, entering by E135's route, a table of PN rates); if E135 hands on none, the
  normalization's two parameters are free;
- A2, APL's inhibition of the claws as a gain in place of the subtraction (E92's form, fit_mb_model's
  `apl_inh_divisive`);
- A3, the per-claw threshold free (E16's form, fit_mb_model's `claw_thr`);
- A4, without the inhibition among KCs (E110's form, the KC-KC scale at 0);
- the control: A0 against flies whose odor labels are permuted at every stage (`permuted_cache`).

Parts (CLI `--part`):
- identity: at the family's own parameters (each S wiring at its held megamat point of E134's step 3), the fit's loss
  reproduces E134's stage table (score_stage_table.csv, group 'megamat_held', the 'common' rows) to rounding when read
  against all flies, and every run reproduces E134's saved step-3 counts exactly (`identity`).
- fit: one arm's fit (`--arm`), Nelder-Mead in log factors of each S wiring's own parameters (`ArmModel.params_for`),
  every evaluation's loss, parts, rates and parameters in fit_history_<arm>.csv, the stages of every evaluation saved
  under outputs/ (`save_eval`), resumable from its history (`fit_arm`).
- read: an arm's fit at fixed parameters (`read_arm`): on the T wirings against the held-out half of the flies; on the
  eight wirings with the validation2 and the mixture inputs; on the APL16 input with APL on and off (about 130 runs
  over the arms; lanes and claims as E134's step 3, `planned_reads`).
- score: the classes (`classify`), the comparison table, realism (`realism_rows`), provenance.
- commands: the lanes' command lines.

The classes (the entry's, fixed beforehand), on T and the held-out half of the flies, against the unrefitted family
there:
- "improves": the megamat stage error is lower by at least IMPROVE_SE times the standard error over wirings, and no
  stage's population error is higher by more than that stage's own spread;
- "holds": it improves, and the KC population error is not higher on validation2, on the mixtures or on the APL16
  panel by more than the spread there;
- "megamat only": it improves and a held-out panel is worse;
- "one stage at another's cost": the megamat stage error is lower and a stage is higher by more than its spread;
- an arm "adds" if it improves on A0 by the same margin;
- no class for a fit outside a band, at a bound, or no better than the permuted control.

Realism (reported, not optimized): for every fit the four bands on the APL16 input (kcq/apl_effect.py's), the megamat
response rate beside Remy's flies' 0.0915, the largest spike counts, and every fitted scale beside the family's; named
in advance: the APL->bouton scale against its bound of three times the baseline, the bouton share of APL's inhibition
against the connectome's tenth, any per-claw threshold above 70 Hz, and in A1 any glomerulus term carried over from
E135. A fit that leaves a band gets a "Realism flag:" post and no class.

Concentrations: the fit sees megamat only, at 10^-3 on both sides. The validation2 read is one panel at matched
concentrations, and the mixtures are used as measured. The APL16 read involves paa, euc and -aPine, whose input is the
ORN data within 1 log10 unit of the flies' concentrations.

Holds (CLAUDE.md): only Pratyush's hemibrain (the family's wirings); no olfsysm change for A0, A2, A3 and A4 (their
forms are options at olfsysm b1e89d7); A1's stage enters by E135's route; no input search; only Yang's processed
products; Remy's processed data without the bouton fly 2022-10-21 fly02 (viki's lock excludes it). Nothing is
published.

Outputs:
- outputs/E136_stage_population_fit/ (local, large): fit_<arm>/ (every evaluation's stages, E134's step-3 layout),
  read/ (the reads' runs).
- <shared>/results/E136_stage_population_fit/ (compact): fit_history_<arm>.csv, fit_<arm>.json (the best parameters
  as absolute values and factors, the loss and its parts, the stopping reason), identity.json, reads_<host>_<lane>.csv,
  claims/, scores.csv, classes.json, realism.csv, comparison_table.md, variants.csv, provenance_<part>.json.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E136_stage_population_fit.py \
        --part fit --arm A0

STATUS: interface stub (2026-11-24, disco). Bodies raise NotImplementedError until the open definitions
(`OPEN_DEFINITIONS`) are ruled and the runner is written.
"""
from __future__ import annotations

import argparse
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E136_stage_population_fit'
SCRIPT = 'experiments/E136_stage_population_fit.py'

# the family's wirings (E83's, as E134 takes them; `check_wirings` checks they are E134's)
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS

# the stages a loss can read (E134's MEGAMAT_ERROR_STAGES); which enter is E134's verdict (`readable_stages`)
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
# the fly halves: 'fit' (the loss's) and 'held' (the reads')
HALVES: Tuple[str, ...] = ('fit', 'held')
# the operating point: APL off at this multiple of the held APL-on rate (the entry's 2.2, the flies' APL16 ratio)
OFF_ON_RATIO = 2.2
MAX_FREE = 8
# the classes' margin: the improvement needed, in standard errors over the wirings
IMPROVE_SE = 2.0
# the read panels and conditions
READ_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')

# the arms: name -> (description, the free parameters, the fixed model kwargs beside the family's). The free
# parameters are factors on each wiring's own parameters (`ArmModel.params_for`); 'thr' the somatic threshold
# (fixed_thr), 'apl' the APL->KC scale with the KC->APL and bouton->APL scales tied to it at each wiring's own ratios
# (E72's TIED), 'aplpn' the APL->bouton scale, 'kckc' the KC-KC scale, 'claw_thr' the per-claw threshold (Hz,
# absolute), 'pn_inhadd' / 'pn_inhsc' the PN normalization's two parameters (A1 when E135 hands on none).
# PROVISIONAL until the open definitions are ruled (OPEN_DEFINITIONS 'free_parameters').
Arm = namedtuple('Arm', 'name description free kws')
ARMS: Dict[str, Arm] = {
    'A0': Arm('A0', 'the family as it is, refitted under the new loss', ('thr', 'apl', 'aplpn', 'kckc'), {}),
    'A1': Arm('A1', "with the antennal-lobe stage that E135 hands on, fixed", ('thr', 'apl', 'aplpn', 'kckc'),
              {'pn_rates': 'E135'}),
    'A2': Arm('A2', "APL's inhibition of the claws as a gain (E92's form)", ('thr', 'apl', 'aplpn', 'kckc'),
              {'apl_inh_divisive': True}),
    'A3': Arm('A3', "the per-claw threshold free (E16's form)", ('thr', 'apl', 'aplpn', 'kckc', 'claw_thr'), {}),
    'A4': Arm('A4', "without the inhibition among KCs (E110's form)", ('thr', 'apl', 'aplpn'), {'kckc_fixed': 0.0}),
}
CONTROL = 'A0_perm'

# what the entry leaves open, for a planner to rule before the fits (the stub's defaults are marked PROVISIONAL)
OPEN_DEFINITIONS: Dict[str, str] = {
    'halves': "which flies make each stage's fixed fit half and held-out half (seed, sizes for the 4- and 5-fly "
              "stages), and whether Remy's same-fly fly may sit in both the claw and the soma halves",
    'operating_point': "how the APL-on rate (0.0915) and the APL-off rate (2.2 x) are held in the fit: re-held by "
                       "bisection at every evaluation (about 25 simulations per wiring) or penalty terms in the loss "
                       "(their scales and weights against the stage error in one-fly units)",
    'free_parameters': "each arm's free parameters (at most 8): whether the tied KC->APL and bouton->APL ratios are "
                       "freed, and A2's and A3's added parameters",
    'start': "the start: each wiring's held megamat point (E134's megamat_held) or its fit's own parameters",
    'optimizer': "the optimizer and its budget: Nelder-Mead in log factors (E72 part D), its initial simplex, maxfev "
                 "(8 to 10 h per fit at about 3 min per evaluation)",
    'claw_reading': "the claw stage under A2 (the divisive gain) and A3 (the per-claw threshold): E90's claw_inputs "
                    "reads the subtractive net drive before any per-claw transfer",
    'permuted_control': "the control's permutation: one odor permutation for all stages and flies, or per fly; its "
                        "seed",
    'bounds': "bounds on the factors (the APL->bouton scale's three times the baseline, the per-claw threshold's "
              "range) and what 'at a bound' means for the classes",
}

ReadPoint = namedtuple('ReadPoint', 'arm panel wiring half')


# ---------------------------------------------------------------------------
# paths and guards

def out_dir() -> Path:
    """<shared>/results/E136_stage_population_fit/."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    """outputs/E136_stage_population_fit/<part>/ (created)."""
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def check_wirings() -> None:
    """S_WIRINGS and T_WIRINGS are E134's (and so E83's)."""
    from experiments import E134_population_baseline as B
    if (S_WIRINGS, T_WIRINGS) != (tuple(B.S_WIRINGS), tuple(B.T_WIRINGS)):
        raise RuntimeError(f'the wirings differ from E134\'s: {S_WIRINGS} {T_WIRINGS}')


def require_gate() -> dict:
    """E134's measure gate (results/E134_population_baseline/measure_gate.json) and its classes (score_classes.json);
    RuntimeError unless the gate's verdict is 'usable' and the stage table is in."""
    raise NotImplementedError


def require_e135() -> dict:
    """E135's class and what it hands on to A1 (its fitted antennal-lobe stage: global parameters, kept glomerulus
    terms, the table of PN rates by its route); RuntimeError while E135 has no class."""
    raise NotImplementedError


def readable_stages(classes: Optional[Mapping] = None) -> Tuple[str, ...]:
    """The stages of LOSS_STAGES that E134 found readable (score_classes.json's stage classes, 'not readable' left
    out), in stage order."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the flies: halves and the permuted control

def fly_halves(seed: int = 0) -> Dict[str, Dict[str, List[str]]]:
    """Per stage of LOSS_STAGES, the flies (viki's fly names) of the fit half and of the held-out half (HALVES), fixed
    before any fit. PROVISIONAL (OPEN_DEFINITIONS 'halves')."""
    raise NotImplementedError


def half_cache(half: str, stages: Sequence[str] = LOSS_STAGES, panel: str = 'megamat') -> 'R.FlyCache':
    """A FlyCache holding only the flies of one half (atlas's `subset_flies` of viki's populations) for the megamat
    stages, so that every reading of the loss (curve, noise, criterion, one-fly units) is the half's own."""
    raise NotImplementedError


def permuted_cache(half: str, seed: int, stages: Sequence[str] = LOSS_STAGES) -> 'R.FlyCache':
    """`half_cache` with the flies' odor labels permuted at every stage (the control fit). PROVISIONAL
    (OPEN_DEFINITIONS 'permuted_control')."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the model: an arm on several wirings

class ArmModel:
    """One arm of the family on a set of wirings and one panel: each wiring built once (E134's `build_family`, with
    the arm's kwargs), run at shared factors on its own parameters with every stage saved (E134's `run_stages`).

    Attributes: arm, panel, wirings, base (per wiring, the starting absolute parameters), ratios (per wiring, the tied
    KC->APL and bouton->APL ratios), built (per wiring, E134's `Built`)."""

    def __init__(self, arm: str, wirings: Sequence[str] = S_WIRINGS, panel: str = 'megamat',
                 start: str = 'held'):
        raise NotImplementedError

    def params_for(self, wiring: str, z: Mapping[str, float]) -> Tuple[Dict[str, float], Dict[str, float]]:
        """The APL-on and APL-off absolute parameters of a wiring at log factors `z` (keys the arm's free parameters;
        0 is the start)."""
        raise NotImplementedError

    def run(self, z: Mapping[str, float]) -> Dict[str, 'R.SavedRun']:
        """Every wiring at `z`, both conditions, every stage saved, as in-memory runs (`in_memory_run`)."""
        raise NotImplementedError


def in_memory_run(tag: str, out: Mapping[str, dict], built, info: Mapping[str, object]) -> 'R.SavedRun':
    """E134's `run_stages` output as atlas's `SavedRun` without writing it: the layout `load_step3_run` gives after
    `save_step3_run` (the settled pre-odor means in place of E90's window, odors abbreviated), checked by
    `R.check_meta`; the test checks a saved and reloaded run reads the same."""
    raise NotImplementedError


def save_eval(arm: str, i: int, runs: Mapping[str, 'R.SavedRun']) -> List[Path]:
    """An evaluation's runs under outputs/.../fit_<arm>/ (E134's step-3 layout, tag '<arm>_e<i>_<wiring>')."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the loss

class OperatingPoint:
    """The operating point's targets and its term in the loss: the APL-on rate at the lock's megamat rate (0.0915, by
    two spikes) and the APL-off rate at OFF_ON_RATIO times that, on the wirings' mean rates. PROVISIONAL
    (OPEN_DEFINITIONS 'operating_point')."""

    def __init__(self, lock: Optional[Mapping] = None):
        raise NotImplementedError

    def term(self, rate_on: float, rate_off: float) -> float:
        raise NotImplementedError


class StageLoss:
    """The fit's loss: the megamat stage error of the wirings' mean views over `stages`, against the flies of
    `flies` (a half's FlyCache), with the operating point's term; the parts, every stage's population error and views,
    each wiring's errors and the rates beside it in the returned dict."""

    def __init__(self, stages: Sequence[str], flies: 'R.FlyCache', operating_point: OperatingPoint,
                 lock: Optional[Mapping] = None):
        raise NotImplementedError

    def __call__(self, runs: Mapping[str, 'R.SavedRun']) -> Dict[str, float]:
        raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

def identity(arm: str = 'A0') -> dict:
    """At z = 0 (each S wiring at its held megamat point of E134's step 3): every run reproduces E134's saved step-3
    counts exactly, and the loss's stage errors against all flies reproduce E134's score_stage_table.csv ('megamat_held',
    the 'common' rows) to rounding; written to identity.json; RuntimeError otherwise."""
    raise NotImplementedError


def fit_arm(arm: str, maxfev: Optional[int] = None, control: bool = False, resume: bool = True) -> dict:
    """One arm's fit on the S wirings against the fit half (the control: A0 against the permuted half), resumable
    from fit_history_<arm>.csv; writes fit_<arm>.json."""
    raise NotImplementedError


def planned_reads(arms: Optional[Sequence[str]] = None) -> List[ReadPoint]:
    """The reads at fixed parameters: each arm's fit (and the unrefitted family) on the T wirings against the held-out
    half, on the eight wirings with the validation2 and mixture inputs, and on the APL16 input."""
    raise NotImplementedError


def read_arm(pt: ReadPoint, lane: str, reclaim: bool = False) -> Optional[dict]:
    """One read, claimed as E134's step-3 points are (mkdir of claims/<tag>/), run, saved and scored; its row in
    reads_<host>_<lane>.csv."""
    raise NotImplementedError


def classify(scores: pd.DataFrame) -> dict:
    """The entry's classes per arm (module docstring), the control's comparison and the "adds" test against A0."""
    raise NotImplementedError


def realism_rows(arm: str) -> pd.DataFrame:
    """Realism of an arm's fit (module docstring): the APL16 bands, the megamat rate beside 0.0915, the largest counts,
    every fitted scale beside the family's and the named-in-advance checks."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('identity', 'fit', 'read', 'score', 'commands'))
    ap.add_argument('--arm', choices=tuple(ARMS) + (CONTROL,))
    ap.add_argument('--maxfev', type=int)
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reclaim', action='store_true')
    a = ap.parse_args(argv)
    raise NotImplementedError(f'--part {a.part}: E136 is an interface stub')


if __name__ == '__main__':
    main()
