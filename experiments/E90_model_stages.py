#!/usr/bin/env python3
"""E90: the model's stages against the flies' stages. The model's PN, bouton and claw stages read as the PN-dendrite,
bouton and claw recordings were made; whether the ester and the alcohol pairs leave the flies before the claws, at the
step from claws to somata with APL on, or at that step without APL; and the step in wild-type flies as a second reading
of the main target's pattern that uses no fly with APL output blocked (TNT).

The entry is item E90 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from atlas
2026-11-13 11:00). C2's family throughout (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72
part B's eight wirings at its own fit (E83's `saved_fit` / `build_simulator`); a population is the mean over wirings of
the pair distances. The 66 pairs of the panel's 12 megamat odors, in groups (GROUPS): the ester group (15 pairs), the
alcohol group (6), ester with alcohol (24), the other pairs (21), and E78's six near pairs.

Parts (CLI `--part`):
- A (no simulations): A1 the flies' stage table (stages_flies.csv; per pair and group, 90 % intervals from resampling
  each stage's flies, N_BOOT resamples): ORNs (E85's per-fly reader, olfdb dataset 1, the model's 36 glomeruli) and PN
  dendrites (GH146, al_analysis's sent_to_remy/2023-10-29 tables, beside the pebbled ORNs of the same table) across
  glomeruli, per fly then averaged; Remy's PN boutons and KC claws from the per-frame matrices over STAGE_WINDOW, as
  measured (the z-scored `stim` matrices, as data/remy_rdm_timeseries/windows/) and with the noise removed through the
  block averages (each fly's window-mean block-average correlations divided by the geometric mean of the two odors'
  reliabilities, the diagonal); the single-trial correlations of rdm_timeseries_explainers.html beside where they decode;
  the somata (Remy's wild type, Yang's control flies, Yang's flies with APL output blocked (TNT); per fly / recording,
  1 - Pearson and 1 - cosine). A2 the two steps per pair (steps_flies.csv): somata with APL on minus claws (Remy's
  somata; Yang's control flies) and somata with APL output blocked (TNT) minus claws, against Yang's per-pair APL effect
  and the error common to all wirings (E72's common_error_pairs.csv), the step's floor (flies resampled at both stages)
  and the agreement of the claw flies' halves (the three splits of two and two); step_target.json. A3 the antennal lobe
  as recorded (al_transfer.csv): PN dendrites against ORNs on the shared glomeruli (the transfer, each odor's summed
  response, the pair distances, per fly), also under a saturating map from dF/F to rate (the bound for the indicator).
- B0 (simulations): identity: each wiring's run with the rates saved reproduces its saved fit's counts exactly, APL on
  and off; on the real wiring the claw inputs formed as wiring x bouton rate minus APL's inhibition, clipped at 0,
  reproduce a run's `claw_sims` exactly (identity.json).
- B (simulations, resumable per wiring): each wiring's PN rates (pn_sims), bouton rates after APL's inhibition
  (bouton_sims) and APL's activity (inh_sims), APL on and off, read as window means over MODEL_WINDOW minus the
  pre-odor means; the stage matrices (glomeruli / boutons / claws x odors) saved under outputs/E90_model_stages/B/.
  B5: the real wiring at the flies' 2 s pulse at E82's re-held B1 point (its saved counts reproduced first), the stages
  read over PULSE_WINDOW.
- Bscore (no simulations): B1 the PN stage against A3 (entries, per-odor totals, pair distances; pn_stage.csv); B2 the
  bouton stage and B3 the claw stage (the drive alone and with CLAW_SPIKE_SHARES of the KC's spike signal) in the three
  READINGS ('rate', 'session_z', 'noisy': session_z through measurement noise at the flies' block-to-block
  reliability, the noise then removed as in A1); B4 the soma stage (raw counts and the calibrated readout, Pearson and
  cosine); the check of the noise removal (noise_check.csv); stages_model.csv, steps_model.csv.
- C (no simulations): the verdict per group, fixed before any run (VERDICT_* below; verdict.json) and the step target's
  usability (step_target.json).
- D1 / D2 (only as the verdict says): D1 the firing-ceiling family (E69's fit 1, E72 part G's fits on P1, P2, M1) read
  the same way beside C2's; D2 the PN stage of E77's screen's normalization (S1) and exponent (S3) settings beside their
  saved main-target scores (d2_pn_settings.csv).
- table: comparison_table.md and the figures.

Outputs: results/E90_model_stages/ (stages_flies.csv, steps_flies.csv, step_target.json, al_transfer.csv,
stages_model.csv, steps_model.csv, noise_check.csv, verdict.json, comparison_table.md, figures, provenance_<part>.json);
stage matrices and runs under outputs/E90_model_stages/ (archived to data/model_outputs/E90_model_stages/).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E90_model_stages.py --part A
    OMP_NUM_THREADS=12 OMP_THREAD_LIMIT=12 nice -n 10 ./venv/bin/python -u experiments/E90_model_stages.py --part B0
    OMP_NUM_THREADS=12 OMP_THREAD_LIMIT=12 nice -n 10 ./venv/bin/python -u experiments/E90_model_stages.py --part B
    ./venv/bin/python -u experiments/E90_model_stages.py --part Bscore
    ./venv/bin/python -u experiments/E90_model_stages.py --part C
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E90_model_stages'

# ---------------------------------------------------------------------------
# the panel's megamat odors, the pairs and the groups

ODORS = tuple(H.MEGAMAT12)
ESTER = ('2h', 'IaA', 'pa', 'eb', 'ep', 'ms')
ALCOHOL = ('1-8ol', '1-5ol', '1-6ol', '6al')
NEAR = (('eb', 'ep'), ('6al', 't2h'), ('IaA', 'pa'), ('2h', 'pa'), ('1-5ol', '1-6ol'), ('2h', 'IaA'))   # E78's SIX
GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near', 'all66')
GROUP_NAMES = {'ester': 'ester group (15 pairs)', 'alcohol': 'alcohol group (6)', 'ester_alcohol': 'ester with alcohol (24)',
               'other': 'the other pairs (21)', 'near': "E78's near pairs (6)", 'all66': 'all 66 megamat pairs'}

# ---------------------------------------------------------------------------
# part A: the flies

STAGE_WINDOW = (0.5, 3.0)               # Remy's boutons and claws: per-frame matrices averaged over this window (s)
N_BOOT = 1000
CI = (5.0, 95.0)                        # 90 % intervals
BOOT_SEED = 0
REMY_VARIANT = {'bouton': 'pn_{kind}_F_zscore', 'claw': 'kc_{kind}_Fc_zscore'}   # the matrices Remy's figures use
SOMA_SOURCES = ('remy', 'ctrl', 'tnt')  # Remy's wild type; Yang's control flies; Yang's APL output blocked (TNT)
# A3's saturating map from dF/F to rate (the bound for the indicator): r = x / (1 - x / (SAT_HEADROOM * x_max)), x_max the
# stage's largest fly-mean entry
SAT_HEADROOM = 1.2

# ---------------------------------------------------------------------------
# part B: the model

MODEL_WINDOW = (0.0, 0.5)               # the 0.5 s pulse: the stages' window (s)
PRE_WINDOW = (-0.5, 0.0)                # the pre-odor window whose mean is subtracted (s)
PULSE = 2.0                             # B5: the flies' pulse (E82)
PULSE_WINDOW = (0.5, 2.0)               # B5: the stages read over this window (s)
READINGS = ('rate', 'session_z', 'noisy')
CLAW_SPIKE_SHARES = (0.0, 0.25, 0.5)    # B3's versions of what a claw's calcium reports
N_NOISE_REPS = 20                       # replicate noise draws per reading (each of the stage's number of flies)
NOISE_SEED = 0
CLAW_IDENTITY_ODORS = ('eb', '1-8ol')   # B0: the odors of the run with claw_sims saved (memory)

# ---------------------------------------------------------------------------
# parts C and D: the verdict, fixed before any run

VERDICT_DIFF = 0.05                     # the population's group mean this far or more from the flies'
VERDICT_SAME_SIDE = 6                   # of the 8 wirings on the same side
VERDICT_EFFECT_ERR = 0.03               # "not resolved" when nothing differs and the group's error in the effect is above
NOISE_CHECK_TOL = 0.03                  # the noise removal within this of the truth for every group (E88's A0)
STEP_HALVES_MIN = 0.6                   # the step target: the claw flies' halves agree at this or more on average
STEP_FLOOR_RATIO = 0.5                  # and its floor at most this times its SD over the 66 pairs
VERDICTS = ('upstream', 'the step with APL on', 'the step without APL', 'both steps', 'not resolved', 'none')


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# pairs and groups

def megamat_pairs() -> pd.MultiIndex:
    """The 66 megamat pairs in canonical (sorted) order, as E83's `_cs`."""
    raise NotImplementedError


def group_masks(pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """Boolean masks over `pairs` for every name in GROUPS."""
    raise NotImplementedError


def group_means(pairs: pd.Series) -> Dict[str, float]:
    """A pair Series' mean over each group's pairs."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A: readers

def remy_matrices(stage: str, kind: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Remy's per-fly per-frame matrices of `stage` ('bouton' / 'claw'), `kind` 'stim' or 'blockavg' (REMY_VARIANT):
    (flies x frames x 17 x 17 distances 1 - r, frame times, odors)."""
    raise NotImplementedError


def remy_window_means(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> np.ndarray:
    """Per fly, the mean over the window's frames of the `stim` matrices (flies x 17 x 17, 1 - r), as
    decode_scripts/window_csvs.py makes data/remy_rdm_timeseries/windows/ (closed window)."""
    raise NotImplementedError


def remy_noise_removed(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> np.ndarray:
    """Per fly, 1 - (window-mean block-average correlation of two odors) / sqrt(product of the two odors' window-mean
    reliabilities, the diagonal); flies x 17 x 17."""
    raise NotImplementedError


def remy_reliability(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> pd.DataFrame:
    """Per fly and odor, the window-mean block-average reliability (1 - the diagonal); flies x 17."""
    raise NotImplementedError


def remy_single_trial(stage: str) -> Optional[dict]:
    """The single-trial correlations of rdm_timeseries_explainers.html (`D.trial`), if they decode: per fly, the
    between-odor correlations averaged over trial pairs and the same-odor ones; None otherwise."""
    raise NotImplementedError


def al_table(kind: str) -> pd.DataFrame:
    """The megamat trial-mean responses of al_analysis's sent_to_remy/2023-10-29 table `kind` ('GH146' / 'pebbled'):
    rows (fly, odor) with odor the abbreviation, columns glomeruli (NaN where a fly did not have it)."""
    raise NotImplementedError


def shared_glomeruli() -> List[str]:
    """The 28 glomeruli the GH146 recordings and the model's input share (the note of 2026-11-07)."""
    raise NotImplementedError


def per_fly_distances(table: pd.DataFrame, glomeruli: Optional[Sequence[str]] = None,
                      measure: str = 'pearson') -> pd.DataFrame:
    """Per fly the 66 pair distances across the fly's glomeruli (those in `glomeruli` it has); flies x pairs."""
    raise NotImplementedError


def orn_fly_distances() -> pd.DataFrame:
    """The ORN flies' (E85's reader, olfdb dataset 1) per-fly 66 pair distances, 1 - Pearson of the trial means across
    the fly's measured glomerul of the model's 36; flies x pairs."""
    raise NotImplementedError


def remy_soma_fly_distances(measure: str = 'pearson') -> pd.DataFrame:
    """Remy's wild-type megamat KC somata, per fly, the 66 pair distances; flies x pairs."""
    raise NotImplementedError


def yang_soma_distances(genotype: str, measure: str = 'pearson') -> pd.DataFrame:
    """Yang's per-recording 66 megamat pair distances (`H.yang_recording_rdms`); recordings x pairs."""
    raise NotImplementedError


def stage_samples() -> Dict[str, pd.DataFrame]:
    """Every fly stage as units-of-resampling x pairs: {stage key: DataFrame}."""
    raise NotImplementedError


def bootstrap_mean(samples: pd.DataFrame, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> np.ndarray:
    """n_boot x pairs means over resampled rows (flies or recordings, with replacement)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A: tables

def part_a1() -> pd.DataFrame:
    """stages_flies.csv: per stage, reading, measure and pair / group, the mean and its 90 % interval."""
    raise NotImplementedError


def part_a2() -> pd.DataFrame:
    """steps_flies.csv and step_target.json."""
    raise NotImplementedError


def part_a3() -> pd.DataFrame:
    """al_transfer.csv."""
    raise NotImplementedError


def part_a() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# B: the model's stages

def claw_inputs(W: np.ndarray, boutons: np.ndarray, inh: np.ndarray, wAPLKC: np.ndarray) -> np.ndarray:
    """The claws' net drive as olfsysm forms it: wiring x bouton rate at t minus the claw's scaled APL weight times APL's
    activity at t - 1, clipped at 0 (claws x time; `boutons` boutons x time, `inh` time)."""
    raise NotImplementedError


class StageRun:
    """One wiring at its fit (C2's family), its simulator built with olfsysm's variables captured (E43's
    `captured_olfsysm_vars`), APL's activity saved; `run(cond)` returns the condition's outputs and its stage matrices
    (PN glomeruli, boutons, claws x odors: window means minus pre-odor means)."""

    def __init__(self, label: str):
        raise NotImplementedError

    def run(self, cond: str) -> dict:
        raise NotImplementedError


def part_b0() -> dict:
    raise NotImplementedError


def part_b(wirings: Optional[Sequence[str]] = None) -> None:
    raise NotImplementedError


def part_b5() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# B: readings

def session_z(values: pd.DataFrame) -> pd.DataFrame:
    """Each unit divided by its SD over a simulated session of the panel's odors (`pyolfsysm.readout.session_zscore`,
    linear amplitude, no noise)."""
    raise NotImplementedError


def noisy_reading(values: pd.DataFrame, reliability: pd.Series, n_flies: int, rng: np.random.Generator) -> dict:
    """`values` (units x odors, session z) through measurement noise at the flies' per-odor block-to-block reliability:
    three trials per odor, per pseudo-fly; returns the as-measured and the noise-removed pair distances (means over the
    pseudo-flies)."""
    raise NotImplementedError


def part_bscore() -> pd.DataFrame:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# C: the verdict

def differs(model_pop: float, model_wirings: Sequence[float], flies: float, lo: float, hi: float) -> dict:
    """The rule for one stage or step, group and reading: |population - flies| >= VERDICT_DIFF, the population outside
    the flies' 90 % interval [lo, hi], and at least VERDICT_SAME_SIDE of the wirings on the population's side."""
    raise NotImplementedError


def group_verdict(table: pd.DataFrame, effect_err: float) -> dict:
    """The verdict of one group from its rows (stage or step, reading, the rule's result)."""
    raise NotImplementedError


def part_c() -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# D

def part_d1() -> pd.DataFrame:
    raise NotImplementedError


def part_d2() -> pd.DataFrame:
    raise NotImplementedError


def comparison_table() -> None:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'B0', 'B', 'B5', 'Bscore', 'C', 'D1', 'D2', 'table'])
    ap.add_argument('--wirings', nargs='*', default=None)
    a = ap.parse_args()
    {'A': part_a, 'B0': part_b0, 'B': lambda: part_b(a.wirings), 'B5': part_b5, 'Bscore': part_bscore, 'C': part_c,
     'D1': part_d1, 'D2': part_d2, 'table': comparison_table}[a.part]()


if __name__ == '__main__':
    main()
