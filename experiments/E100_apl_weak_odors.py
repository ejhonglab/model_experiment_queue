#!/usr/bin/env python3
"""E100: the flies' APL at the weak odors and at the panel's concentrations. Yang's further APL recordings (set W, the
January 2025 weak-odor set; set S, the 2022 series with pa, paa and ms from 10^-9 to 10^-1; set H, the 2h and ma series of
2021-12 to 2022-01; set P, the one fly with the 16 panel odors) read for how strongly and where in the calyx APL responds,
against the best family's APL at those stimuli and against the connectome-predicted maps. A diagnostic: it moves no loss.

The entry is item E100 of experiment_priorities.md (proposed by hal 2026-11-15 13:44; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-15 12:58, at Tom's request in hal's terminal). C2's family (the breadth-2
APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit, nothing re-tuned. The
recordings: notes/yang_apl_recordings_inventory.md in the shared directory.

Two dependencies on other hosts' work (not duplicated here):
- E98 (the strength of APL's response in the flies against the model's APL; atlas; experiments/E98_apl_response_strength.py).
  E100 reuses its reader of the 2021 maps, its ratios, its indicator family (`indicator_ratio`, `model_ratio_set`,
  `classify`, `comparison_verdict`) and its fixed-parameter runner (part B1, `ExtendedRun`). Part A2 is read from E98's
  saved B1 runs (atlas was asked to add 4h-ol at 10^-6 and farn at 10^-2 to B1), by `read_e98_b1`; E100 runs no model of
  its own for A2. Until B1's outputs exist, A2 raises NotImplementedError naming E98's B1.
- viki's copy of Yang's processed APL recordings to data/yang_apl_recordings/ in the shared directory (recordings.csv,
  raw_index.csv, README.md). Parts B to E read only that copy (never Yang's folders) and only after Yang has answered
  hal's questions (the entry's Requires). Until the copy's README is in place, their readers follow the interface below
  (`STIM_TABLE_COLUMNS`, `RECORDING_COLUMNS`), tested on made-up data only.

Definitions (the entry's, fixed before any reading; E98's hold where not changed here):
- A recording's calyx response to a presentation: the calyx-mean peak dF/F inside the outline (Yang's peak maps,
  `*_dff_peak5frames_ave.tif`, the same product as the 2021 set's). Set W's two changes to E98's definitions: the
  solvent's response is the mean of the recording's own solvent presentations (one per block, last in the block;
  `recording_solvent`), and part B gives every response raw, minus the recording's solvent and minus its block's solvent
  (`block_solvent`); where a response can be at or below the solvent's, the mean over flies is arithmetic, and the table
  says so (`fly_ratio_w`).
- Concentration matching (Tom's rule): a stimulus without ORN input at its own concentration takes the input of the same
  odor within MATCH_LOG_UNITS (1 log10 unit, inclusive); each session's estimate (SESSIONS, the diagnostics of the lab's
  three ORN sessions) is kept separate (`input_estimates`). paa at 10^-3 takes the panel's input at 10^-4, as everywhere.
- W1 (the weak odors), six ratios over paa at 10^-3 (W1_RATIOS). Per ratio and session estimate, E98's comparison verdict
  over the eight wirings (6 of 8); a ratio is inside if it is inside under any session's estimate, and flatter or steeper
  only if so under all (`ratio_over_estimates`); W1 is flatter, inside or steeper if at least W1_MIN_RATIOS (4) of the six
  are, otherwise not decided (`w1_verdict`). Ratios over farn at 10^-2 (exact input) are given beside and decide nothing.
- W2 (set S), ms over pa and paa over pa, all at 10^-3: decided only with W2_MIN_FLIES (4) or more flies, otherwise
  described (`w2_verdict`).
- W3 (the maps against the connectome), on the pairs A4 keeps (power at least A4_POWER_MIN, 0.8): a pair passes at a share
  s if the prediction exceeds the label-permuted null's 95th percentile (W3_NULL_PCT) in more than s of the recordings in at
  least W3_MIN_CALYCES (5) of the 8 calyces. W3 passes if at least W3_MIN_PAIRS (3) pairs pass at W3_PASS_SHARE (1/2); it
  fails if no pair passes at W3_FAIL_SHARE (1/4); otherwise not decided (`w3_verdict`).

Parts (CLI `--part`):
- 0, the identities; a failure stops the experiment and is reported (identity.json). Those that need no new recording,
  run now (`part_0 --no-recordings`): the PN stage on the extended input reproduces E79's cached panel PN rates
  (connectome_pn_rates.csv) exactly; the 1-log matcher reproduces data/yang_weak_odor_kc/orn_input_availability.csv's
  statuses; the calyx-response reader reproduces E79's calyx means on one 2021 recording (E79_maps_across_baselines.csv).
  After viki's copy: the reader reproduces Yang's calyx-response tables from her maps and outlines; each recording's pixel
  size and channel list come from its Experiment.xml. E98's own part 0 checks that B1 reproduces E90's saved runs.
- A (before any new recording is read): A1, the inventory as a table (from the copy's stimulus and settings files; written
  by viki's copy, read here); A2, the family at the stimuli of sets W, S and H that have ORN input, read from E98's B1
  (`read_e98_b1`): APL's window and pre-odor means, its evoked activity and the KCs' response rates with APL on and acutely
  off, on the eight wirings, under each session's estimate (a2_family.csv); A3, the connectome-predicted maps (E79's
  bouton-driven prediction with C2's PN stage at the weak stimuli, eight calyces) for A3_PAIRS, with the anatomy-alone and
  label-permuted controls (a3_pn_rates.csv, a3_inputs.csv, a3_maps.csv, a3_agreement.csv, a3_anatomy.csv); A4, the power
  of W3 from the 2021 maps' reliability (a4_power.csv, a4_power_by_calyx.csv; a pair under A4_POWER_MIN leaves W3's list,
  a4_w3_pairs.json).
- B (after Yang's answer): per recording and stimulus the calyx response raw, minus the recording's solvent and minus the
  block's solvent; block-to-block agreement; drift; sensitivity (the planes apart; without recordings whose solvent
  response exceeds SOLVENT_LARGE_SHARE of their largest response); E98's saturating indicator beside the family for set
  S; the table per weak odor (the flies' APL, the family's APL, the flies' KCs, the family's KCs with APL on and off).
- C (after Yang's answer): E79's analyses 1 to 3 on sets W, S and H (repeat reliability; the common map's share; the
  remainder's reliability and length; whether the solvent's map is the common map; the dorsal-ventral centroid at
  matched magnitude; the common map's share at the weak odors against the panel's concentrations).
- D (after Yang's answer): E79's analysis 4 on A3's pairs (A4's list), the orientation fixed as far as Yang has given it
  (posterior view, dorsal up, "L" the fly's right; both medial edges until she answers: D_FLIPS), slabs within
  D_SLAB_WITHIN_UM of the plane's nominal depth; set S's pa minus ms at 10^-3 beside.
- E (only if Yang has set P processed or agrees to its processing): APL's calyx response to each of the 16 panel odors in
  one fly against the family's profile and E98's proportional alternatives. A lead, no verdict.
- score: the verdict (W1, W2, W3; verdict.json), fixed above.

Deviations from the entry (each also in the note):
- A2 runs nothing: it reads E98's B1 (the entry allows either; the task given to this script says so).
- A4's "true wiring-based map" is another calyx's predicted map (the seven other calyces; the model's own hemibrain boutons
  and hemibrain R are not truths for each other), at a slab within D_SLAB_WITHIN_UM of a nominal plane depth drawn from
  A4_PLANE_DEPTHS_UM and at a blur drawn from E79's three, plus independent noise per grid cell at the reliability of the
  2021 set's measured difference maps (the median split-half r of E79's pa minus ms and paa minus ms, read from
  connectome_measured.csv: 0.59 to 0.69), the scored map being the mean of A4_N_PRESENTATIONS (2) presentations as in 2021.
  Power is the share of simulated recordings, pooled over scored calyces, truths and draws, whose best score over slabs,
  blurs and D_FLIPS exceeds the null's 95th percentile of the permutations' maxima over the same search; the lower end of
  the reliability range decides (a4_power.csv gives both, and the share per scored calyx).
- A3's "contrast against the label-permuted null": per calyx, pair and session estimate, the SD over grid cells of the
  predicted difference map (each odor's map over its mean inside the outline, as E79) and the agreement between calyces of
  the predicted difference maps (mean Pearson r over calyx pairs, the left calyces mirrored), each against the same in
  N_PERM permutations of the glomerulus labels; the anatomy-alone control is the r of each odor's predicted map with the
  unweighted PN-to-APL synapse map (which predicts no difference between odors).

CLI (from the repository root):
    ./venv/bin/python -u -m experiments.E100_apl_weak_odors --part 0 --no-recordings
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E100_apl_weak_odors --part A3
    ./venv/bin/python -u -m experiments.E100_apl_weak_odors --part A4
    ./venv/bin/python -u -m experiments.E100_apl_weak_odors --part A2     # after E98's B1
    ./venv/bin/python -u -m experiments.E100_apl_weak_odors --part B      # after Yang's answer and viki's copy
    --dry-run lists what a part would read and write, and exits.

Outputs: results/E100_apl_weak_odors/ in the shared directory (the tables above, identity.json, verdict.json,
provenance_<part>.json); the PN-stage run of A3 under outputs/E100_apl_weak_odors/, archived to
data/model_outputs/E100_apl_weak_odors/.
"""
from __future__ import annotations

import argparse
import json
import socket
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E100_apl_weak_odors'
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's / E98's

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any reading)

Stim = Tuple[str, float]                 # (odor abbreviation, log10 concentration)
MATCH_LOG_UNITS = 1.0                    # Tom's rule, inclusive
SESSIONS = ('kiwi-control', 'megamat', 'validation2')   # the lab's three ORN sessions, each with its diagnostics
DIAG_PANEL = 'glomeruli_diagnostics'
SOLVENT = 'solvent'                      # the stimulus table's label of the solvent (paraffin oil: 'pfo' read as it)
SOLVENT_LABELS = ('solvent', 'pfo', 'paraffin oil')
SOLVENT_LARGE_SHARE = 1.0 / 3.0          # B's sensitivity: recordings whose solvent response exceeds this share of their
                                         # largest response are left out
MIN_RESPONSE = 0.02                      # E98's: the denominator of a ratio must exceed this (dF/F)
N_BOOT = 10000
SEED = 0
CI = (5.0, 95.0)

# set W's stimuli (the inventory), set S's and the denominators
SET_W_STIMULI: Tuple[Stim, ...] = (('2h', -7.0), ('2h', -6.0), ('ma', -7.0), ('4h-ol', -6.0), ('farn', -3.0),
                                   ('farn', -2.0), ('2-but', -6.0), ('ms', -5.5), ('ms', -4.5), ('paa', -6.0),
                                   ('paa', -5.0), ('paa', -3.0))
SET_S_STIMULI: Tuple[Stim, ...] = tuple([('pa', float(c)) for c in range(-9, -1)] + [('paa', float(c)) for c in range(-8, 0)]
                                        + [('ms', float(c)) for c in range(-8, 0)])
W1_DENOMINATOR: Stim = ('paa', -3.0)
W1_BESIDE_DENOMINATOR: Stim = ('farn', -2.0)
W1_RATIOS: Tuple[Stim, ...] = (('2h', -6.0), ('ma', -7.0), ('4h-ol', -6.0), ('farn', -2.0), ('2-but', -6.0), ('paa', -5.0))
W1_MIN_RATIOS = 4                        # of the six
W2_RATIOS: Tuple[Tuple[str, Stim, Stim], ...] = (('ms/pa', ('ms', -3.0), ('pa', -3.0)),
                                                 ('paa/pa', ('paa', -3.0), ('pa', -3.0)))
W2_MIN_FLIES = 4
VERDICT_MIN_WIRINGS = 6                  # E98's, of the 8 wirings

# A3 / A4 / W3: pairs of weak stimuli that drive different single glomeruli (each odor's map over its mean, a minus b)
A3_PAIRS: Tuple[Tuple[Stim, Stim], ...] = ((('2h', -6.0), ('ma', -7.0)), (('2h', -6.0), ('farn', -2.0)),
                                           (('ma', -7.0), ('farn', -2.0)), (('4h-ol', -6.0), ('2h', -6.0)),
                                           (('2-but', -6.0), ('ma', -7.0)))
N_PERM = 1000                            # E79's label permutations
A4_POWER_MIN = 0.8
A4_N_DRAWS = 20                          # noise draws per (scored calyx, truth calyx, reliability)
A4_N_PRESENTATIONS = 2                   # presentations averaged in a scored map (the 2021 set's)
A4_PLANE_DEPTHS_UM = (12.0, 15.0, 18.0, 21.0, 24.0)    # set W's planes (the inventory)
W3_NULL_PCT = 95.0
W3_PASS_SHARE = 0.5
W3_FAIL_SHARE = 0.25
W3_MIN_CALYCES = 5                       # of the 8
W3_MIN_PAIRS = 3
D_FLIPS = ('none', 'horizontal')         # posterior view, dorsal up; both medial edges until Yang answers
D_SLAB_WITHIN_UM = 6.0

# the interface assumed for viki's copy until its README is in place (one row per presentation / per recording)
STIM_TABLE_COLUMNS = ('recording', 'frame', 'block', 'odor', 'log10_conc')
RECORDING_COLUMNS = ('recording', 'set', 'date', 'fly', 'side', 'plane_um', 'pulse_s', 'peak_maps', 'roi', 'stimuli',
                     'experiment_xml')
COPY_DIR_NAME = 'yang_apl_recordings'


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def copy_dir() -> Path:
    """viki's copy of Yang's processed APL recordings (data/yang_apl_recordings/ in the shared directory)."""
    return H.shared_dir() / 'data' / COPY_DIR_NAME


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(path.read_text()) if path.exists() else None


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(),
                                         time=time.strftime('%Y-%m-%d %H:%M'), **extra))
    _write_json(path, d)


def stim_label(s: Stim) -> str:
    """'<odor> @ <log10 conc>' with the concentration as the lab's tables write it ('2h @ -6', 'ms @ -5.5')."""
    raise NotImplementedError


def parse_stim(label: str) -> Stim:
    """The inverse of `stim_label` ('2h @ -6' -> ('2h', -6.0))."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# concentrations and the ORN input

def concentration_match(c: float, available: Sequence[float], max_units: float = MATCH_LOG_UNITS) -> Optional[float]:
    """The available concentration nearest to `c` if within `max_units` log10 units (inclusive), else None; of two equally
    near, the higher (the stronger input, as E49)."""
    raise NotImplementedError


def session_inputs(sessions: Sequence[str] = SESSIONS) -> Dict[str, pd.DataFrame]:
    """Each session's ORN spike deltas (olfsysm's pyolfsysm/data/internal/<session>_mean_est_spike_deltas.csv; glomerulus x
    (panel, 'odor @ conc'))."""
    raise NotImplementedError


def input_estimates(stim: Stim, inputs: Optional[Dict[str, pd.DataFrame]] = None,
                    panel: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Every estimate of a stimulus's ORN input, one row per session: the session's diagnostics column of the same odor
    within MATCH_LOG_UNITS (`concentration_match`), and the panel's own input (the APL16 input, `panel`, as session
    'panel') when the odor is a panel odor within the same rule; columns session, panel_name, column, input_conc,
    log10_distance, status ('exact' or 'substitute'). Empty when no input is within the rule."""
    raise NotImplementedError


def extended_input(stims: Sequence[Stim], panel: Optional[pd.DataFrame] = None
                   ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The APL16 input with one added column per (stimulus, estimate) (panel 'E100', column '<stim label> [<session>]'),
    restricted to the APL16 input's 36 glomeruli (glomeruli a session lacks filled with 0), and the table of the added
    columns (stimulus, session, source column, status, the glomeruli filled with 0, the input total and the top glomerulus
    with its delta)."""
    raise NotImplementedError


def availability_identity() -> dict:
    """Part 0: `input_estimates` reproduces the statuses of data/yang_weak_odor_kc/orn_input_availability.csv (exact /
    substitute / none for every listed monomolecular stimulus)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the flies: the stimulus table, the shuffled blocks, the solvent, the ratios (parts B to E; viki's copy)

def read_stim_table(path: Optional[Path] = None) -> pd.DataFrame:
    """The presentations of every recording of viki's copy (STIM_TABLE_COLUMNS), solvent labels mapped to SOLVENT."""
    raise NotImplementedError


def shuffled_blocks(frames: pd.DataFrame) -> pd.DataFrame:
    """One recording's presentations (frame order) with the shuffled-block design read: `block` (taken from the table, or,
    where absent, a new block after each solvent presentation, the solvent being last in its block), `position` (1-based
    in the block), `is_solvent`, `rep` (the stimulus's presentation index, 0-based) and the flags `solvent_last` (the
    block ends with its one solvent), `complete` (the block holds every stimulus of the recording once) and `duplicate` (a
    stimulus presented twice in one block)."""
    raise NotImplementedError


def recording_solvent(resp: pd.DataFrame) -> float:
    """The recording's solvent response: the mean of its solvent presentations' responses (one per block). `resp` has
    is_solvent and response."""
    raise NotImplementedError


def block_solvent(resp: pd.DataFrame) -> pd.Series:
    """Each block's solvent response (its own solvent presentation; NaN for a block without one), by block."""
    raise NotImplementedError


def solvent_subtracted(resp: pd.DataFrame) -> pd.DataFrame:
    """Per presentation (odors only) the response raw, minus the recording's solvent ('minus_recording') and minus its
    block's solvent ('minus_block'), with the stimulus means over presentations of each."""
    raise NotImplementedError


def large_solvent(resp: pd.DataFrame, share: float = SOLVENT_LARGE_SHARE) -> bool:
    """Whether the recording's solvent response exceeds `share` of its largest stimulus response (mean over presentations
    of the raw responses)."""
    raise NotImplementedError


def fly_ratio_w(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The flies' ratio from per-recording numerator and denominator responses (columns fly, num, den; solvent already
    subtracted): recordings whose denominator is at most MIN_RESPONSE left out; the geometric mean over flies (each fly's
    recordings on the log scale) when every ratio is positive, else the arithmetic mean (each fly's recordings averaged
    arithmetically), with 'mean' saying which; the 90 % interval from `n_boot` bootstrap draws of flies; the counts."""
    raise NotImplementedError


def calyx_responses(stack: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Per presentation (frame of a peak-map stack), the mean dF/F inside the calyx outline."""
    raise NotImplementedError


def experiment_xml(path: Path) -> dict:
    """A ThorImage Experiment.xml: the pixel size (um), the enabled channels (by name) and the frame size."""
    raise NotImplementedError


def reader_identity_2021(recording: Optional[str] = None) -> dict:
    """Part 0: on one 2021 recording (default the first of recordings.csv), the calyx means of `calyx_responses` with the
    order-matched, first-in-block and all-paraffin-oil baselines reproduce E79's E79_maps_across_baselines.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the verdict rules

def ratio_over_estimates(per_estimate: Dict[str, str]) -> str:
    """One ratio's class from its comparison verdict under each session's input estimate: 'inside' if inside under any;
    'flatter' or 'steeper' only if that under all; otherwise 'not decided'."""
    raise NotImplementedError


def w1_verdict(ratio_classes: Dict[str, str], min_ratios: int = W1_MIN_RATIOS) -> str:
    """W1 from the six ratios' classes: 'flatter', 'inside' or 'steeper' if at least `min_ratios` are; otherwise
    'not decided'."""
    raise NotImplementedError


def w2_verdict(comparison: str, n_flies: int, min_flies: int = W2_MIN_FLIES) -> str:
    """A W2 ratio's verdict: its comparison verdict with `min_flies` or more flies, otherwise 'described'."""
    raise NotImplementedError


def w3_pair_calyces(exceed: pd.DataFrame, share: float) -> int:
    """The number of calyces in which the prediction exceeds the null's 95th percentile in more than `share` of the
    recordings (`exceed`: one row per calyx and recording, column exceeds, bool)."""
    raise NotImplementedError


def w3_verdict(exceed: pd.DataFrame, pairs: Optional[Sequence[str]] = None) -> dict:
    """W3 from one row per pair, calyx and recording (columns pair, calyx, recording, exceeds) on the pairs A4 kept:
    'passes' if at least W3_MIN_PAIRS pairs pass at W3_PASS_SHARE in at least W3_MIN_CALYCES calyces, 'fails' if no pair
    does at W3_FAIL_SHARE, else 'not decided'; with the per-pair counts."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A2: the family at the weak stimuli, read from E98's B1

def read_e98_b1(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """E98's saved B1 runs: one row per wiring, condition ('on', 'acute_off') and B1 stimulus, with APL's window mean
    ('win'), pre-odor mean ('pre') and the KCs' response rate. Read from outputs/E98_apl_response_strength/B/ or its
    archive (data/model_outputs/E98_apl_response_strength/B/), files `<wiring>_<cond>_apl.csv` (index stimulus, columns
    win, pre; E90's `_save_stage` layout, which E98's stub names as its source) and `<wiring>_<cond>_responses.parquet`
    (KC x stimulus); E98's own loader is used instead once it has one. Raises NotImplementedError naming E98's B1 while
    they are not there."""
    raise NotImplementedError


def part_a2() -> dict:
    """A2: a2_family.csv (per wiring, stimulus of sets W, S and H with ORN input, session estimate and condition: APL's
    window, pre-odor and evoked activity and the KCs' response rate) and the family's W1 and W2 ratios through E98's
    indicator family (a2_model_ratios.csv)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A3: the connectome-predicted maps

def pn_rates_extended(stims: Sequence[Stim] = (), rebuild: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """C2's PN stage (E79's `model_pn_rates` recipe: C2's simulator, APL off, the 0-0.5 s mean of pn_sims minus the
    -0.5-0 s mean) on the extended input (`extended_input`): glomerulus x column evoked rates (the 16 panel odors by
    abbreviation, the added columns by their labels) and the added columns' table; cached in a3_pn_rates.csv and
    a3_inputs.csv."""
    raise NotImplementedError


def pn_identity(rates: pd.DataFrame) -> dict:
    """Part 0: the 16 panel odors' evoked PN rates of the extended run equal E79's connectome_pn_rates.csv (largest
    absolute difference)."""
    raise NotImplementedError


def pair_maps(Wt: np.ndarray, perms: np.ndarray, G: np.ndarray, inside: np.ndarray) -> np.ndarray:
    """The predicted difference maps (permutations x 144) of a pair from the column x 2 weights `Wt` (the pair's two
    stimuli's PN rates per glomerulus), each permutation's column order and the columns' grids G: each stimulus's map over
    its mean inside the outline, a minus b (E79's `predicted_comparisons` for a pair)."""
    raise NotImplementedError


def part_a3(n_perm: int = N_PERM) -> dict:
    """A3 (module docstring): a3_maps.csv (per calyx, pair, session estimate, slab and blur: the SD of the predicted
    difference map and its percentile in the label-permuted null), a3_agreement.csv (per pair and estimate, the mean r
    between calyces' whole-calyx difference maps, against the null), a3_anatomy.csv (each stimulus's map against the
    anatomy-alone map) and a3_grids.npz (the predicted grids for A4 and part D)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A4: the power of W3

def measured_reliability() -> Dict[str, float]:
    """The 2021 set's measured difference-map reliability: the median over recordings of E79's split-half r of pa minus ms
    and paa minus ms (connectome_measured.csv)."""
    raise NotImplementedError


def noise_sd(signal: np.ndarray, reliability: float, n_presentations: int = A4_N_PRESENTATIONS) -> float:
    """The SD of independent per-cell noise added to a true map `signal` (its inside cells) so that two single
    presentations would agree at Pearson r `reliability`, for the mean of `n_presentations` presentations:
    sd(signal) sqrt((1 - r) / r / n)."""
    raise NotImplementedError


def detect(pred: Dict[Tuple[float, float], np.ndarray], inside: Dict[Tuple[float, float], np.ndarray], meas: np.ndarray,
           meas_inside: np.ndarray, flips: Sequence[str] = D_FLIPS) -> Tuple[float, float, bool]:
    """One recording against one calyx: `pred` maps (slab depth, blur) to the predicted grids (row 0 the true labels, the
    rest permutations; permutations x 144), `inside` to their inside flags. The best r over slabs, blurs and flips, the
    95th percentile of the permutations' maxima over the same search, and whether the first exceeds the second."""
    raise NotImplementedError


def power(pair: str, grids: dict, reliabilities: Sequence[float], n_draws: int = A4_N_DRAWS,
          seed: int = SEED) -> pd.DataFrame:
    """A4 for one pair: simulated recordings (a truth calyx's predicted map plus noise at each reliability) scored against
    every other calyx's prediction (`detect`); one row per scored calyx, truth calyx, reliability and draw."""
    raise NotImplementedError


def part_a4() -> dict:
    """A4 (module docstring): a4_power.csv, a4_power_by_calyx.csv, a4_w3_pairs.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# parts 0, B to E, score

def part_0(recordings: bool = False) -> dict:
    """Part 0 (module docstring): identity.json. With `recordings` False, only the identities that need no new
    recording."""
    raise NotImplementedError


def part_b() -> dict:
    """Part B (module docstring), from viki's copy, after Yang's answer."""
    raise NotImplementedError


def part_c() -> dict:
    """Part C (module docstring), from viki's copy, after Yang's answer."""
    raise NotImplementedError


def part_d() -> dict:
    """Part D (module docstring), from viki's copy, after Yang's answer."""
    raise NotImplementedError


def part_e() -> dict:
    """Part E (module docstring), only if set P is processed."""
    raise NotImplementedError


def score() -> dict:
    """The verdict (W1, W2, W3) from parts A2, A4, B and D: verdict.json."""
    raise NotImplementedError


PARTS = ('0', 'A2', 'A3', 'A4', 'B', 'C', 'D', 'E', 'score')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--no-recordings', action='store_true', help="part 0: the identities that need no new recording")
    ap.add_argument('--n-perm', type=int, default=N_PERM, help='A3: label permutations')
    ap.add_argument('--dry-run', action='store_true', help='say what the part reads and writes, and exit')
    a = ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
