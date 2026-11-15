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
  hal's questions (the entry's Requires). Their readers follow the copy's layout (its README.md and recordings.csv;
  `read_recordings`, `recording_paths`, `yaml_stimuli`, `read_recording`), tested on
  made-up data shaped like it and on the 12 recordings Yang released (data/yang_apl_calyx_maps/); `read_recording` refuses
  every other recording until YANG_RELEASED is set after her answer.

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
  run now (`part_0 --no-recordings`): the PN stage on the extended input reproduces, for the 16 panel odors, E90's saved
  real-wiring PN window means exactly and, read with the step times pyolfsysm used before olfsysm 17cef71, E79's cached
  panel PN rates (connectome_pn_rates.csv) exactly (`pn_identity`); the 1-log matcher reproduces data/yang_weak_odor_kc/orn_input_availability.csv's
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
  the reliability range decides (a4_power.csv gives both, and the share per scored calyx). Beside, as a sensitivity that
  decides nothing: the reliability of a map with half the response at the same noise (`scaled_reliability`), since set
  W's calyx responses are 0.13 to 0.64 of paa at 10^-3 (the planner's probes). Truth and prediction take the
  same session's input estimate, and a pair's power is its lowest over the three sessions' estimates (the entry fixes no
  rule over estimates for W3; part D will score each estimate apart, and W3 passes or fails only if it does so under
  every estimate, otherwise it is not decided, as W1's outside classes).
  As a control that decides nothing, each calyx is also scored against
  its own prediction plus noise (`power_own_calyx`): the power if every fly's calyx were wired like the scored one.
- A3 reads the PN stage with pyolfsysm's current step times (olfsysm 17cef71 on), as E90 does; E79's cached panel rates
  were read with the earlier ones, which put the 0-0.5 s window a fraction of a step later, so its evoked rates differ
  from A3's by up to 0.021 Hz (0.3 %); part 0 checks both readings.
- Each session's estimate is a column of that session's table at the concentration nearest the stimulus's; where a
  session holds that concentration in two panels (megamat's ms at 10^-3 in the diagnostics and in its own panel), both
  are kept as separate estimates. The PN stage runs on the 36 glomeruli of the panel's input, as everywhere; the weak
  odors' most-driven glomeruli (VM5d, DM4, VM2, DC3, VM5v, VM7d) are among them.
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
import zlib
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
A4_RESPONSE_SCALE = 0.5                # A4's sensitivity (described, decides nothing): the reliability of a map whose
                                       # response is this fraction of the 2021 maps' at the same noise
W3_NULL_PCT = 95.0
W3_PASS_SHARE = 0.5
W3_FAIL_SHARE = 0.25
W3_MIN_CALYCES = 5                       # of the 8
W3_MIN_PAIRS = 3
D_FLIPS = ('none', 'horizontal')         # posterior view, dorsal up; both medial edges until Yang answers
D_SLAB_WITHIN_UM = 6.0

# viki's copy (data/yang_apl_recordings/, README.md there): recordings.csv (one row per recording; set, stem, plane_depth_um,
# pulse_s, nblock, nstim, peak_map_copied, shared_path, ...) and per recording `<stem>.yaml` (odor_pin, odor_stim_seq, nblock,
# nstim), `<stem>[_C1]_<reg>_dff_peak5frames_ave.tif` (one peak dF/F map per presentation in odor_stim_seq order), `<stem>*.roi`,
# `<stem>*_roi_peak_dff.csv` (Yang's calyx peak per presentation) and `<stem>_Experiment.xml`
STIM_TABLE_COLUMNS = ('recording', 'frame', 'block', 'odor', 'log10_conc')
E100_SETS = ('W', 'S', 'H')
# Yang has not released the recordings beyond the 12 of data/yang_apl_calyx_maps/ (hal's questions, the entry's Requires):
# `read_recording` refuses any other recording until this is set True, after her answer
YANG_RELEASED = False
SOLVENT_PIN = 37
# Yang's odor names (her yaml odor_pin / odor_pin_full) -> abbreviations
ODOR_ABBREV = {'paraffin oil': SOLVENT, 'pfo': SOLVENT, 'phenylacetaldehyde': 'paa', 'pentyl acetate': 'pa',
               'methyl salicylate': 'ms', '2-heptanone': '2h', 'methyl acetate': 'ma', 'farnesol': 'farn',
               '4-heptanol': '4h-ol', '2-butanone': '2-but', 'geranyl acetate': 'ga'}
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


def _e79():
    """E79's two modules (E79_connectome_maps, E79_apl_calyx_constraints), with their shared-directory paths repointed to
    `H.shared_dir()` (and Yang's Claude's connectome tables beside it, Tom/claude_APL_modeling_YZ4TC) on hosts where the
    tensor mount is not /mnt/tensor, the path they were written with (hal mounts it at /local/tensor)."""
    from data_prep import apl_calyx_maps as ACM
    from experiments import E79_apl_calyx_constraints as DS
    from experiments import E79_connectome_maps as CM
    shared = H.shared_dir()
    if not Path(ACM.SHARED).exists() and Path(ACM.SHARED) != shared:
        old = Path(ACM.SHARED)

        def re(x):
            x = Path(x)
            return shared / x.relative_to(old) if x.is_relative_to(old) else x
        for mod, names in ((ACM, ('DATA', 'OUT')), (DS, ('DATA', 'OUT')),
                           (CM, ('OUT', 'C2_ON_COUNTS', 'COMMON_ERROR', 'ORN_INPUT'))):
            for n in names:
                setattr(mod, n, re(getattr(mod, n)))
            mod.SHARED = shared
        yang_apl = shared.parent / Path(CM.YANG_APL).relative_to(old.parent)
        if yang_apl.exists():
            CM.YANG_APL = yang_apl
    return CM, DS


def stim_label(s: Stim) -> str:
    """'<odor> @ <log10 conc>' with the concentration as the lab's tables write it ('2h @ -6', 'ms @ -5.5')."""
    return f'{s[0]} @ {float(s[1]):g}'


def parse_stim(label: str) -> Stim:
    """The inverse of `stim_label` ('2h @ -6' -> ('2h', -6.0))."""
    odor, conc = label.split(' @ ')
    return odor.strip(), float(conc)


def pair_label(a: Stim, b: Stim) -> str:
    return f'{stim_label(a)} - {stim_label(b)}'


# ---------------------------------------------------------------------------
# concentrations and the ORN input

def concentration_match(c: float, available: Sequence[float], max_units: float = MATCH_LOG_UNITS) -> Optional[float]:
    """The available concentration nearest to `c` if within `max_units` log10 units (inclusive), else None; of two equally
    near, the higher (the stronger input, as E49)."""
    av = sorted({float(x) for x in available}, reverse=True)
    if not av:
        return None
    best = min(av, key=lambda x: abs(x - float(c)))      # min keeps the first, the higher, of a tie
    return best if abs(best - float(c)) <= max_units + 1e-9 else None


def orn_delta_dir() -> Path:
    """olfsysm's pyolfsysm/data/internal (the installed pyolfsysm's, found without importing it)."""
    import importlib.util
    return next(Path(p) for p in importlib.util.find_spec('pyolfsysm').submodule_search_locations
                if Path(p).is_dir()) / 'data' / 'internal'


def session_inputs(sessions: Sequence[str] = SESSIONS) -> Dict[str, pd.DataFrame]:
    """Each session's ORN spike deltas (olfsysm's pyolfsysm/data/internal/<session>_mean_est_spike_deltas.csv; glomerulus x
    (panel, 'odor @ conc'))."""
    out = {}
    for s in sessions:
        d = pd.read_csv(orn_delta_dir() / f'{s}_mean_est_spike_deltas.csv', header=[0, 1], index_col=0)
        d.index.name = 'glomerulus'
        out[s] = d.astype(float)
    return out


def _single_columns(df: pd.DataFrame, odor: str) -> List[Tuple[str, str, float]]:
    """(panel, column, log10 conc) of the single-odor columns of `odor` in one session's table."""
    out = []
    for panel, col in df.columns:
        if ' @ ' not in col or '+' in col.split(' @ ')[0]:
            continue
        try:
            o, c = parse_stim(col)
        except ValueError:
            continue
        if o == odor:
            out.append((panel, col, c))
    return out


def input_estimates(stim: Stim, inputs: Optional[Dict[str, pd.DataFrame]] = None) -> pd.DataFrame:
    """Every estimate of a stimulus's ORN input: per session (each of the lab's three ORN sessions, SESSIONS), the
    session's columns of the same odor at the concentration nearest to the stimulus's (`concentration_match`, within
    MATCH_LOG_UNITS, inclusive), one row per panel holding it (the diagnostics and the session's own panel, e.g. megamat's
    ms at 10^-3 in both; each a separate measurement, kept separate); columns stimulus, session, panel_name, column,
    input_conc, log10_distance, status ('exact' or 'substitute'), estimate (the column's label in the extended input).
    Empty when no session has input within the rule. The APL16 panel's input is among them (megamat's megamat panel and
    validation2's validation2 panel)."""
    inputs = session_inputs() if inputs is None else inputs
    rows = []
    for s, df in inputs.items():
        cols = _single_columns(df, stim[0])
        c = concentration_match(stim[1], [x[2] for x in cols])
        if c is None:
            continue
        for panel, col, cc in cols:
            if abs(cc - c) > 1e-9:
                continue
            dist = abs(cc - float(stim[1]))
            rows.append(dict(stimulus=stim_label(stim), session=s, panel_name=panel, column=col, input_conc=cc,
                             log10_distance=round(dist, 3), status='exact' if dist < 1e-9 else 'substitute',
                             estimate=estimate_label((stim[0], cc), s, panel)))
    return pd.DataFrame(rows, columns=['stimulus', 'session', 'panel_name', 'column', 'input_conc', 'log10_distance',
                                       'status', 'estimate'])


def estimate_label(source: Stim, session: str, panel: str) -> str:
    """The extended input's column for one measured input (`source`: the odor and the concentration of the session's
    column): '<odor> [<session>/<panel>] @ <conc>' (one ' @ ', as fit_mb_model's odor-name parsing wants), the panel
    shortened ('diag' for the diagnostics). Stimuli that the 1-log rule gives the same column share it."""
    p = 'diag' if panel == DIAG_PANEL else panel
    return f'{source[0]} [{session}/{p}] @ {float(source[1]):g}'


def extended_input(stims: Sequence[Stim], panel: Optional[pd.DataFrame] = None,
                   inputs: Optional[Dict[str, pd.DataFrame]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The APL16 input with one added column per (stimulus, estimate) (panel 'E100', column `estimate_label`), restricted
    to the APL16 input's 36 glomeruli (glomeruli a session lacks filled with 0), and the table of the added columns
    (`input_estimates`'s columns, the glomeruli filled with 0 and dropped, the input total over the 36 glomeruli, the top
    glomerulus with its delta, and the input vector itself as one column per glomerulus, `orn_<glomerulus>`)."""
    panel = H.apl16_orn_deltas() if panel is None else panel
    inputs = session_inputs() if inputs is None else inputs
    gl = list(panel.index)
    cols, rows = {}, []
    for st in stims:
        for r in input_estimates(st, inputs).to_dict('records'):
            if r['estimate'] in cols:
                rows.append(dict(next(x for x in rows if x['estimate'] == r['estimate']), **r))
                continue
            src = inputs[r['session']]
            v = src[(r['panel_name'], r['column'])].reindex(gl)
            filled = sorted(g for g in gl if g not in src.index or pd.isna(src[(r['panel_name'], r['column'])].get(g)))
            v = v.fillna(0.0).astype(float)
            cols[r['estimate']] = v
            top = v.idxmax()
            rows.append(dict(r, glomeruli_filled_0=';'.join(filled),
                             glomeruli_dropped=';'.join(sorted(set(src.index) - set(gl))), input_total=float(v.sum()),
                             top_glomerulus=top, top_delta=float(v[top]),
                             **{f'orn_{g}': float(v[g]) for g in gl}))
    extra = pd.DataFrame(cols, index=panel.index)
    extra.columns = pd.MultiIndex.from_tuples([('E100', c) for c in extra.columns], names=panel.columns.names)
    orn = pd.concat([panel, extra], axis='columns')
    orn.index.name = panel.index.name
    return orn, pd.DataFrame(rows)


def availability_identity() -> dict:
    """Part 0: `input_estimates` reproduces data/yang_weak_odor_kc/orn_input_availability.csv's spike-delta statuses and
    closest concentrations (exact / substitute / none) for every listed monomolecular stimulus."""
    av = pd.read_csv(H.shared_dir() / 'data' / 'yang_weak_odor_kc' / 'orn_input_availability.csv')
    av = av[av['kind'] == 'monomolecular']
    inputs = session_inputs()
    rows = []
    for r in av.to_dict('records'):
        st = (str(r['abbrev']), float(r['log10_conc']))
        est = input_estimates(st, inputs)
        status = ('exact' if (est['status'] == 'exact').any() else 'substitute') if len(est) else 'none'
        est = est.assign(_neg=-est['input_conc']).sort_values(['log10_distance', '_neg'])
        closest = stim_label((st[0], float(est['input_conc'].iloc[0]))) if len(est) else ''
        want_closest = r['spike_deltas_closest'] if r['spike_deltas_status'] != 'none' else ''
        rows.append(dict(stimulus=r['odor'], status=status, want_status=r['spike_deltas_status'], closest=closest,
                         want_closest=want_closest,
                         ok=bool(status == r['spike_deltas_status'] and closest == want_closest)))
    df = pd.DataFrame(rows)
    return dict(passed=bool(df['ok'].all()), n=len(df), failures=df[~df['ok']].to_dict('records'))


# ---------------------------------------------------------------------------
# the flies: the stimulus table, the shuffled blocks, the solvent, the ratios (parts B to E; viki's copy)

def read_recordings(sets: Optional[Sequence[str]] = E100_SETS) -> pd.DataFrame:
    """viki's recordings.csv (`sets` None: every set), with `folder` (the recording's copied folder, absolute) and
    `released` (one of the 12 recordings of data/yang_apl_calyx_maps/, or YANG_RELEASED)."""
    path = copy_dir() / 'recordings.csv'
    if not path.exists():
        raise NotImplementedError(f"viki's copy of Yang's APL recordings is not in place ({path} missing)")
    r = pd.read_csv(path)
    if sets is not None:
        r = r[r['set'].isin(sets)].copy()
    r['folder'] = [str(H.shared_dir() / sp) if isinstance(sp, str) else '' for sp in r['shared_path']]
    rel = set(released_2021())
    r['released'] = [bool(YANG_RELEASED or st in rel) for st in r['stem']]
    return r.reset_index(drop=True)


def released_2021() -> List[str]:
    """The stems of the 12 recordings Yang released (data/yang_apl_calyx_maps/2021_APL_18um/recordings.csv)."""
    _, DS = _e79()
    return [f.split('_moco')[0] for f in pd.read_csv(DS.DATA / 'recordings.csv')['file']]


def odor_abbrev_yang(name: str) -> Tuple[str, float]:
    """One of Yang's odor_pin entries ('phenylacetaldehyde(-6.294)', '2h(-7)', 'PFO', 'paraffin oil') as (abbreviation,
    log10 concentration); the solvent as (SOLVENT, NaN); an unknown name is kept as written."""
    import re
    n = ' '.join(str(name).split())
    m = re.fullmatch(r'(?P<o>.+?)\s*\((?P<c>[-+0-9.eE]+)\)', n)
    odor, conc = (m.group('o'), float(m.group('c'))) if m else (n, np.nan)
    key = odor.strip().lower()
    if key in ('paraffin oil', 'pfo', 'solvent'):
        return SOLVENT, np.nan
    return ODOR_ABBREV.get(key, odor.strip()), conc


def yaml_stimuli(path: Path, recording: Optional[str] = None) -> pd.DataFrame:
    """The presentations of one recording from Yang's stimulus yaml: frame k is the k-th entry of odor_stim_seq (viki's
    README; the peak-map stack has one frame per presentation in that order), block k // nstim, position in the block
    k % nstim + 1; with the pin, Yang's odor name (odor_pin, 1-based pins), the abbreviation and log10 concentration
    (`odor_abbrev_yang`; the solvent as SOLVENT). Columns STIM_TABLE_COLUMNS plus position, pin, odor_name."""
    import yaml
    y = yaml.safe_load(Path(path).read_text())
    pins = list(y['odor_pin'])
    seq = [int(x) for x in y['odor_stim_seq']]
    nstim = int(y['nstim'])
    if len(seq) != int(y['nblock']) * nstim:
        raise ValueError(f"{path}: {len(seq)} presentations for nblock x nstim = {y['nblock']} x {nstim}")
    rows = []
    for k, pin in enumerate(seq):
        name = pins[pin - 1]
        odor, conc = odor_abbrev_yang(name)
        rows.append(dict(recording=recording or Path(path).stem.replace('_copy', ''), frame=k, block=k // nstim,
                         position=k % nstim + 1, pin=pin, odor_name=' '.join(str(name).split()), odor=odor,
                         log10_conc=conc))
    return pd.DataFrame(rows)


def recording_paths(row) -> dict:
    """One recording's files in viki's copy (the 12 released 2021 recordings' peak maps, outlines and yaml from
    data/yang_apl_calyx_maps/): the yaml, the APL peak-map stack (channel C1 or a single channel; the copied variant), the
    calyx outline matching the map's registration tag ('<stem>_<tag>.roi', else '<stem>.roi'), Yang's ROI table for the
    same tag, the Experiment.xml. Missing files are None."""
    stem, folder = row['stem'], Path(row['folder'])
    _, DS = _e79()
    if stem in released_2021():
        maps = sorted(DS.DATA.glob(f'{stem}_*dff_peak5frames_ave.tif'))
        yml = DS.DATA / f'{stem}.yaml'
        rois = [DS.DATA / f'{stem}_moco.roi']
    else:
        maps = sorted(f for f in folder.glob(f'{stem}_*dff_peak5frames_ave.tif') if '_C2' not in f.name)
        yml = next((f for f in (folder / f'{stem}.yaml', folder / f'{stem}_copy.yaml') if f.exists()), None)
        rois = []
    pm = maps[0] if maps else None
    tag = pm.name[len(stem) + 1:].split('_dff_peak5frames_ave')[0] if pm is not None else ''
    rois += [folder / f'{stem}_{tag}.roi', folder / f'{stem}.roi']
    roi = next((f for f in rois if f.exists()), None)
    tab = next((f for f in (folder / f'{stem}_{tag}_roi_peak_dff.csv',) if f.exists()), None)
    xml = next(iter(sorted(folder.glob(f'{stem}*Experiment.xml'))), None)
    return dict(yaml=yml, peak_maps=pm, tag=tag, roi=roi, roi_table=tab, experiment_xml=xml)


def read_recording(row, maps: bool = False) -> dict:
    """One recording read for parts B to D: its presentations (`yaml_stimuli`, `shuffled_blocks`' flags where the design
    has the solvent last) with each presentation's calyx response (`calyx_responses` of the APL peak maps inside the
    outline) and Yang's own value beside (`yang_peak`, from her ROI table, in presentation order), the outline mask, the
    pixel size and channels (`experiment_xml`), and with `maps` the stack. Refuses a recording Yang has not released
    (`YANG_RELEASED`)."""
    import tifffile
    if not (YANG_RELEASED or row['stem'] in released_2021()):
        raise PermissionError(f"{row['stem']}: Yang has not released this recording (E100's entry, Requires); set "
                              "YANG_RELEASED after her answer")
    _, DS = _e79()
    f = recording_paths(row)
    if f['peak_maps'] is None or f['yaml'] is None or f['roi'] is None:
        raise FileNotFoundError(f"{row['stem']}: peak maps, yaml or outline missing ({f})")
    stack = tifffile.imread(f['peak_maps']).astype(float)
    frames = yaml_stimuli(f['yaml'], row['stem'])
    if len(frames) != stack.shape[0]:
        raise ValueError(f"{row['stem']}: {stack.shape[0]} peak maps for {len(frames)} presentations")
    mask = DS.ACM.roi_mask(DS.ACM.read_imagej_polygon(f['roi']), stack.shape[1:])
    frames['response'] = calyx_responses(stack, mask)
    frames['is_solvent'] = frames['odor'] == SOLVENT
    frames['stim'] = [SOLVENT if sv else stim_label((o, c)) for sv, o, c in
                      zip(frames['is_solvent'], frames['odor'], frames['log10_conc'])]
    if f['roi_table'] is not None:
        t = pd.read_csv(f['roi_table'])
        if len(t) == len(frames):
            frames['yang_peak'] = t['roi_peak_dff'].to_numpy(dtype=float)
    out = dict(stem=row['stem'], set=row['set'], fly=f"{row['date']}_{row['fly']}", side=row['hemisphere'],
               plane_um=row['plane_depth_um'], pulse_s=row['pulse_s'], frames=frames, mask=mask, files=f,
               xml=experiment_xml(f['experiment_xml']) if f['experiment_xml'] is not None else None)
    if maps:
        out['stack'] = stack
    return out


def yang_table_identity(stems: Optional[Sequence[str]] = None) -> dict:
    """Part 0: on the released 2021 recordings (default all 12), the calyx-mean of each peak map inside the outline
    against Yang's ROI table (roi_peak_dff), per presentation: the largest absolute difference, the largest at odor and at
    solvent presentations, and Pearson r; passed when the largest difference is at most 1e-6 dF/F."""
    recs = read_recordings(sets=None)
    stems = released_2021() if stems is None else list(stems)
    rows = []
    for _, row in recs[recs['stem'].isin(stems)].iterrows():
        rec = read_recording(row)
        fr = rec['frames']
        if 'yang_peak' not in fr:
            rows.append(dict(recording=row['stem'], n=0))
            continue
        d = (fr['response'] - fr['yang_peak']).abs()
        rows.append(dict(recording=row['stem'], n=len(fr), max_abs_diff=float(d.max()),
                         max_abs_diff_odor=float(d[~fr['is_solvent']].max()),
                         max_abs_diff_solvent=float(d[fr['is_solvent']].max()),
                         pearson_r=float(np.corrcoef(fr['response'], fr['yang_peak'])[0, 1]),
                         median_ratio=float((fr['response'] / fr['yang_peak']).median())))
    df = pd.DataFrame(rows)
    worst = float(df['max_abs_diff'].max()) if 'max_abs_diff' in df else np.inf
    return dict(passed=bool(worst <= 1e-6), n_recordings=int(df['n'].gt(0).sum()), max_abs_diff_dff=worst,
                min_pearson_r=float(df['pearson_r'].min()) if 'pearson_r' in df else np.nan,
                per_recording=df.to_dict('records'))


def xml_identity(stems: Optional[Sequence[str]] = None) -> dict:
    """Part 0: each released 2021 recording's pixel size and channel list read from its Experiment.xml equal
    recordings.csv's (pixel_size_um, channels_xml)."""
    recs = read_recordings(sets=None)
    stems = released_2021() if stems is None else list(stems)
    rows = []
    for _, row in recs[recs['stem'].isin(stems)].iterrows():
        f = recording_paths(row)
        if f['experiment_xml'] is None:
            rows.append(dict(recording=row['stem'], ok=False, why='no Experiment.xml'))
            continue
        x = experiment_xml(f['experiment_xml'])
        ok = (abs(x['pixel_um'] - float(row['pixel_size_um'])) < 1e-6
              and '+'.join(x['channels']) == str(row['channels_xml']))
        rows.append(dict(recording=row['stem'], ok=bool(ok), pixel_um=x['pixel_um'], channels='+'.join(x['channels'])))
    return dict(passed=bool(rows and all(r['ok'] for r in rows)), per_recording=rows)


def read_stim_table(path: Optional[Path] = None) -> pd.DataFrame:
    """A table of presentations (STIM_TABLE_COLUMNS, e.g. the 2021 set's frames.csv), solvent labels mapped to SOLVENT;
    the copy's own stimulus tables are read per recording by `yaml_stimuli`."""
    _, DS = _e79()
    path = DS.DATA / 'frames.csv' if path is None else Path(path)
    t = pd.read_csv(path)
    t['odor'] = t['odor'].astype(str).str.strip()
    t.loc[t['odor'].str.lower().isin(SOLVENT_LABELS), 'odor'] = SOLVENT
    return t


def shuffled_blocks(frames: pd.DataFrame) -> pd.DataFrame:
    """One recording's presentations (frame order) with the shuffled-block design read: `block` (taken from the table, or,
    where absent, a new block after each solvent presentation, the solvent being last in its block), `position` (1-based
    in the block), `is_solvent`, `stim` (the stimulus label, SOLVENT for the solvent), `rep` (the stimulus's presentation
    index, 0-based) and the flags `solvent_last` (the block ends with its one solvent), `complete` (the block holds every
    stimulus of the recording once) and `duplicate` (a stimulus presented twice in one block)."""
    f = frames.sort_values('frame').reset_index(drop=True).copy()
    f['is_solvent'] = f['odor'].astype(str).str.lower().isin(SOLVENT_LABELS)
    f['stim'] = [SOLVENT if sv else stim_label((o, c)) for sv, o, c in zip(f['is_solvent'], f['odor'], f['log10_conc'])]
    if 'block' not in f.columns or f['block'].isna().all():
        # a block ends with its solvent: the block of a presentation is the number of solvents before it
        f['block'] = np.concatenate([[0], np.cumsum(f['is_solvent'].to_numpy()[:-1])]).astype(int)
    f['position'] = f.groupby('block').cumcount() + 1
    f['rep'] = f.groupby('stim').cumcount()
    stims = set(f.loc[~f['is_solvent'], 'stim'])
    flags = {}
    for b, g in f.groupby('block'):
        n_solv = int(g['is_solvent'].sum())
        last = bool(g['is_solvent'].iloc[-1])
        odors = g.loc[~g['is_solvent'], 'stim']
        flags[b] = dict(solvent_last=bool(last and n_solv == 1), complete=bool(set(odors) == stims and odors.is_unique),
                        duplicate=bool(odors.duplicated().any()))
    for k in ('solvent_last', 'complete', 'duplicate'):
        f[k] = f['block'].map({b: v[k] for b, v in flags.items()}).astype(bool)
    return f


def recording_solvent(resp: pd.DataFrame) -> float:
    """The recording's solvent response: the mean of its solvent presentations' responses (one per block). `resp` has
    is_solvent and response."""
    v = resp.loc[resp['is_solvent'], 'response']
    return float(v.mean()) if len(v) else np.nan


def block_solvent(resp: pd.DataFrame) -> pd.Series:
    """Each block's solvent response (the mean of its solvent presentations, one in the design; NaN for a block without
    one), by block."""
    blocks = pd.Index(sorted(resp['block'].unique()), name='block')
    return resp[resp['is_solvent']].groupby('block')['response'].mean().reindex(blocks)


def solvent_subtracted(resp: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per presentation, per stimulus): the odors' responses raw, minus the recording's solvent ('minus_recording') and
    minus their block's solvent ('minus_block'); per stimulus the mean over presentations of each, the number of
    presentations and the block-to-block SD of the raw response."""
    rs = recording_solvent(resp)
    bs = block_solvent(resp)
    o = resp[~resp['is_solvent']].copy()
    o['raw'] = o['response'].astype(float)
    o['minus_recording'] = o['raw'] - rs
    o['minus_block'] = o['raw'] - o['block'].map(bs).astype(float)
    g = o.groupby('stim')
    per_stim = pd.DataFrame(dict(raw=g['raw'].mean(), minus_recording=g['minus_recording'].mean(),
                                 minus_block=g['minus_block'].mean(), n_presentations=g.size(),
                                 sd_over_blocks=g['raw'].std(ddof=1)))
    per_stim['solvent_recording'] = rs
    return o, per_stim


def large_solvent(resp: pd.DataFrame, share: float = SOLVENT_LARGE_SHARE) -> bool:
    """Whether the recording's solvent response exceeds `share` of its largest stimulus response (mean over presentations
    of the raw responses)."""
    largest = resp[~resp['is_solvent']].groupby('stim')['response'].mean().max()
    return bool(recording_solvent(resp) > share * largest)


def fly_ratio_w(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The flies' ratio from per-recording numerator and denominator responses (columns fly, num, den; solvent already
    subtracted): recordings whose denominator is at most MIN_RESPONSE left out; the geometric mean over flies (each fly's
    recordings on the log scale) when every ratio is positive, else the arithmetic mean (each fly's recordings averaged
    arithmetically), with 'mean' saying which; the 90 % interval from `n_boot` bootstrap draws of flies (E79's
    `fly_bootstrap` resampling, seed `seed`); the counts."""
    d = per_recording[per_recording['den'] > MIN_RESPONSE].copy()
    d['ratio'] = d['num'] / d['den']
    n_rec, n_out = len(d), int(len(per_recording) - len(d))
    if n_rec == 0:
        return dict(ratio=np.nan, lo=np.nan, hi=np.nan, mean='none', n_flies=0, n_recordings=0, n_left_out=n_out)
    geometric = bool((d['ratio'] > 0).all())
    if geometric:
        per_fly = np.log(d['ratio']).groupby(d['fly']).mean().to_numpy()
    else:
        per_fly = d['ratio'].groupby(d['fly']).mean().to_numpy()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(per_fly), size=(n_boot, len(per_fly)))
    boot = per_fly[idx].mean(axis=1)
    est, lo, hi = per_fly.mean(), *np.percentile(boot, CI)
    if geometric:
        est, lo, hi = np.exp(est), np.exp(lo), np.exp(hi)
    return dict(ratio=float(est), lo=float(lo), hi=float(hi), mean='geometric' if geometric else 'arithmetic',
                n_flies=int(len(per_fly)), n_recordings=int(n_rec), n_left_out=n_out)


def calyx_responses(stack: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Per presentation (frame of a peak-map stack), the mean dF/F inside the calyx outline."""
    return np.asarray(stack, dtype=float)[:, mask].mean(axis=1)


def experiment_xml(path: Path) -> dict:
    """A ThorImage Experiment.xml: the pixel size (um, LSM's pixelSizeUM), the enabled channels (Wavelength names, and the
    LSM's channel bit mask where written) and the frame size (pixelX, pixelY)."""
    import xml.etree.ElementTree as ET
    root = ET.parse(path).getroot()
    lsm = root.find('.//LSM')
    if lsm is None:
        raise ValueError(f'{path}: no LSM element')
    chans = [w.get('name') for w in root.findall('.//Wavelengths/Wavelength')]
    return dict(pixel_um=float(lsm.get('pixelSizeUM')), channels=chans,
                channel_mask=int(lsm.get('channel')) if lsm.get('channel') is not None else None,
                pixel_x=int(lsm.get('pixelX')), pixel_y=int(lsm.get('pixelY')))


def reader_identity_2021(recording: Optional[str] = None) -> dict:
    """Part 0: on one 2021 recording (default the first of recordings.csv), the calyx means of `calyx_responses` with the
    order-matched, first-in-block and all-paraffin-oil baselines reproduce E79's E79_maps_across_baselines.csv (the
    per-stimulus mean of the presentations' calyx means, each presentation minus its baseline map's calyx mean)."""
    import tifffile
    _, DS = _e79()
    recs = pd.read_csv(DS.DATA / 'recordings.csv')
    rec = recs.iloc[0] if recording is None else recs[recs['file'].str.startswith(recording)].iloc[0]
    name = rec['file'].split('_moco')[0]
    stack = tifffile.imread(DS.DATA / rec['file']).astype(float)
    frames = DS.add_positions(pd.read_csv(DS.DATA / 'frames.csv').query('recording == @name')
                              .set_index('frame').sort_index())
    mask = DS.ACM.roi_mask(DS.ACM.read_imagej_polygon(DS.DATA / rec['roi']), stack.shape[1:])
    per = calyx_responses(stack, mask)
    b = DS.blank_frames(frames)
    base = dict(all_pfo=per[b['all']].mean(), first_in_block=per[b['first_in_block']].mean())
    by_pos = {2: per[b['pfo_pos2'][0]], 3: per[b['pfo_pos3'][0]]}
    ref = pd.read_csv(DS.OUT / 'E79_maps_across_baselines.csv').query('recording == @name').set_index('stimulus')
    odor = frames[(frames['odor'] != 'pfo') & ~frames['doubtful']]
    worst = 0.0
    n = 0
    for (o, c), g in odor.groupby(['odor', 'log10_conc']):
        lab = DS.stim_label((o, float(c)))
        if lab not in ref.index:
            continue
        got = dict(all_pfo=float(np.mean([per[i] - base['all_pfo'] for i in g.index])),
                   first_in_block=float(np.mean([per[i] - base['first_in_block'] for i in g.index])),
                   order_matched=float(np.mean([per[i] - by_pos[p] for i, p in zip(g.index, g['position'])])))
        for k, v in got.items():
            worst = max(worst, abs(v - float(ref.loc[lab, f'mean_{k}'])))
        n += 1
    return dict(passed=bool(n > 0 and worst <= 1e-9), recording=name, n_stimuli=n, max_abs_diff_dff=worst)


# ---------------------------------------------------------------------------
# the verdict rules

def ratio_over_estimates(per_estimate: Dict[str, str]) -> str:
    """One ratio's class from its comparison verdict under each session's input estimate: 'inside' if inside under any;
    'flatter' or 'steeper' only if that under all; otherwise 'not decided'."""
    v = list(per_estimate.values())
    if not v:
        return 'not decided'
    if 'inside' in v:
        return 'inside'
    for k in ('flatter', 'steeper'):
        if all(x == k for x in v):
            return k
    return 'not decided'


def w1_verdict(ratio_classes: Dict[str, str], min_ratios: int = W1_MIN_RATIOS) -> str:
    """W1 from the six ratios' classes: 'flatter', 'inside' or 'steeper' if at least `min_ratios` are; otherwise
    'not decided'."""
    v = list(ratio_classes.values())
    for k in ('flatter', 'inside', 'steeper'):
        if sum(x == k for x in v) >= min_ratios:
            return k
    return 'not decided'


def w2_verdict(comparison: str, n_flies: int, min_flies: int = W2_MIN_FLIES) -> str:
    """A W2 ratio's verdict: its comparison verdict with `min_flies` or more flies, otherwise 'described'."""
    return comparison if int(n_flies) >= min_flies else 'described'


def w3_pair_calyces(exceed: pd.DataFrame, share: float) -> int:
    """The number of calyces in which the prediction exceeds the null's 95th percentile in more than `share` of the
    recordings (`exceed`: one row per calyx and recording, column exceeds, bool)."""
    sh = exceed.groupby('calyx')['exceeds'].mean()
    return int((sh > share).sum())


def w3_verdict(exceed: pd.DataFrame, pairs: Optional[Sequence[str]] = None) -> dict:
    """W3 from one row per pair, calyx and recording (columns pair, calyx, recording, exceeds) on the pairs A4 kept
    (`pairs`, default every pair in `exceed`): 'passes' if at least W3_MIN_PAIRS pairs pass at W3_PASS_SHARE in at least
    W3_MIN_CALYCES calyces, 'fails' if no pair does at W3_FAIL_SHARE, else 'not decided' ('no pair' when A4 kept none);
    with the per-pair counts."""
    pairs = sorted(exceed['pair'].unique()) if pairs is None else list(pairs)
    rows = []
    for p in pairs:
        e = exceed[exceed['pair'] == p]
        rows.append(dict(pair=p, n_calyces_over_half=w3_pair_calyces(e, W3_PASS_SHARE),
                         n_calyces_over_quarter=w3_pair_calyces(e, W3_FAIL_SHARE)))
    n_pass = sum(r['n_calyces_over_half'] >= W3_MIN_CALYCES for r in rows)
    n_quarter = sum(r['n_calyces_over_quarter'] >= W3_MIN_CALYCES for r in rows)
    if not pairs:
        verdict = 'no pair'           # A4 left every pair out: W3 can neither pass nor fail
    else:
        verdict = 'passes' if n_pass >= W3_MIN_PAIRS else 'fails' if n_quarter == 0 else 'not decided'
    return dict(verdict=verdict, n_pairs=len(pairs), n_pairs_pass=int(n_pass), n_pairs_over_quarter=int(n_quarter),
                pairs=rows)


# ---------------------------------------------------------------------------
# A2: the family at the weak stimuli, read from E98's B1

E98 = 'E98_apl_response_strength'
B1_CONDS = ('on', 'acute_off')


def b1_dirs() -> List[Path]:
    """Where E98's B1 runs may be: its local outputs and its archive, with and without a B1 subfolder."""
    roots = [H.HERE / 'outputs' / E98 / 'B', H.shared_dir() / 'data' / 'model_outputs' / E98 / 'B']
    return [d for r in roots for d in (r / 'B1', r)]


def parse_b1_label(label: str) -> dict:
    """A B1 stimulus label read as odor, log10 concentration and source: '<odor> @ <conc>' (a panel odor's abbreviation
    with its concentration, or the panel's abbreviation alone: source 'panel'), '<odor> [<source>] @ <conc>' (E100's own
    estimate labels), '<odor> @ <conc> [<source>]' or '<odor> @ <conc> (<source>)'. Unparsed labels give odor None."""
    import re
    lab = str(label).strip()
    m = re.fullmatch(r'(?P<o>[^\[\]@()]+?)\s*\[(?P<s>[^\]]+)\]\s*@\s*(?P<c>[-+0-9.eE]+)', lab)
    if m is None:
        m = re.fullmatch(r'(?P<o>[^\[\]@()]+?)\s*@\s*(?P<c>[-+0-9.eE]+)\s*(?:[\[(](?P<s>[^\])]+)[\])])?', lab)
    if m is not None:
        return dict(label=lab, odor=m.group('o').strip(), conc=float(m.group('c')),
                    source=(m.group('s') or 'panel').strip())
    from kcq.harness import APL16_ODORS, ORN_INPUT_CONC
    if lab in APL16_ODORS:
        return dict(label=lab, odor=lab, conc=float(ORN_INPUT_CONC[lab]), source='panel')
    return dict(label=lab, odor=None, conc=np.nan, source=None)


def read_e98_b1(wirings: Sequence[str] = ALL_WIRINGS, dirs: Optional[Sequence[Path]] = None) -> pd.DataFrame:
    """E98's saved B1 runs: one row per wiring, condition ('on', 'acute_off') and B1 stimulus, with APL's window mean
    ('win'), pre-odor mean ('pre'), evoked activity ('evoked', win minus pre), the KCs' response rate ('response_rate', the
    fraction of KCs responding) and the parsed label (`parse_b1_label`). Read from outputs/E98_apl_response_strength/B/
    or its archive (data/model_outputs/E98_apl_response_strength/B/; `dirs` overrides), files `<wiring>_<cond>_apl.csv`
    (index stimulus, columns win, pre; E90's `_save_stage` layout, which E98's stub names as its source; a 'B1_' prefix
    is accepted) and `<wiring>_<cond>_responses.parquet` (KC x stimulus). E98's own loader (`load_b1`) is used instead once
    it has one. Raises NotImplementedError naming E98's B1 while they are not there."""
    if dirs is None:
        try:
            from experiments import E98_apl_response_strength as E98M
            loader = getattr(E98M, 'load_b1', None)
        except ImportError:
            loader = None
        if loader is not None:
            return loader(wirings)
        dirs = b1_dirs()
    rows = []
    found = set()
    for w in wirings:
        for cond in B1_CONDS:
            for d in dirs:
                f = next((d / f'{pre}{w}_{cond}_apl.csv' for pre in ('', 'B1_')
                          if (d / f'{pre}{w}_{cond}_apl.csv').exists()), None)
                if f is None:
                    continue
                apl = pd.read_csv(f, index_col=0)
                rf = f.with_name(f.name.replace('_apl.csv', '_responses.parquet'))
                resp = pd.read_parquet(rf) if rf.exists() else None
                rlab = None if resp is None else [c[-1] if isinstance(c, tuple) else str(c) for c in resp.columns]
                for stim, a in apl.iterrows():
                    rate = np.nan
                    if resp is not None and str(stim) in rlab:
                        rate = float(resp.iloc[:, rlab.index(str(stim))].astype(float).mean())
                    rows.append(dict(wiring=w, cond=cond, stimulus=str(stim), win=float(a['win']), pre=float(a['pre']),
                                     evoked=float(a['win'] - a['pre']), response_rate=rate, file=str(f),
                                     **{f'b1_{k}': v for k, v in parse_b1_label(stim).items() if k != 'label'}))
                found.add((w, cond))
                break
    if not rows:
        raise NotImplementedError(
            "E98's B1 runs (experiments/E98_apl_response_strength.py --part B --arm B1, atlas) are not there yet: none of "
            + ', '.join(str(d) for d in dirs) + " holds <wiring>_<on|acute_off>_apl.csv; A2 reads them once archived")
    df = pd.DataFrame(rows)
    df.attrs['missing'] = sorted(set((w, c) for w in wirings for c in B1_CONDS) - found)
    return df


def a2_stimuli() -> Tuple[Stim, ...]:
    """A2's stimuli: set W's, set S's and the panel's pa, ms and paa at 10^-3 (set H's concentrations are read from
    viki's copy when it is in place)."""
    out: List[Stim] = []
    for s_ in SET_W_STIMULI + SET_S_STIMULI:
        if s_ not in out:
            out.append(s_)
    return tuple(out)


def match_b1(b1: pd.DataFrame, stims: Sequence[Stim]) -> pd.DataFrame:
    """Per stimulus with ORN input (`input_estimates`) and estimate, the B1 stimulus with the same odor and source
    concentration and the same source (the session, or the panel for the APL16 input's own columns: megamat's megamat
    panel and validation2's validation2 panel); one row per stimulus, estimate and B1 stimulus (none: b1_stimulus NaN)."""
    lab = b1.drop_duplicates('stimulus')[['stimulus', 'b1_odor', 'b1_conc', 'b1_source']]
    rows = []
    inputs = session_inputs()
    for st in stims:
        for e in input_estimates(st, inputs).to_dict('records'):
            panel_own = (e['session'], e['panel_name']) in (('megamat', 'megamat'), ('validation2', 'validation2'))
            cand = lab[(lab['b1_odor'] == st[0]) & np.isclose(lab['b1_conc'].astype(float), e['input_conc'])]
            hit = cand[cand['b1_source'].astype(str).str.contains(e['session'], regex=False)]
            if hit.empty and panel_own:
                hit = cand[cand['b1_source'] == 'panel']
            rows.append(dict(e, b1_stimulus=hit['stimulus'].iloc[0] if len(hit) else np.nan))
    return pd.DataFrame(rows)


def part_a2() -> dict:
    """A2: a2_family.csv (per wiring, stimulus of sets W and S with ORN input, session estimate and condition: APL's
    window, pre-odor and evoked activity and the KCs' response rate, read from E98's B1) and a2_unmatched.csv (the
    stimuli with input that B1 does not carry). The family's W1 and W2 ratios through E98's indicator family are formed
    in `score`, beside the flies' intervals."""
    b1 = read_e98_b1()
    m = match_b1(b1, a2_stimuli())
    fam = m.dropna(subset=['b1_stimulus']).merge(b1, left_on='b1_stimulus', right_on='stimulus',
                                                 suffixes=('', '_b1')).drop(columns=['stimulus_b1'])
    out = out_dir()
    fam.to_csv(out / 'a2_family.csv', index=False)
    m[m['b1_stimulus'].isna()].to_csv(out / 'a2_unmatched.csv', index=False)
    return dict(n_rows=len(fam), n_unmatched=int(m['b1_stimulus'].isna().sum()), missing_runs=b1.attrs.get('missing'))


# ---------------------------------------------------------------------------
# A3: the connectome-predicted maps

def a3_stimuli() -> Tuple[Stim, ...]:
    """A3's stimuli: those of A3_PAIRS, then W1's paa at 10^-5 (its input described beside; no map)."""
    out: List[Stim] = []
    for a, b in A3_PAIRS:
        for s_ in (a, b):
            if s_ not in out:
                out.append(s_)
    return tuple(out) + (('paa', -5.0),)


def _column_labels(df: pd.DataFrame) -> List[str]:
    """Output column labels: the panel odors by abbreviation, the added columns by their full estimate label."""
    out = []
    for c in df.columns:
        o = c[-1] if isinstance(c, tuple) else str(c)
        out.append(o if '[' in o else H.odor_abbrev(c))
    assert len(set(out)) == len(out), out
    return out


def pn_rates_extended(stims: Optional[Sequence[Stim]] = None, rebuild: bool = False
                      ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """C2's PN stage (E79's `model_pn_rates` recipe: C2's simulator, APL off, the 0-0.5 s mean of pn_sims minus the
    -0.5-0 s mean) on the extended input (`extended_input` of `stims`, default `a3_stimuli()`): glomerulus x column evoked
    rates (the 16 panel odors by abbreviation, the added columns by their estimate labels) and the added columns' table;
    cached in a3_pn_rates.csv (with a3_pn_rates_total.csv) and a3_inputs.csv in the results folder. The PN stage does not
    depend on the wiring or on APL, so one run serves the eight calyces (C2's real-wiring fit, as E79)."""
    path, ipath = out_dir() / 'a3_pn_rates.csv', out_dir() / 'a3_inputs.csv'
    if path.exists() and ipath.exists() and not rebuild:
        return pd.read_csv(path, index_col=0), pd.read_csv(ipath)
    from experiments import E77_input_screen as S
    from pyolfsysm.mb_model import get_sim_var, get_time_index
    stims = a3_stimuli() if stims is None else tuple(stims)
    orn, table = extended_input(stims)
    P = S.c2_parts()
    t0 = time.time()
    # C2's own breadth: the APL->KC weights do not reach the PN stage, and handle_breadth would refit for a new input
    sim, info = S.build_c2_simulator(orn=orn, breadth=P['breadth'])
    r, _ = sim(dict(info['params'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    t = np.asarray(get_time_index(mp, rv))
    pn = get_sim_var(rv.pn, 'pn_sims')
    gl = list(dict.fromkeys(info['wPNKC'].columns.get_level_values('glomerulus')))
    if np.asarray(pn[0]).shape[0] != len(gl):
        raise RuntimeError(f'{np.asarray(pn[0]).shape[0]} PN rows for {len(gl)} glomeruli in the wiring')
    win = (t >= S.PN_WINDOW[0]) & (t <= S.PN_WINDOW[1])
    pre = (t >= -0.5) & (t < 0.0)
    tot = np.stack([np.asarray(pn[int(j)])[:, win].mean(axis=1) for j in order], axis=1)
    base = np.stack([np.asarray(pn[int(j)])[:, pre].mean(axis=1) for j in order], axis=1)
    labels = _column_labels(r)
    total = pd.DataFrame(tot, index=pd.Index(gl, name='glomerulus'), columns=labels)
    evoked = pd.DataFrame(tot - base, index=total.index, columns=labels)
    # the same run read with the step times pyolfsysm's get_time_index gave before 2026-11-13 (olfsysm 17cef71), which
    # E79's cached connectome_pn_rates.csv was read with (part 0's identity)
    tl = np.asarray(get_time_index(mp, rv, legacy_linspace=True))
    wl, pl = (tl >= S.PN_WINDOW[0]) & (tl <= S.PN_WINDOW[1]), (tl >= -0.5) & (tl < 0.0)
    legacy = pd.DataFrame(np.stack([np.asarray(pn[int(j)])[:, wl].mean(axis=1) - np.asarray(pn[int(j)])[:, pl].mean(axis=1)
                                    for j in order], axis=1), index=total.index, columns=labels)
    total.to_csv(out_dir() / 'a3_pn_rates_total.csv')
    legacy.to_csv(out_dir() / 'a3_pn_rates_legacy_time.csv')
    evoked.to_csv(path)
    table.to_csv(ipath, index=False)
    orn.to_csv(local_dir('A3') / 'orn_extended.csv')
    _write_json(out_dir() / 'provenance_A3_pn_rates.json',
                dict(script=H.script_path(__file__), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                     seconds=round(time.time() - t0, 1), stimuli=[stim_label(x) for x in stims],
                     family="C2 (E67's rate0.068_aplboth_k2_cap3), real wiring, APL off (APL->KC and APL->bouton scales 0)",
                     params=dict(info['params'], wAPLKC_scale=0.0, wAPLPN_scale=0.0), kws=info['kws'],
                     pn_window=list(S.PN_WINDOW), pre_window=[-0.5, 0.0],
                     model_provenance=H.model_provenance(orn_deltas=orn)))
    return evoked, table


def pn_identity(rates: Optional[pd.DataFrame] = None) -> dict:
    """Part 0: the 16 panel odors' PN stage in the extended run against the saved ones: (a) its 0-0.5 s window means equal
    E90's saved real-wiring stage runs (`real_on_pn_win.parquet`, the same family on the 16-odor input with the current
    code); (b) read with the legacy step times (a3_pn_rates_legacy_time.csv), its evoked rates equal E79's cached
    connectome_pn_rates.csv. Each the largest absolute difference (Hz), passed at 1e-9; and, for the record, the
    difference between the evoked rates A3 uses (current step times) and E79's cache."""
    CM, _ = _e79()
    rates = pd.read_csv(out_dir() / 'a3_pn_rates.csv', index_col=0) if rates is None else rates
    ref = pd.read_csv(CM.OUT / 'connectome_pn_rates.csv', index_col=0)

    def maxdiff(a, b):
        a = a.reindex(index=b.index, columns=b.columns)
        return float(np.abs(a.to_numpy(dtype=float) - b.to_numpy(dtype=float)).max()) if not a.isna().any().any() \
            else np.inf
    out = dict(n_glomeruli=int(ref.shape[0]), n_odors=int(ref.shape[1]),
               evoked_vs_e79_cache_hz=maxdiff(rates, ref))
    e90 = H.shared_dir() / 'data' / 'model_outputs' / 'E90_model_stages' / 'B' / 'real_on_pn_win.parquet'
    tot = out_dir() / 'a3_pn_rates_total.csv'
    if e90.exists() and tot.exists():
        out['window_vs_e90_hz'] = maxdiff(pd.read_csv(tot, index_col=0), pd.read_parquet(e90))
    leg = out_dir() / 'a3_pn_rates_legacy_time.csv'
    if leg.exists():
        out['legacy_evoked_vs_e79_cache_hz'] = maxdiff(pd.read_csv(leg, index_col=0), ref)
    checks = [out.get('window_vs_e90_hz', np.inf), out.get('legacy_evoked_vs_e79_cache_hz', np.inf)]
    out['passed'] = bool(all(x <= 1e-9 for x in checks))
    return out


def pair_maps(Wt: np.ndarray, perms: np.ndarray, G: np.ndarray, inside: np.ndarray) -> np.ndarray:
    """The predicted difference maps (permutations x 144) of a pair from the column x 2 weights `Wt` (the pair's two
    stimuli's PN rates per glomerulus), each permutation's column order and the columns' grids G: each stimulus's map over
    its mean inside the outline, a minus b (E79's `predicted_comparisons` for a pair)."""
    ma = Wt[:, 0][perms] @ G
    mb = Wt[:, 1][perms] @ G
    with np.errstate(invalid='ignore', divide='ignore'):
        return (ma / ma[:, inside].mean(axis=1, keepdims=True) - mb / mb[:, inside].mean(axis=1, keepdims=True))


def _pair_columns(rates: pd.DataFrame, inputs: pd.DataFrame) -> Dict[Tuple[str, str], Tuple[str, str]]:
    """(pair label, session) -> the two rate columns of the pair's stimuli under that session's estimate (the session's
    diagnostics; a pair is scored under a session only when both stimuli have an estimate there)."""
    out = {}
    for a, b in A3_PAIRS:
        for sess in SESSIONS:
            ca = inputs[(inputs['stimulus'] == stim_label(a)) & (inputs['session'] == sess)]['estimate']
            cb = inputs[(inputs['stimulus'] == stim_label(b)) & (inputs['session'] == sess)]['estimate']
            if len(ca) and len(cb):
                out[(pair_label(a, b), sess)] = (str(ca.iloc[0]), str(cb.iloc[0]))
    return out


def calyx_grids(calyx: str, parts: dict, rates: pd.DataFrame, n_perm: int = N_PERM) -> dict:
    """One calyx's bouton-driven column grids (E79's `slab_column_grids`) for the whole calyx and every slab with at least
    E79's share of the APL calyx synapses, the anatomy-alone grids (unweighted PN-to-APL synapses), the glomerulus names
    and the label permutations (row 0 the true labels; seeded by the calyx's position in the calyx list, so the same for
    every pair, slab and blur)."""
    CM, _ = _e79()
    frame = CM.calyx_frame(parts['apl'])
    centers, a0 = CM.slab_centers(parts['apl'])
    n_apl = len(parts['apl'])
    points, names, _ = CM._alternative_points(calyx, parts, 'bouton', rates, None)
    anat_points, _, _ = CM._alternative_points(calyx, parts, 'anatomy', rates, None)
    seed = SEED + 1000 * (1 + sorted(CALYCES).index(calyx) if calyx in CALYCES else 0)
    rng = np.random.default_rng(seed)
    perms = np.vstack([np.arange(len(names))] + [rng.permutation(len(names)) for _ in range(n_perm)])
    slabs = {}
    for center in [None] + list(centers):
        if center is not None:
            lo, hi = center - CM.SLAB_THICK_UM / 2, center + CM.SLAB_THICK_UM / 2
            if ((parts['apl']['A'] >= lo) & (parts['apl']['A'] < hi)).sum() < CM.SLAB_MIN_SHARE * n_apl:
                continue
        got = CM.slab_column_grids(points, len(names), parts['apl'], frame, center)
        if got is None:
            continue
        anat = CM.slab_column_grids(anat_points, 1, parts['apl'], frame, center)
        depth = CM.WHOLE_CALYX if center is None else float(center - a0)
        slabs[depth] = dict(G=got[0], inside=got[1], anat=None if anat is None else anat[0])
    return dict(names=names, perms=perms, slabs=slabs, side=CM.calyx_side(calyx), seed=seed)


CALYCES = ('flywire_L', 'flywire_R', 'malecns_L', 'malecns_R', 'banc_L', 'banc_R', 'hemibrain_R', 'model_hemibrain_R')
DUPLICATE_CALYCES = frozenset({'hemibrain_R', 'model_hemibrain_R'})   # the same fly's calyx, two bouton segmentations


def _wt(rates: pd.DataFrame, names: Sequence[str], cols: Sequence[str]) -> np.ndarray:
    return rates.reindex(list(names)).fillna(0.0)[list(cols)].to_numpy(dtype=float)


def _masked_r(x: np.ndarray, y: np.ndarray, inside: np.ndarray) -> float:
    CM, _ = _e79()
    return float(CM.masked_pearson(np.asarray(x, float)[None, :], np.asarray(y, float), inside)[0])


def load_a3_inputs(n_perm: int = N_PERM) -> dict:
    """The calyces (E79's `load_calyces`), the extended PN rates and A3's inputs table, and each calyx's grids."""
    CM, _ = _e79()
    rates, inputs = pn_rates_extended()
    calyces = CM.load_calyces()
    missing = sorted(set(CALYCES) - set(calyces))
    if missing:
        raise RuntimeError(f'calyces missing from E79: {missing}')
    t0 = time.time()
    grids = {}
    for c in CALYCES:
        grids[c] = calyx_grids(c, calyces[c], rates, n_perm)
        print(f'{time.strftime("%H:%M:%S")} {c}: {len(grids[c]["names"])} glomeruli, {len(grids[c]["slabs"])} slabs '
              f'({time.time() - t0:.0f} s)', flush=True)
    return dict(rates=rates, inputs=inputs, grids=grids, columns=_pair_columns(rates, inputs))


def part_a3(n_perm: int = N_PERM, state: Optional[dict] = None) -> dict:
    """A3 (module docstring): a3_maps.csv (per calyx, pair, session estimate, slab and blur: the SD of the predicted
    difference map over the grid cells inside the outline and its percentile in the label-permuted null, and the r of
    the difference map with the anatomy-alone map), a3_agreement.csv (per pair, estimate and blur, the mean r between the
    seven dataset calyces' whole-calyx difference maps, opposite sides mirrored, against the same with each calyx's own
    label permutation), a3_anatomy.csv (each stimulus's whole-calyx map against the anatomy-alone map), a3_sessions.csv
    (the r between session estimates' whole-calyx difference maps), a3_pred_grids.npz (the predicted grids, true labels)
    and a3_summary.json."""
    CM, _ = _e79()
    state = load_a3_inputs(n_perm) if state is None else state
    rates, grids, columns = state['rates'], state['grids'], state['columns']
    rows, anat_rows, keep = [], [], {}
    for c, g in grids.items():
        for depth, sl in g['slabs'].items():
            for blur, G in sl['G'].items():
                ins = sl['inside']
                anat = None if sl['anat'] is None else sl['anat'][blur][0]
                for (pair, sess), cols in columns.items():
                    D = pair_maps(_wt(rates, g['names'], cols), g['perms'], G, ins)
                    sd = np.nanstd(D[:, ins], axis=1)
                    keep[(c, pair, sess, depth, blur)] = (D[0].astype(np.float32), ins)
                    rows.append(dict(calyx=c, pair=pair, session=sess, slab_depth_um=depth, blur_um=blur,
                                     n_cells=int(ins.sum()), sd=float(sd[0]),
                                     null_pct_sd=float((sd[1:] < sd[0]).mean() * 100),
                                     null_median_sd=float(np.median(sd[1:])),
                                     r_with_anatomy=np.nan if anat is None else _masked_r(D[0], anat, ins)))
                if depth == CM.WHOLE_CALYX and anat is not None:
                    for est in sorted({x for cols in columns.values() for x in cols}):
                        m = _wt(rates, g['names'], [est])[:, 0] @ G
                        anat_rows.append(dict(calyx=c, estimate=est, blur_um=blur, r_with_anatomy=_masked_r(m, anat, ins)))
    maps = pd.DataFrame(rows)
    anat = pd.DataFrame(anat_rows)
    # agreement between the seven dataset calyces' whole-calyx difference maps (right side as reference)
    agr_rows = []
    data_calyces = [c for c in CALYCES if not c.startswith('model_')]
    for (pair, sess), cols in columns.items():
        for blur in CM.BLURS_UM:
            D = {}
            for c in data_calyces:
                g = grids[c]
                sl = g['slabs'][CM.WHOLE_CALYX]
                d = pair_maps(_wt(rates, g['names'], cols), g['perms'], sl['G'][blur], sl['inside'])
                ins = sl['inside']
                if g['side'] == 'L':
                    d, ins = CM.flip_grid(d, 'horizontal'), CM.flip_grid(ins, 'horizontal')
                D[c] = (d, ins)
            rs = []
            for i, a in enumerate(data_calyces):
                for b in data_calyces[i + 1:]:
                    ins = D[a][1] & D[b][1]
                    X, Y = D[a][0][:, ins], D[b][0][:, ins]
                    X = X - X.mean(axis=1, keepdims=True)
                    Y = Y - Y.mean(axis=1, keepdims=True)
                    with np.errstate(invalid='ignore', divide='ignore'):
                        rs.append((X * Y).sum(axis=1) / np.sqrt((X ** 2).sum(axis=1) * (Y ** 2).sum(axis=1)))
            m = np.nanmean(np.vstack(rs), axis=0)
            agr_rows.append(dict(pair=pair, session=sess, blur_um=blur, n_calyx_pairs=len(rs), mean_r=float(m[0]),
                                 null_pct=float((m[1:] < m[0]).mean() * 100), null_p95=float(np.quantile(m[1:], 0.95)),
                                 null_median=float(np.median(m[1:]))))
    agr = pd.DataFrame(agr_rows)
    # session estimates against each other (whole calyx, blur 4 um)
    ses_rows = []
    for c in CALYCES:
        for a, b in A3_PAIRS:
            pl = pair_label(a, b)
            ss = [s_ for s_ in SESSIONS if (c, pl, s_, CM.WHOLE_CALYX, 4.0) in keep]
            for i, s1 in enumerate(ss):
                for s2 in ss[i + 1:]:
                    x, ins = keep[(c, pl, s1, CM.WHOLE_CALYX, 4.0)]
                    y, _ = keep[(c, pl, s2, CM.WHOLE_CALYX, 4.0)]
                    ses_rows.append(dict(calyx=c, pair=pl, session_a=s1, session_b=s2, r=_masked_r(x, y, ins)))
    ses = pd.DataFrame(ses_rows)
    out = out_dir()
    maps.to_csv(out / 'a3_maps.csv', index=False)
    anat.to_csv(out / 'a3_anatomy.csv', index=False)
    agr.to_csv(out / 'a3_agreement.csv', index=False)
    ses.to_csv(out / 'a3_sessions.csv', index=False)
    keys = sorted(keep)
    np.savez_compressed(out / 'a3_pred_grids.npz', grids=np.stack([keep[k][0] for k in keys]),
                        inside=np.stack([keep[k][1] for k in keys]),
                        keys=np.array(['|'.join(map(str, k)) for k in keys]))
    whole = maps[maps['slab_depth_um'] == CM.WHOLE_CALYX]
    summ = dict(
        n_perm=int(n_perm), calyces=list(CALYCES), pairs=[pair_label(a, b) for a, b in A3_PAIRS],
        sessions=list(SESSIONS),
        whole_calyx_sd_null_pct_median=whole.groupby(['pair', 'blur_um'])['null_pct_sd'].median().unstack().round(1)
        .to_dict(),
        agreement=agr.groupby(['pair', 'blur_um'])[['mean_r', 'null_pct']].median().round(3).reset_index()
        .to_dict('records'),
        session_r_median=ses.groupby('pair')['r'].median().round(3).to_dict() if len(ses) else {},
        seeds={c: grids[c]['seed'] for c in CALYCES})
    _write_json(out / 'a3_summary.json', summ)
    return summ


# ---------------------------------------------------------------------------
# A4: the power of W3

def measured_reliability() -> Dict[str, float]:
    """The 2021 set's measured difference-map reliability: the median over recordings of E79's split-half r of pa minus ms
    and paa minus ms (connectome_measured.csv, E79's main measured baseline: every paraffin oil presentation)."""
    CM, _ = _e79()
    m = pd.read_csv(CM.OUT / 'connectome_measured.csv')
    return {c: float(m.loc[m['comparison'] == c, 'split_half_r'].median()) for c in ('paa-ms', 'pa-ms')}


def noise_sd(signal: np.ndarray, reliability: float, n_presentations: int = A4_N_PRESENTATIONS) -> float:
    """The SD of independent per-cell noise added to a true map `signal` (its inside cells) so that two single
    presentations would agree at Pearson r `reliability`, for the mean of `n_presentations` presentations:
    sd(signal) sqrt((1 - r) / r / n)."""
    r = float(reliability)
    return float(np.std(signal) * np.sqrt((1.0 - r) / r / n_presentations))


def scaled_reliability(r: float, f: float) -> float:
    """The split-half reliability of a map whose response is `f` times that of a map with reliability `r`, at the same
    noise (the signal's variance scales with f^2): f^2 r / (f^2 r + 1 - r)."""
    return float(f * f * r / (f * f * r + 1.0 - r))


def detect(pred: Dict[Tuple[float, float], np.ndarray], inside: Dict[Tuple[float, float], np.ndarray], meas: np.ndarray,
           meas_inside: np.ndarray, flips: Sequence[str] = D_FLIPS) -> Tuple[float, float, bool]:
    """One recording against one calyx: `pred` maps (slab depth, blur) to the predicted grids (row 0 the true labels, the
    rest permutations; permutations x 144), `inside` to their inside flags. The best r over slabs, blurs and flips, the
    95th percentile of the permutations' maxima over the same search, and whether the first exceeds the second."""
    CM, _ = _e79()
    mx = None
    for k, P in pred.items():
        for flip in flips:
            y = CM.flip_grid(meas, flip)
            ins = inside[k] & CM.flip_grid(meas_inside, flip)
            if ins.sum() < 3:
                continue
            r = CM.masked_pearson(P, y, ins)
            mx = r if mx is None else np.fmax(mx, r)
    if mx is None or not np.isfinite(mx[0]):
        return np.nan, np.nan, False
    p95 = float(np.nanquantile(mx[1:], W3_NULL_PCT / 100))
    return float(mx[0]), p95, bool(mx[0] > p95)


def _depth_slabs(slabs: dict, nominal: float, within: float = D_SLAB_WITHIN_UM) -> List[float]:
    CM, _ = _e79()
    return [d for d in slabs if d != CM.WHOLE_CALYX and abs(d - nominal) <= within + 1e-9]


def power(pair: str, session: str, state: dict, reliabilities: Sequence[float], n_draws: int = A4_N_DRAWS,
          seed: int = SEED, self_truth: bool = False) -> pd.DataFrame:
    """A4 for one pair under one session's input estimate: simulated recordings (a truth calyx's predicted map at a slab
    within D_SLAB_WITHIN_UM of a nominal plane depth from A4_PLANE_DEPTHS_UM and a blur from E79's, plus noise at each
    reliability, `noise_sd`) scored against every other calyx's prediction over its slabs within D_SLAB_WITHIN_UM of the
    nominal depth, E79's blurs and D_FLIPS (`detect`, the scored calyx mirrored into the truth's side); one row per scored
    calyx, truth calyx, reliability and draw. With `self_truth`, each calyx is also scored against its own truth (a
    control, flagged `self_truth`: the power without any difference between flies' calyces)."""
    CM, _ = _e79()
    rates, grids = state['rates'], state['grids']
    cols = state['columns'][(pair, session)]
    rng = np.random.default_rng([seed, zlib.crc32(f'{pair}|{session}'.encode())])
    # every calyx's predicted grids (all permutations) at the slabs any nominal depth can reach
    lo, hi = min(A4_PLANE_DEPTHS_UM) - D_SLAB_WITHIN_UM, max(A4_PLANE_DEPTHS_UM) + D_SLAB_WITHIN_UM
    P, INS = {}, {}
    for c, g in grids.items():
        Wt = _wt(rates, g['names'], cols)
        for d, sl in g['slabs'].items():
            if d == CM.WHOLE_CALYX or not (lo - 1e-9 <= d <= hi + 1e-9):
                continue
            for blur, G in sl['G'].items():
                P[(c, d, blur)] = pair_maps(Wt, g['perms'], G, sl['inside']).astype(np.float32)
                INS[(c, d)] = sl['inside']
    rows = []
    for t in CALYCES:
        t_side = grids[t]['side']
        for rel in reliabilities:
            for k in range(n_draws):
                nominal = float(rng.choice(A4_PLANE_DEPTHS_UM))
                cands = [d for d in _depth_slabs(grids[t]['slabs'], nominal) if (t, d) in INS]
                if not cands:
                    continue
                d_t = float(rng.choice(cands))
                b_t = float(rng.choice(sorted(grids[t]['slabs'][d_t]['G'])))
                truth = P[(t, d_t, b_t)][0].astype(float)
                ins_t = INS[(t, d_t)]
                sig = truth[ins_t]
                meas = truth + rng.normal(0.0, noise_sd(sig, rel), size=truth.shape)
                for c in CALYCES:
                    if (c == t and not self_truth) or {c, t} == DUPLICATE_CALYCES:
                        continue
                    mirror = grids[c]['side'] != t_side
                    pred, ins = {}, {}
                    for d in _depth_slabs(grids[c]['slabs'], nominal):
                        for blur in sorted(grids[c]['slabs'][d]['G']):
                            pr = P[(c, d, blur)]
                            pred[(d, blur)] = CM.flip_grid(pr, 'horizontal') if mirror else pr
                            ins[(d, blur)] = CM.flip_grid(INS[(c, d)], 'horizontal') if mirror else INS[(c, d)]
                    if not pred:
                        continue
                    best, p95, ex = detect(pred, ins, meas, ins_t)
                    rows.append(dict(pair=pair, session=session, scored_calyx=c, truth_calyx=t, self_truth=c == t,
                                     reliability=rel,
                                     draw=k, nominal_depth_um=nominal, truth_depth_um=d_t, truth_blur_um=b_t,
                                     best_r=best, null_p95=p95, exceeds=ex))
    return pd.DataFrame(rows)


def part_a4(state: Optional[dict] = None, n_draws: int = A4_N_DRAWS) -> dict:
    """A4 (module docstring): a4_power_draws.csv.gz (every simulated recording), a4_power.csv (per pair, session estimate
    and reliability: the share of simulated recordings that exceed the null's 95th percentile, pooled over scored calyces,
    truths and draws, and the number of scored calyces in which more than half / a quarter do, W3's counts),
    a4_power_by_calyx.csv and a4_w3_pairs.json (the pairs kept: power at least A4_POWER_MIN at the lower reliability under
    every session's estimate)."""
    state = load_a3_inputs() if state is None else state
    rel = measured_reliability()
    rels = sorted(rel.values())
    weak = scaled_reliability(rels[0], A4_RESPONSE_SCALE)
    t0 = time.time()
    draws = []
    for (pair, sess) in sorted(state['columns']):
        draws.append(power(pair, sess, state, rels + [weak], n_draws, self_truth=True))
        x = draws[-1]
        print(f'{time.strftime("%H:%M:%S")} A4 {pair} [{sess}]: share exceeding '
              f'{x[~x.self_truth].groupby("reliability")["exceeds"].mean().round(3).to_dict()}, own calyx as truth '
              f'{x[x.self_truth].groupby("reliability")["exceeds"].mean().round(3).to_dict()} ({time.time() - t0:.0f} s)',
              flush=True)
    D_all = pd.concat(draws, ignore_index=True)
    selfp = (D_all[D_all['self_truth']].groupby(['pair', 'session', 'reliability'])
             .agg(power_own_calyx=('exceeds', 'mean')).reset_index())
    D = D_all[~D_all['self_truth']]
    by_cal = (D.groupby(['pair', 'session', 'reliability', 'scored_calyx'])
              .agg(share=('exceeds', 'mean'), n=('exceeds', 'size'), median_best_r=('best_r', 'median'),
                   median_null_p95=('null_p95', 'median')).reset_index())
    pw = (D.groupby(['pair', 'session', 'reliability'])
          .agg(power=('exceeds', 'mean'), n=('exceeds', 'size'), median_best_r=('best_r', 'median'),
               median_null_p95=('null_p95', 'median')).reset_index())
    cnt = by_cal.groupby(['pair', 'session', 'reliability'])['share']
    pw = pw.merge(cnt.apply(lambda x: int((x > W3_PASS_SHARE).sum())).rename('n_calyces_over_half').reset_index())
    pw = pw.merge(cnt.apply(lambda x: int((x > W3_FAIL_SHARE).sum())).rename('n_calyces_over_quarter').reset_index())
    pw = pw.merge(selfp, how='left')
    role = {rels[0]: 'decisive (2021, paa minus ms)', rels[1]: 'beside (2021, pa minus ms)',
            weak: f'sensitivity: response x {A4_RESPONSE_SCALE:g}'}
    for t_ in (by_cal, pw):
        t_.insert(3, 'role', t_['reliability'].map(role))
    out = out_dir()
    D_all.to_csv(out / 'a4_power_draws.csv.gz', index=False)
    by_cal.to_csv(out / 'a4_power_by_calyx.csv', index=False)
    pw.to_csv(out / 'a4_power.csv', index=False)
    low = pw[pw['reliability'] == rels[0]]
    decisive = low.groupby('pair')['power'].min()
    kept = [p for p in decisive.index if decisive[p] >= A4_POWER_MIN]
    sens = pw[pw['reliability'] == weak].groupby('pair')['power'].min()
    res = dict(reliabilities=rel, decisive_reliability=rels[0], sensitivity_reliability=weak,
               sensitivity_power={p: round(float(v), 4) for p, v in sens.items()},
               own_calyx_power={p: round(float(v), 4) for p, v in
                                pw[pw['reliability'] == rels[0]].groupby('pair')['power_own_calyx'].min().items()},
               n_draws=int(n_draws), n_presentations=A4_N_PRESENTATIONS,
               plane_depths_um=list(A4_PLANE_DEPTHS_UM), slab_within_um=D_SLAB_WITHIN_UM, flips=list(D_FLIPS),
               power_min=A4_POWER_MIN, decisive_power={p: round(float(v), 4) for p, v in decisive.items()},
               kept=kept, left_out=[p for p in decisive.index if p not in kept],
               w3_can_pass=len(kept) >= W3_MIN_PAIRS, seconds=round(time.time() - t0, 1))
    _write_json(out / 'a4_w3_pairs.json', res)
    return res


# ---------------------------------------------------------------------------
# parts 0, B to E, score

def part_0(recordings: bool = False) -> dict:
    """Part 0 (module docstring): identity.json. With `recordings` False, only the identities that need no new
    recording (the PN stage, the 1-log matcher, the calyx-response reader on a 2021 recording)."""
    rates, _ = pn_rates_extended()
    res = dict(pn_stage=pn_identity(rates), availability=availability_identity(), reader_2021=reader_identity_2021())
    if recordings:
        # on the 12 released 2021 recordings only, until Yang releases the others (YANG_RELEASED)
        res['yang_tables'] = yang_table_identity()
        res['experiment_xml'] = xml_identity()
    res['passed'] = all(v['passed'] for v in res.values() if isinstance(v, dict) and v.get('passed') is not None)
    path = out_dir() / 'identity.json'
    old = _read_json(path) or {}
    old.update(res)
    _write_json(path, old)
    return res


def _waiting(part: str) -> None:
    raise NotImplementedError(
        f"part {part} reads Yang's recordings from viki's copy ({copy_dir()}) only after Yang has answered hal's "
        "questions (the entry's Requires; YANG_RELEASED); its readers are written and tested on made-up data and the 12 "
        "released 2021 recordings")


YANG_L_IS_FLY_RIGHT = True    # the entry: posterior view, dorsal up, "L" files the fly's right (until Yang answers)


def recording_table(frames: pd.DataFrame) -> dict:
    """Part B for one recording (`read_recording`'s frames, or made-up ones with block, is_solvent, stim, response):
    per presentation and per stimulus the responses raw, minus the recording's solvent and minus the block's solvent
    (`solvent_subtracted`); the recording's solvent response and whether it is large (`large_solvent`); the agreement
    between blocks (the mean over pairs of blocks of the Pearson r of their odor responses, raw); the drift (the slope per
    block of the solvent's response and of the mean odor response)."""
    per, stim = solvent_subtracted(frames)
    piv = per.pivot_table(index='stim', columns='block', values='raw')
    rs = [piv[a].corr(piv[b]) for i, a in enumerate(piv.columns) for b in piv.columns[i + 1:]]
    sol = frames[frames['is_solvent']].groupby('block')['response'].mean()
    odor = per.groupby('block')['raw'].mean()

    def slope(x):
        return float(np.polyfit(x.index.to_numpy(dtype=float), x.to_numpy(dtype=float), 1)[0]) if len(x) >= 2 else np.nan
    return dict(per_presentation=per, per_stim=stim, solvent=recording_solvent(frames), large_solvent=large_solvent(frames),
                block_agreement_r=float(np.nanmean(rs)) if rs else np.nan, n_blocks=int(frames['block'].nunique()),
                solvent_drift_per_block=slope(sol), odor_drift_per_block=slope(odor))


def flies_ratios(per_stim: pd.DataFrame, numerators: Sequence[Stim], denominator: Stim,
                 reading: str = 'minus_recording', n_boot: int = N_BOOT) -> pd.DataFrame:
    """The flies' ratios of each numerator over `denominator` from per-recording stimulus means (columns recording, fly,
    stim and `reading`; `recording_table`'s per_stim with the recording and fly added): per recording the two responses,
    then `fly_ratio_w` over flies; one row per numerator with the reading, the mean used and the counts."""
    rows = []
    den_lab = stim_label(denominator)
    for num in numerators:
        lab = stim_label(num)
        d = per_stim.pivot_table(index=['fly', 'recording'], columns='stim', values=reading)
        if lab not in d or den_lab not in d:
            rows.append(dict(numerator=lab, denominator=den_lab, reading=reading, n_flies=0))
            continue
        pr = d[[lab, den_lab]].dropna().reset_index().rename(columns={lab: 'num', den_lab: 'den'})
        rows.append(dict(numerator=lab, denominator=den_lab, reading=reading, **fly_ratio_w(pr, n_boot=n_boot)))
    return pd.DataFrame(rows)


def measured_pair_grid(stack: np.ndarray, frames: pd.DataFrame, mask: np.ndarray, a: str, b: str,
                       solvent: str = 'recording') -> Optional[dict]:
    """Part D's measured map of a pair (stimulus labels `a`, `b`) in one recording: each stimulus's mean peak map minus the
    solvent's map (the mean of the recording's solvent presentations, or with `solvent` 'block' each presentation minus its
    block's own), over its calyx mean, blurred E79's MEAS_BLUR_PX, on E79's 12 x 12 grid over the outline's bounding box,
    a minus b (image axes as recorded); with the inside cells, the calyx means and the split-half r (the pair's map from
    the even blocks against the odd blocks'). None when a stimulus is missing."""
    from scipy.ndimage import gaussian_filter
    CM, _ = _e79()
    f = frames.reset_index(drop=True)
    if not ((f['stim'] == a).any() and (f['stim'] == b).any()):
        return None
    sol = f[f['is_solvent']]
    rec_sol = stack[sol.index].mean(axis=0)
    blk_sol = {bk: stack[g.index].mean(axis=0) for bk, g in sol.groupby('block')}
    M_, inside = CM.cell_matrix(mask)

    def stim_map(lab, blocks=None):
        g = f[(f['stim'] == lab) & (True if blocks is None else f['block'].isin(blocks))]
        if g.empty:
            return None, np.nan
        reps = [stack[i] - (rec_sol if solvent == 'recording' else blk_sol.get(bk, rec_sol)) for i, bk in
                zip(g.index, g['block'])]
        m = np.mean(reps, axis=0)
        return m, float(m[mask].mean())

    def grid(lab, blocks=None):
        m, mean = stim_map(lab, blocks)
        if m is None or not mean:
            return None, mean
        return np.asarray(M_.T @ (gaussian_filter(m, CM.MEAS_BLUR_PX) / mean).ravel()).ravel(), mean
    ga, ma = grid(a)
    gb, mb = grid(b)
    if ga is None or gb is None:
        return None
    blocks = sorted(f['block'].unique())
    halves = [blocks[0::2], blocks[1::2]]
    hg = []
    for hb in halves:
        x, _ = grid(a, hb)
        y, _ = grid(b, hb)
        hg.append(None if x is None or y is None else x - y)
    sh = (_masked_r(hg[0], hg[1], inside) if all(h is not None for h in hg) else np.nan)
    return dict(grid=ga - gb, inside=inside, mean_a=ma, mean_b=mb, split_half_r=sh)


def score_recording(meas: dict, side: str, plane_um: float, pair: str, session: str, state: dict) -> pd.DataFrame:
    """Part D for one recording and pair under one session's estimate: per calyx (E79's eight), the best r over its slabs
    within D_SLAB_WITHIN_UM of the plane's nominal depth, E79's blurs and D_FLIPS, against the label-permuted null's 95th
    percentile over the same search (`detect`); the calyx's prediction mirrored into the recording's hemisphere (Yang's
    "L" read as the fly's right while YANG_L_IS_FLY_RIGHT). One row per calyx with W3's `exceeds`."""
    CM, _ = _e79()
    rates, grids = state['rates'], state['grids']
    cols = state['columns'][(pair, session)]
    fly_side = side
    if YANG_L_IS_FLY_RIGHT:
        fly_side = {'L': 'R', 'R': 'L'}.get(side, side)
    rows = []
    for c in CALYCES:
        g = grids[c]
        Wt = _wt(rates, g['names'], cols)
        mirror = g['side'] != fly_side
        pred, ins = {}, {}
        for d in _depth_slabs(g['slabs'], float(plane_um)):
            sl = g['slabs'][d]
            for blur, G in sl['G'].items():
                pr = pair_maps(Wt, g['perms'], G, sl['inside'])
                pred[(d, blur)] = CM.flip_grid(pr, 'horizontal') if mirror else pr
                ins[(d, blur)] = CM.flip_grid(sl['inside'], 'horizontal') if mirror else sl['inside']
        best, p95, ex = detect(pred, ins, meas['grid'], meas['inside']) if pred else (np.nan, np.nan, False)
        rows.append(dict(calyx=c, pair=pair, session=session, n_slabs=len({k[0] for k in pred}), best_r=best,
                         null_p95=p95, exceeds=ex))
    return pd.DataFrame(rows)


def w3_over_estimates(per_estimate: Dict[str, str]) -> str:
    """W3 over the session estimates: 'passes' or 'fails' only if so under every estimate, otherwise 'not decided'."""
    v = list(per_estimate.values())
    for k in ('passes', 'fails'):
        if v and all(x == k for x in v):
            return k
    return 'not decided'


def part_b() -> dict:
    """Part B (module docstring), from viki's copy, after Yang's answer (`recording_table` per recording,
    `flies_ratios` for W1)."""
    if not YANG_RELEASED:
        _waiting('B')
    raise NotImplementedError('part B: the run over sets W, S and H is written after Yang answers')


def part_c() -> dict:
    """Part C (module docstring), from viki's copy, after Yang's answer (E79's analyses 1 to 3 on the stimulus maps)."""
    _waiting('C')


def part_d() -> dict:
    """Part D (module docstring), from viki's copy, after Yang's answer (`measured_pair_grid`, `score_recording`,
    `w3_verdict` per estimate, `w3_over_estimates`)."""
    _waiting('D')


def part_e() -> dict:
    """Part E (module docstring), only if set P is processed."""
    _waiting('E')


def score() -> dict:
    """The verdict (W1, W2, W3) from parts A2, A4, B and D: verdict.json."""
    _waiting('score')


PARTS = ('0', 'A2', 'A3', 'A4', 'B', 'C', 'D', 'E', 'score')
READS = {
    '0': "E79's connectome_pn_rates.csv and E79_maps_across_baselines.csv, data/yang_weak_odor_kc/"
         "orn_input_availability.csv, one 2021 recording of data/yang_apl_calyx_maps/, A3's a3_pn_rates.csv",
    'A2': "E98's B1 runs (outputs/E98_apl_response_strength/B/ or its archive), the lab's three ORN sessions",
    'A3': "C2's simulator (E77's build_c2_simulator) on the extended input, E79's calyces (Yang's Claude's connectome "
          "tables in Tom/claude_APL_modeling_YZ4TC and Pratyush's boutons)",
    'A4': "A3's inputs (rebuilt in the process), E79's connectome_measured.csv",
}
WRITES = {
    '0': 'identity.json', 'A2': 'a2_family.csv, a2_unmatched.csv',
    'A3': 'a3_pn_rates.csv, a3_pn_rates_total.csv, a3_inputs.csv, a3_maps.csv, a3_agreement.csv, a3_anatomy.csv, '
          'a3_sessions.csv, a3_pred_grids.npz, a3_summary.json; outputs/E100_apl_weak_odors/A3/orn_extended.csv',
    'A4': 'a4_power_draws.csv.gz, a4_power.csv, a4_power_by_calyx.csv, a4_w3_pairs.json',
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS + ('A3A4',),
                    help="A3A4 runs A3 and then A4 in one process (the calyces' grids built once)")
    ap.add_argument('--no-recordings', action='store_true', help="part 0: the identities that need no new recording")
    ap.add_argument('--n-perm', type=int, default=N_PERM, help='A3 / A4: label permutations')
    ap.add_argument('--n-draws', type=int, default=A4_N_DRAWS, help='A4: noise draws per scored and truth calyx')
    ap.add_argument('--rebuild-pn', action='store_true', help='A3: rerun the PN stage even if a3_pn_rates.csv exists')
    ap.add_argument('--dry-run', action='store_true', help='say what the part reads and writes, and exit')
    a = ap.parse_args()
    parts = ('A3', 'A4') if a.part == 'A3A4' else (a.part,)
    if a.dry_run:
        for p_ in parts:
            print(f'part {p_}: reads {READS.get(p_, "viki copy of Yang recordings (after Yang answers)")}; writes '
                  f'{WRITES.get(p_, "not written yet")} under {out_dir()}')
        return
    if a.part == '0':
        provenance('0', recordings=not a.no_recordings)
        print(json.dumps(part_0(recordings=not a.no_recordings), indent=1, default=float), flush=True)
    elif a.part == 'A2':
        provenance('A2')
        print(json.dumps(part_a2(), indent=1, default=float), flush=True)
    elif a.part in ('A3', 'A4', 'A3A4'):
        if a.rebuild_pn:
            pn_rates_extended(rebuild=True)
        state = load_a3_inputs(a.n_perm)
        if 'A3' in parts:
            provenance('A3', n_perm=a.n_perm)
            print(json.dumps(part_a3(a.n_perm, state), indent=1, default=float), flush=True)
        if 'A4' in parts:
            provenance('A4', n_perm=a.n_perm, n_draws=a.n_draws)
            print(json.dumps(part_a4(state, a.n_draws), indent=1, default=float), flush=True)
    else:
        {'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e, 'score': score}[a.part]()


if __name__ == '__main__':
    main()
