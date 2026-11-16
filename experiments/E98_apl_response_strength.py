#!/usr/bin/env python3
"""E98: the strength of APL's response in the flies against the model's APL. In the best family APL's activity is
nearly the same for every odor of the panel, and two fifths of it is there before the odor; whether Yang's 2021 calyx
maps, read as each fly's APL response by odor and concentration, agree with that APL wherever the two can be compared,
and which recordings would tell a flat APL from one that follows the KCs' response. A diagnostic: it moves no loss.

The entry is item E98 of experiment_priorities.md (proposed by atlas 2026-11-15 12:45; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-15 12:26). C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the
0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit` / `build_simulator`), nothing
re-tuned. The flies: Yang's 2021 single-plane APL calyx maps (data/yang_apl_calyx_maps/2021_APL_18um/ in the shared
directory; 12 recordings, 8 flies; GCaMP6f in APL; pa, paa and ms as concentration series; 7 recordings with a 1 s
pulse, 5 with 0.5 s), read with E79's reader, baselines, blank and fly bootstrap (experiments/E79_apl_calyx_constraints.py).

Definitions (the entry's, fixed before any reading):
- A recording's response to a stimulus (the carrying reading, MAIN_READING; the entry's (C) as amended 2026-11-15 15:24,
  which struck "with the recording's position-matched blank subtracted"): the calyx-mean peak dF/F inside the outline
  under E79's order-matched baseline, the mean of the two presentations, nothing further subtracted. A recording enters
  a ratio only if both responses exceed MIN_RESPONSE (0.02).
- A concentration that was not imaged: linear interpolation of the response on log10 concentration between its two
  imaged neighbours, each recording inside its own imaged range only; for a concentration that no recording brackets
  (V3 to V5), the nearest imaged concentration stands in when it is within MATCH_LOG_UNITS (1, inclusive; Tom's rule)
  of the one asked for (`value_at`). 20210805_3FRfun1 and 20210823_4FRfun2 (no ms at 10^-2) are therefore out of V1.
- The flies' ratio for two stimuli: the geometric mean over the 8 flies (a fly's recordings averaged on the log scale),
  with a 90 % interval from N_BOOT (10,000) bootstrap draws of flies (E79's `fly_bootstrap`, seed SEED).
- The model's APL for a stimulus: the mean of APL's activity in the odor window (a; E90's MODEL_WINDOW, 0-0.5 s) and its
  pre-odor mean (a0), as E90's stage runs save them (`<wiring>_on_apl.csv`).
- The model's ratio through an indicator of exponent n (INDICATOR_NS, 1 to 2.5 in steps of 0.25), in two forms
  (`indicator_ratio`): 'baseline', (a1^n - a0^n) / (a2^n - a0^n), the pre-odor level in the indicator's baseline; and
  'evoked', ((a1 - a0) / (a2 - a0))^n, only the evoked part seen. GCaMP6f's Hill coefficient is about 2.3 (a literature
  value: Chen et al. 2013, Nature 499:295).
- A wiring is inside a comparison if any of its ratios over n and the two forms lies in the flies' interval; flatter if
  all of them are nearer 1 than the interval or on the other side of 1; steeper if all are beyond the interval's far end
  (`classify`). A comparison is inside, flatter or steeper if at least VERDICT_MIN_WIRINGS (6) of the 8 wirings are;
  otherwise not decided (`comparison_verdict`). The linear reading (n 1, evoked) is given beside every verdict.

Parts (CLI `--part`):
- 0, the identities (identity.json; part A refuses to run until the reader's pass): the reader reproduces E79's
  per-stimulus calyx means under its three baselines (E79_maps_across_baselines.csv) and its blank table
  (E79_baseline_noise.csv, E79_blank_vs_odor_maps.csv); Hallem's 10^-2 rows through olfsysm's hc_data.csv equal E38's
  (`hallem_identity`; olfdb's copy of the table checked beside when present); the extended-input run reproduces each
  wiring's saved run on the 16 panel odors with 0 count differences and E90's saved APL window and pre-odor means (real
  wiring: mean over odors 19.548, pre-odor share 0.423) ('reproduces', which part B requires); a blank stimulus gives
  an evoked APL of zero ('blank_zero': it does not with E90's pre-odor window, which includes the KC layer's start-up
  from zero at -0.5 s; 0.51 on the real wiring, -0.016 over SETTLED_PRE_WINDOW; part C's pre-odor level, `with_a0`,
  waits on a decision) (`--no-runs`: the reader and Hallem only).
- A (no simulations): A1 the flies' responses per recording and stimulus under E79's three baselines, with and without
  the blank subtracted, the two presentations' agreement, the block order and the drift of the first-in-block paraffin
  oil response (a1_responses.csv, a1_agreement.csv, a1_pfo_drift.csv); A2 the dose-response per odor, recording and fly:
  the fold change per log unit, the lowest concentration above the blank, the values at the model inputs' concentrations
  (a2_dose_response.csv, a2_at_model_concs.csv, a2_medians.csv); A3 the fixed list of ratios (RATIOS: V1 to V5 and the
  entry's probe rows) with their intervals, under every reading and the sensitivity subsets (SUBSETS: the 1 s and the
  0.5 s recordings apart; without the two recordings with a large blank and 20210805_3FRfun1) (a3_ratios.csv,
  a3_per_fly.csv); A4 the family's saved runs (E90's stage runs): APL's window mean, pre-odor mean and evoked activity
  per odor on the eight wirings, the profile over the 16 odors against each odor's input total and PN total, and the
  model's ratios of the fixed list that the panel's odors give, through the indicator family (a4_apl.csv,
  a4_profile.csv, a4_model_ratios.csv, a4_preview.csv); A5 what other APLs would show for the same ratios: one in
  proportion to the KCs' response rates in each genotype, one in proportion to the PN dendrites' totals (the megamat
  odors, E90's tables), E37's saved table for the spike-driven fits (a5_alternatives.csv); A6 the power table: for every
  pair of panel odors the ratio by the family's APL and by a proportional APL, and the number of flies that separates the
  two (a6_power.csv, the five that need the fewest in a_summary.json). a_summary.json.
- B (simulations at fixed parameters, APL on and acutely off, each saved with its counts, responses, APL means and PN
  stage; `--root` puts the runs under another folder, e.g. the archive's from_hal/ for hal's lane):
  - B1 (`run_b1`, `b1_specs`, `b1_input`): the extended input on the eight wirings (the real wiring first; each run checked
    by part 0's run identities), on the input's 36 glomeruli, 52 stimuli with fixed labels: the 16 panel odors; paa at
    10^-5 from the validation2 session's diagnostics (V2), the megamat and kiwi-control sessions' beside; ms at 10^-3
    from the three sessions' diagnostics and the control panel (V1's other estimates of ms's input); 2h, t2h, 2-but and
    1-6ol at 10^-6, 2h at 10^-5, ma at 10^-7 and ea at 10^-8 from the kiwi-control session (LOW_SESSION, E49's source);
    E100's 12 (4h-ol at 10^-6 and farn at 10^-2 from all three sessions, 2h and 2-but at 10^-6 and ma at 10^-7 from
    the other two; read by none of E98's parts); Hallem's pa and ms at 10^-2, -4, -6 and -8 and paa at 10^-2 (olfsysm's
    pyolfsysm/data/hc_data.csv, `hallem_rows`, on the glomeruli it shares with the input; part 0 checks its 10^-2 rows
    against E38's drosolf table, `hallem_identity`); the blank. The saved layout and its reader, `load_b1` (which E100
    uses), are described in B1_LAYOUT. With B1 the KCs' response rate to the blank, APL on and off, against the flies'
    solvent.
    E101's part A3 ran B1's 52 stimuli on all eight wirings (hal, 2026-11-15) through `ExtendedRun`; `load_b1_all` takes
    the wirings E98 has not run from its archived runs (data/model_outputs/E101_apl_vs_kc_claws/A3/, E98's layout; every
    file's sha256 and E101's commits recorded), after checking that they equal E98's own runs exactly where both exist
    (the real wiring).
  - B2 to B4 on the real wiring (`run_point_b`, `point_spec`, POINT_NAMES), each on B1's stimuli (the runs of
    2026-11-15 on the 46 without Hallem's rows below 10^-2, which came after) and checked on the 16
    panel odors against its source's saved run (0 differences): B2 C2 with the 1 s pulse (E80's C2_pulse1_fixed; APL's
    mean also over 0-1 s, 'win_pulse'), E93's fuller input F at its held point B_F_real (`e93_input`: B1's stimuli on
    46 glomeruli), E95's refit at the connectome's bouton share with the uniform component (f1); B3 E36's KC-spike-driven
    fit and E37's spike-only point at 40 Hz with the gain; B4 E77's screen points S1_inhadd1000 and S1_inhadd3141 (weaker
    PN normalization) at their held readings.
- C (no simulations; the entry's (C) as amended 2026-11-15 15:24; `verdict_core`): every comparison read under the
  carrying reading, its three partners (PARTNERS: the first-in-block baseline with no blank; the two recordings without
  ms at 10^-2 read at ms 10^-3.431, V1 only; without the two recordings with a large blank and 20210805_3FRfun1) and
  the entry's literal reading (the blank removed twice, decides nothing); a category is held when the carrying reading
  and every partner give it (`held_category`), otherwise reading-dependent (not decided). V1 (paa over ms near 10^-3,
  the lab's input, the hull of the two readings of paa's concentration, under each estimate of ms's input, V1_ESTIMATES,
  the panel's own first; `v1_from_estimates`) and V2 (paa at 10^-5 over paa at 10^-4, in the series as recorded); the
  pulse sets' permutation test (`pulse_permutation`, the 56 splits of the 8 flies into 3 and 5; disagreeing sets leave
  the comparison not decided); V2's descriptions (the order factor, the shuffled recording, the pa 10^-6.771 repeat, the
  paraffin oil drift; no gate); the pulse-matched readings; V3 to V5 (Hallem's input) described; the outcomes (a) to
  (e) and (h) (verdict.json, c_comparisons.csv, c_wirings.csv with each wiring's distance from the nearer end).
- table: comparison_table.md and the figures.

Deviations from the entry (each also in the note):
- Resolved by the entry's amendment (2026-11-15 15:24): the first version's main reading subtracted the position-matched
  blank from the order-matched responses, removing the position effect twice (the first-in-block baseline with the blank
  subtracted equals the order-matched baseline without it, `blank_identity`); the carrying reading is now the
  order-matched one without it, and the literal reading is shown beside and decides nothing. The one-log rule on the
  series' design range only (two recordings out of V1) is now the entry's rule; A3's subset 'all_recording_ranges' is
  one of the partners.
- E79's reader leaves out the doubled pa 10^-6.771 block of 20210805_3FRfun1 and 20210823_4FRfun2 (Yang, 2026-11-13: both
  are pa 10^-6.771 repeats); kept out here so that the reader reproduces E79, and no ratio of the fixed list reads it.
- B1's low-concentration stimuli come from the kiwi-control session's diagnostics (E49's source for ma and ea), with the
  other sessions' beside where the entry asks for them; paa at 10^-5 from validation2's, as the entry says.
- Hallem's conversion: the rows on the 16 glomeruli Hallem and the input share, 0 on the other 20 (E38's
  'hallem_only'), read from olfsysm's hc_data.csv (its sha256 and the olfsysm commit in identity.json and the runs'
  rows), DM3 from Or47a (the file gives dm3 to Or33b too).
- E93's F input is read from its saved file (inputs/F.csv, both hashes checked): E93's `input_of('F')` now rebuilds it
  within 1.4e-14, not exactly, and the identity against B_F_real needs the exact input.
- E37's saved APL table has no pre-odor level and is the mean (and peak) over the whole time after odor onset, not over
  0-0.5 s, so its ratios in A5 are of the window values, not of the evoked activity.

Outputs: results/E98_apl_response_strength/ (the tables above, identity.json, provenance_<part>.json, comparison_table.md,
figures); counts, responses and APL means of part B under outputs/E98_apl_response_strength/, archived to
data/model_outputs/E98_apl_response_strength/.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part 0 --no-runs
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part A
    env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E98_apl_response_strength --part 0 --wirings real --threads 4      # B1 on the real wiring
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E98_apl_response_strength --part B --arm B1 --wirings P1,P2,M1,P3,P4,M2,M3 --threads 16 \\
        --root <shared>/data/model_outputs/E98_apl_response_strength/from_hal                # another host's lane
    env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E98_apl_response_strength --part B --arm B2,B3,B4 --threads 4
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part 0 --wirings all8 --threads 1 ... (reads saved runs)
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part C --reading <reading> --subset <all|all_recording_ranges>
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E98_apl_response_strength'
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's ALL_WIRINGS (S then T)

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any reading)

MIN_RESPONSE = 0.02                    # a recording enters a ratio only if both responses exceed this (dF/F)
MATCH_LOG_UNITS = 1.0                  # Tom's rule: a concentration within one log10 unit (inclusive) stands in
N_BOOT = 10000                         # bootstrap draws of flies
SEED = 0
CI = (5.0, 95.0)                       # 90 % intervals
INDICATOR_NS = tuple(float(x) for x in np.arange(1.0, 2.5 + 1e-9, 0.25))
FORMS = ('baseline', 'evoked')
VERDICT_MIN_WIRINGS = 6                # of the 8 wirings
BASELINES = ('order_matched', 'first_in_block', 'all_pfo')          # E79's three
MAIN_READING = ('order_matched', False)            # (baseline, blank subtracted): the carrying reading (amended 15:24)
READINGS = (MAIN_READING, ('order_matched', True), ('first_in_block', False), ('first_in_block', True),
            ('all_pfo', False), ('all_pfo', True))
LARGE_BLANK = 0.05                     # a recording's blank calyx mean above this is large (E79: 0.09 and 0.18, the rest ~0)
NOT_IN_18UM_LIST = ('20210805_3FRfun1',)
SUBSETS = ('all', 'pulse_1s', 'pulse_0.5s', 'without_large_blank_and_3FRfun1', 'all_recording_ranges')
POWER_MIN_CONC = -5.0                  # A6's between-fly SD: A3's across-odor ratios with both stimuli at this or above
POWER_Z = (1.6449, 0.8416)             # A6: a two-sided 90 % interval excluding the other value with 80 % power

# the fixed list of ratios: (id, numerator, denominator, group); a stimulus is (odor, log10 concentration)
Stim = Tuple[str, float]
RATIOS: Tuple[Tuple[str, Stim, Stim, str], ...] = (
    ('V1_paa-4', ('paa', -4.0), ('ms', -3.0), 'V1'),
    ('V1_paa-2.954', ('paa', -2.954), ('ms', -3.0), 'V1'),
    ('V2', ('paa', -5.0), ('paa', -4.0), 'V2'),
    ('V3_ms-4', ('ms', -4.0), ('ms', -2.0), 'V3'),
    ('V3_ms-6', ('ms', -6.0), ('ms', -2.0), 'V3'),
    ('V3_ms-8', ('ms', -8.0), ('ms', -2.0), 'V3'),
    ('V4', ('pa', -6.0), ('pa', -4.0), 'V4'),
    ('V5_-4', ('pa', -4.0), ('ms', -4.0), 'V5'),
    ('V5_-6', ('pa', -6.0), ('ms', -6.0), 'V5'),
    ('V5_paa-2', ('paa', -2.0), ('ms', -2.0), 'V5'),
    ('probe_pa-4.386', ('pa', -4.386), ('ms', -3.0), 'probe'),
    ('probe_-4.863', ('pa', -4.863), ('ms', -4.863), 'probe'),
)
# the planner's probes (the entry's table; E79's order-matched means, no blank subtracted)
PROBES: Dict[str, Tuple[float, float, float]] = {
    'V1_paa-4': (0.85, 0.70, 0.97), 'V1_paa-2.954': (1.14, 1.03, 1.25), 'probe_pa-4.386': (1.14, 1.05, 1.24),
    'V2': (0.53, 0.46, 0.61), 'probe_-4.863': (1.96, 1.61, 2.39)}
# the model inputs' concentrations (log10) of the three odors: the lab's input (megamat pa and ms at -3, validation2 paa
# at -4), the KC panel's (paa at -2.954), B1's added stimuli and Hallem's series
MODEL_CONCS: Dict[str, Tuple[float, ...]] = {'pa': (-3.0, -4.0, -6.0, -8.0, -2.0), 'ms': (-3.0, -2.0, -4.0, -6.0, -8.0),
                                             'paa': (-4.0, -2.954, -5.0, -2.0)}
# the panel's inputs of the three odors (odor -> log10 concentration of the lab's input)
PANEL_INPUT_CONC = {'pa': -3.0, 'ms': -3.0, 'paa': -4.0}


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(path.read_text()) if path.exists() else None


_PROCESS_COMMITS: Dict[str, object] = {}


def process_commits(strict: bool = False) -> Dict[str, object]:
    """The commits of the code this process loaded, read once (`main` reads them before any work), so that a run's rows
    carry them even if the checkout moves while it runs; `strict` refuses a dirty checkout of this repository or of
    pyolfsysm (part B's runs)."""
    if not _PROCESS_COMMITS:
        _PROCESS_COMMITS.update({k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    if strict and any(_PROCESS_COMMITS.get(k) for k in ('queue_dirty', 'pyolfsysm_dirty')):
        raise RuntimeError(f'part B runs from clean checkouts only: {_PROCESS_COMMITS}')
    return dict(_PROCESS_COMMITS)


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                                         **extra))
    _write_json(path, d)


# ---------------------------------------------------------------------------
# the flies: responses, concentrations, ratios

CONC_TOL = 1e-6                        # two log10 concentrations closer than this are the same


def reading_label(reading: Tuple[str, bool]) -> str:
    return f"{reading[0]}{'_blank' if reading[1] else ''}"


def recording_info() -> pd.DataFrame:
    """recordings.csv indexed by recording name, with 'fly' (E79's `<date>_fly<n>`) and 'pulse_s'."""
    from experiments import E79_apl_calyx_constraints as E79
    r = pd.read_csv(E79.DATA / 'recordings.csv')
    r.index = [f.split('_moco')[0] for f in r['file']]
    r['fly'] = [f'{d}_fly{f}' for d, f in zip(r['date'], r['fly'])]
    r['pulse_s'] = r['odor_pulse_s'].astype(float)
    return r


def load_recordings():
    """E79's `load_recordings` (the 12 recordings, their stacks, frames and calyx masks)."""
    from experiments import E79_apl_calyx_constraints as E79
    return E79.load_recordings()


def stimulus_responses(recs, pulses: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """Per recording, stimulus, baseline (BASELINES) and blank (subtracted or not): the calyx means of the two
    presentations ('rep1', 'rep2', in frame order) and of their mean ('response'; E79's `stimulus_maps`, doubtful
    presentations left out), the recording's blank calyx mean ('blank'), its fly, side, pulse (`pulses`, default
    recordings.csv's) and the stimulus's block."""
    from experiments import E79_apl_calyx_constraints as E79
    if pulses is None:
        pulses = recording_info()['pulse_s'].to_dict()
    rows = []
    for rec in recs:
        fr, st, mask = rec.frames, rec.stack, rec.mask
        blank = float(E79.position_matched_blank(st, fr)[mask].mean())
        ok = fr[(fr['odor'] != 'pfo') & ~fr['doubtful']]
        block = {(o, float(c)): int(g['block'].min()) for (o, c), g in ok.groupby(['odor', 'log10_conc'])}
        for bl in BASELINES:
            for (o, c), (m, reps) in E79.stimulus_maps(st, fr, bl).items():
                if len(reps) != 2:
                    raise RuntimeError(f'{rec.name} {o} {c:g}: {len(reps)} presentations')
                vals = dict(response=float(m[mask].mean()), rep1=float(reps[0][mask].mean()), rep2=float(reps[1][mask].mean()))
                for sub in (False, True):
                    rows.append(dict(recording=rec.name, fly=rec.fly, side=rec.side, pulse_s=float(pulses[rec.name]),
                                     odor=o, conc=float(c), stimulus=E79.stim_label((o, c)), block=block[(o, float(c))],
                                     baseline=bl, blank_subtracted=sub, blank=blank,
                                     **{k: v - blank if sub else v for k, v in vals.items()}))
    return pd.DataFrame(rows)


def _reading(resp: pd.DataFrame, reading: Tuple[str, bool]) -> pd.DataFrame:
    return resp[(resp['baseline'] == reading[0]) & (resp['blank_subtracted'].astype(bool) == bool(reading[1]))]


def blank_identity(resp: pd.DataFrame) -> float:
    """The largest |first-in-block with the blank subtracted - order-matched without it| over recordings and stimuli
    (zero up to rounding: the two are the same reading)."""
    a = _reading(resp, ('first_in_block', True)).set_index(['recording', 'stimulus'])['response']
    b = _reading(resp, ('order_matched', False)).set_index(['recording', 'stimulus'])['response']
    if not a.index.sort_values().equals(b.index.sort_values()):
        raise RuntimeError('the two readings have different stimuli')
    return float((a - b.reindex(a.index)).abs().max())


def value_at(concs: Sequence[float], values: Sequence[float], c: float,
             design_range: Optional[Tuple[float, float]] = None,
             design_concs: Optional[Sequence[float]] = None) -> Tuple[float, str]:
    """A recording's response at log10 concentration `c` from its imaged concentrations `concs` and responses `values`:
    ('imaged') if imaged; ('interpolated') linear on log10 concentration between its two imaged neighbours inside its
    imaged range; ('substituted') the nearest imaged end when `c` is outside the series' `design_range` (default the
    recording's own range) and within MATCH_LOG_UNITS of that end, inclusive; otherwise NaN ('outside', or 'missing'
    when `c` is inside the design range but outside the recording's own, or, with `design_concs` (the series' imaged
    concentrations, part D), when a design concentration between the two neighbours was not imaged in this recording: a
    stimulus whose block is left out was not imaged)."""
    x = np.asarray(concs, dtype=float)
    y = np.asarray(values, dtype=float)
    if x.size == 0:
        return np.nan, 'missing'
    o = np.argsort(x)
    x, y = x[o], y[o]
    if np.any(np.diff(x) < CONC_TOL):
        raise ValueError(f'repeated concentrations {list(x)}')
    c = float(c)
    hit = np.flatnonzero(np.abs(x - c) <= CONC_TOL)
    if hit.size:
        return float(y[hit[0]]), 'imaged'
    lo, hi = float(x[0]), float(x[-1])
    if lo < c < hi:
        if design_concs is not None:
            a, b = float(x[x < c].max()), float(x[x > c].min())
            if any(a + CONC_TOL < float(d) < b - CONC_TOL for d in design_concs):
                return np.nan, 'missing'
        return float(np.interp(c, x, y)), 'interpolated'
    dlo, dhi = design_range if design_range is not None else (lo, hi)
    if dlo - CONC_TOL <= c <= dhi + CONC_TOL:
        return np.nan, 'missing'
    m = concentration_match(c, [lo, hi])
    if m is None:
        return np.nan, 'outside'
    return float(y[0] if m == lo else y[-1]), 'substituted'


def concentration_match(c: float, imaged: Sequence[float], max_units: float = MATCH_LOG_UNITS) -> Optional[float]:
    """The imaged concentration nearest to `c` if it is within `max_units` log10 units (inclusive), else None."""
    x = np.asarray(imaged, dtype=float)
    if x.size == 0:
        return None
    d = np.abs(x - float(c))
    i = int(np.argmin(d))
    return float(x[i]) if d[i] <= max_units + CONC_TOL else None


def fly_ratio(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The flies' ratio from per-recording log ratios (columns 'fly', 'log_ratio'; natural log): per fly the mean on
    the log scale, the geometric mean over flies, its 90 % interval from `n_boot` bootstrap draws of flies (E79's
    `fly_bootstrap`), the between-fly SD of the log2 ratio and the counts."""
    from experiments import E79_apl_calyx_constraints as E79
    d = per_recording.dropna(subset=['log_ratio'])
    per_fly = d.groupby('fly')['log_ratio'].mean()
    out = dict(n_flies=int(len(per_fly)), n_recordings=int(len(d)), ratio=np.nan, lo90=np.nan, hi90=np.nan,
               sd_log2_flies=np.nan, n_flies_above_1=int((per_fly > 0).sum()))
    if len(per_fly) == 0:
        return out
    mean, lo, hi = E79.fly_bootstrap(per_fly.to_numpy(), n_boot=n_boot, seed=seed)
    out.update(ratio=float(np.exp(mean[0])), lo90=float(np.exp(lo[0])), hi90=float(np.exp(hi[0])))
    if len(per_fly) > 1:
        out['sd_log2_flies'] = float(per_fly.std(ddof=1) / np.log(2))
    return out


def design_ranges(resp: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Each odor's concentration range over all recordings (the series' design range)."""
    g = resp.groupby('odor')['conc']
    return {o: (float(g.min()[o]), float(g.max()[o])) for o in g.groups}


def subset_recordings(resp: pd.DataFrame, subset: str) -> List[str]:
    """The recordings of a subset (SUBSETS): 'all'; 'pulse_1s' / 'pulse_0.5s'; 'without_large_blank_and_3FRfun1' (the
    recordings whose blank's calyx mean exceeds LARGE_BLANK, two, and NOT_IN_18UM_LIST left out); 'all_recording_ranges'
    (all, with the one-log rule applied to each recording's own range: `flies_ratios`)."""
    r = resp.drop_duplicates('recording').set_index('recording')
    if subset in ('all', 'all_recording_ranges'):
        keep = r.index
    elif subset == 'pulse_1s':
        keep = r.index[np.isclose(r['pulse_s'], 1.0)]
    elif subset == 'pulse_0.5s':
        keep = r.index[np.isclose(r['pulse_s'], 0.5)]
    elif subset == 'without_large_blank_and_3FRfun1':
        large = r.index[r['blank'] > LARGE_BLANK]
        if len(large) != 2:
            raise RuntimeError(f'expected two recordings with a large blank, found {list(large)}')
        keep = [x for x in r.index if x not in set(large) | set(NOT_IN_18UM_LIST)]
    else:
        raise ValueError(f'unknown subset {subset!r}')
    return sorted(keep)


def _stim(s: Stim) -> str:
    return f'{s[0]} {s[1]:g}'


def flies_ratios(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING, subset: str = 'all',
                 ratios=RATIOS, ranges: Optional[Dict[str, Tuple[float, float]]] = None,
                 recordings: Optional[Sequence[str]] = None,
                 design: Optional[Dict[str, Sequence[float]]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: every ratio of `ratios` under one reading and subset of the recordings: (one row per ratio, per recording
    rows with each response, how it was read and whether the recording entered). The one-log rule stands in a
    recording's end for a concentration outside the series' design range (`value_at`), or, in subset
    'all_recording_ranges', outside the recording's own range (ms at 10^-3 then reads ms at 10^-3.431 in the two
    recordings without ms at 10^-2). `ranges` gives the design ranges, `design` the design's concentrations per odor
    (part D passes the 12's: their ranges, and a reading that needs a design concentration the recording did not image is
    missing, `value_at`), `recordings` an explicit list of recordings in place of `subset`'s."""
    if design is not None and ranges is None:
        ranges = {o: (min(v), max(v)) for o, v in design.items()}
    dr = (design_ranges(resp) if ranges is None else dict(ranges)) if subset != 'all_recording_ranges' else {}
    keep = set(subset_recordings(resp, subset)) if recordings is None else set(recordings)
    r = _reading(resp, reading)
    r = r[r['recording'].isin(keep)]
    rec_rows, rows = [], []
    for rid, num, den, group in ratios:
        per = []
        for rec, g in r.groupby('recording'):
            row = dict(ratio_id=rid, reading=reading_label(reading), subset=subset, recording=rec, fly=g['fly'].iloc[0],
                       pulse_s=float(g['pulse_s'].iloc[0]))
            for side, (o, c) in (('num', num), ('den', den)):
                go = g[g['odor'] == o]
                v, how = value_at(go['conc'], go['response'], c, dr.get(o), None if design is None else design.get(o))
                row[f'{side}_value'], row[f'{side}_how'] = v, how
            row['entered'] = bool(np.isfinite(row['num_value']) and np.isfinite(row['den_value'])
                                  and row['num_value'] > MIN_RESPONSE and row['den_value'] > MIN_RESPONSE)
            row['log_ratio'] = float(np.log(row['num_value'] / row['den_value'])) if row['entered'] else np.nan
            per.append(row)
        pr = pd.DataFrame(per)
        rec_rows.append(pr)
        fr = fly_ratio(pr[pr['entered']])
        probe = PROBES.get(rid, (np.nan,) * 3)
        rows.append(dict(ratio_id=rid, group=group, numerator=_stim(num), denominator=_stim(den),
                         reading=reading_label(reading), subset=subset, **fr,
                         n_recordings_in_subset=int(len(pr)), n_recordings_left_out=int((~pr['entered']).sum()),
                         num_how=';'.join(sorted(set(pr['num_how']))), den_how=';'.join(sorted(set(pr['den_how']))),
                         probe_ratio=probe[0], probe_lo90=probe[1], probe_hi90=probe[2]))
    return pd.DataFrame(rows), pd.concat(rec_rows, ignore_index=True)


def presentation_se(g: pd.DataFrame) -> float:
    """A recording's typical standard error of a stimulus's response, the median over its stimuli of |rep1 - rep2| / 2."""
    return float(np.median(np.abs(g['rep1'] - g['rep2']) / 2.0))


def dose_response(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A2: per recording and odor the fold change of the response per log10 unit (exp of the slope of ln response on
    log10 concentration over the concentrations above MIN_RESPONSE, two or more needed), the linear slope (dF/F per log10
    unit, all concentrations), the lowest concentration above the blank (the lowest whose response, the blank
    subtracted, exceeds both MIN_RESPONSE and twice the recording's `presentation_se`); and the values at MODEL_CONCS
    (`value_at`)."""
    from scipy.stats import spearmanr
    dr = design_ranges(resp)
    r = _reading(resp, reading)
    rows, at_rows = [], []
    for rec, g in r.groupby('recording'):
        se = presentation_se(g)
        for o, go in g.groupby('odor'):
            go = go.sort_values('conc')
            x, y = go['conc'].to_numpy(float), go['response'].to_numpy(float)
            above = y > MIN_RESPONSE
            fold = float(np.exp(np.polyfit(x[above], np.log(y[above]), 1)[0])) if above.sum() >= 2 else np.nan
            thr = max(MIN_RESPONSE, 2.0 * se)
            over = np.flatnonzero(y > thr)
            rows.append(dict(recording=rec, fly=go['fly'].iloc[0], pulse_s=float(go['pulse_s'].iloc[0]), odor=o,
                             reading=reading_label(reading), n_concs=int(len(x)), n_above=int(above.sum()),
                             fold_per_log_unit=fold, slope_dff_per_log_unit=float(np.polyfit(x, y, 1)[0]),
                             spearman_conc=float(spearmanr(x, y)[0]), presentation_se=se, threshold=thr,
                             lowest_above_blank=float(x[over[0]]) if over.size else np.nan,
                             all_above_from_lowest=bool(over.size and (y[over[0]:] > thr).all()),
                             response_min=float(y.min()), response_max=float(y.max())))
            for c in MODEL_CONCS.get(o, ()):
                v, how = value_at(x, y, c, dr[o])
                at_rows.append(dict(recording=rec, fly=go['fly'].iloc[0], pulse_s=float(go['pulse_s'].iloc[0]), odor=o,
                                    conc=c, reading=reading_label(reading), value=v, how=how))
    return pd.DataFrame(rows), pd.DataFrame(at_rows)


def fly_geomean(per_recording: pd.DataFrame, col: str) -> dict:
    """`fly_ratio` of a positive per-recording quantity (its log averaged within flies)."""
    d = per_recording[['fly', col]].rename(columns={col: 'v'})
    d = d[np.isfinite(d['v']) & (d['v'] > 0)]
    return fly_ratio(pd.DataFrame(dict(fly=d['fly'], log_ratio=np.log(d['v']))))


# ---------------------------------------------------------------------------
# the model: APL's activity, the indicator, the comparison

E90_EXPERIMENT = 'E90_model_stages'


def e90_dir(tag: str) -> Path:
    """The folder holding E90's stage run `tag` (`E90_model_stages.load_stage`'s order: the local outputs first, then
    the model-output archive)."""
    for d in (H.HERE / 'outputs' / E90_EXPERIMENT / 'B', H.shared_dir() / 'data' / 'model_outputs' / E90_EXPERIMENT / 'B'):
        if (d / f'{tag}_apl.csv').exists():
            return d
    raise FileNotFoundError(f"E90's stage run {tag} not found")


def input_totals(orn: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per odor the lab input's total over its glomeruli: the sum of the ORN deltas ('input_total') and of their
    positive parts ('input_total_pos')."""
    orn = H.apl16_orn_deltas() if orn is None else orn
    t = pd.DataFrame(dict(input_total=orn.sum(axis=0), input_total_pos=orn.clip(lower=0).sum(axis=0)))
    t.index = [H.odor_abbrev(c) for c in orn.columns]
    return t


def family_apl(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """A4: per wiring and odor, from E90's stage runs (APL on): APL's window mean (a), pre-odor mean (a0), evoked
    activity (a - a0), the PN stage's evoked total (sum over glomeruli of the window minus the pre-odor means), the
    input's total (sum over the 36 glomeruli of the lab's ORN deltas) and the KCs' response rates with APL on and
    acutely off (two or more spikes, the runs' responses)."""
    tot = input_totals()
    rows = []
    for w in wirings:
        d = e90_dir(f'{w}_on')
        apl = pd.read_csv(d / f'{w}_on_apl.csv', index_col=0)
        pn = pd.read_parquet(d / f'{w}_on_pn_win.parquet') - pd.read_parquet(d / f'{w}_on_pn_pre.parquet')
        rate = {c: H.abbrev_columns(pd.read_parquet(e90_dir(f'{w}_{cc}') / f'{w}_{cc}_responses.parquet').astype(float)).mean()
                for c, cc in (('on', 'on'), ('off', 'acute_off'))}
        for o in apl.index:
            rows.append(dict(wiring=w, odor=str(o), a=float(apl.loc[o, 'win']), a0=float(apl.loc[o, 'pre']),
                             evoked=float(apl.loc[o, 'win'] - apl.loc[o, 'pre']), pn_total=float(pn[str(o)].sum()),
                             input_total=float(tot.loc[o, 'input_total']), input_total_pos=float(tot.loc[o, 'input_total_pos']),
                             kc_rate_on=float(rate['on'][o]), kc_rate_off=float(rate['off'][o]), source=str(d)))
    return pd.DataFrame(rows)


def profile_summary(apl: pd.DataFrame) -> pd.DataFrame:
    """A4's profile per wiring: APL's mean window activity, its CV over the odors (sample SD, as E95's A4), the pre-odor
    share, the evoked activity's CV and range, ms's evoked activity over the odors' mean, and the evoked activity's
    Spearman ρ and Pearson r with the input's total, the PN stage's total and the KCs' response rates."""
    from scipy.stats import pearsonr, spearmanr
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        e = g['evoked']
        row = dict(wiring=w, a_mean=float(g['a'].mean()), a_cv=float(g['a'].std() / g['a'].mean()),
                   a0=float(g['a0'].mean()), a0_share=float(g['a0'].mean() / g['a'].mean()), evoked_mean=float(e.mean()),
                   evoked_cv=float(e.std() / e.mean()), evoked_min=float(e.min()), evoked_min_odor=str(e.idxmin()),
                   evoked_max=float(e.max()), evoked_max_odor=str(e.idxmax()), ms_over_mean=float(e['ms'] / e.mean()),
                   ms_rank_from_lowest=int(e.rank()['ms']))
        for col in ('input_total', 'input_total_pos', 'pn_total', 'kc_rate_on', 'kc_rate_off'):
            row[f'spearman_evoked_{col}'] = float(spearmanr(e, g[col])[0])
            row[f'pearson_evoked_{col}'] = float(pearsonr(e, g[col])[0])
        row['pn_total_cv'] = float(g['pn_total'].std(ddof=0) / g['pn_total'].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def indicator_ratio(a1, a2, a0, n: float, form: str):
    """The model's ratio of stimulus 1 over stimulus 2 seen through an indicator of exponent `n`: 'baseline',
    (a1^n - a0^n) / (a2^n - a0^n); 'evoked', ((a1 - a0) / (a2 - a0))^n. Equal to the evoked ratio at n 1."""
    a1, a2, a0 = (np.asarray(x, dtype=float) for x in (a1, a2, a0))
    if form == 'baseline':
        return (a1 ** n - a0 ** n) / (a2 ** n - a0 ** n)
    if form == 'evoked':
        return ((a1 - a0) / (a2 - a0)) ** n
    raise ValueError(f'unknown form {form!r}')


def model_ratio_set(a1: float, a2: float, a0: float, ns: Sequence[float] = INDICATOR_NS,
                    forms: Sequence[str] = FORMS) -> pd.DataFrame:
    """One wiring's ratios over `ns` and `forms` (columns n, form, ratio)."""
    return pd.DataFrame([dict(n=float(n), form=f, ratio=float(indicator_ratio(a1, a2, a0, n, f))) for f in forms for n in ns])


def classify(ratios: Sequence[float], lo: float, hi: float) -> str:
    """One wiring against the flies' interval [lo, hi]: 'inside' if any ratio lies in it; 'flatter' if all are nearer 1
    than the interval or on the other side of 1; 'steeper' if all are beyond the interval's far end (from 1); 'mixed'
    otherwise (an interval that holds 1 has no 'flatter': nothing is nearer 1 than it)."""
    r = np.asarray(ratios, dtype=float)
    if not np.all(np.isfinite(r) & (r > 0)):
        raise ValueError(f'ratios must be positive and finite: {ratios}')
    x = np.log(r)
    L, U = np.log(float(lo)), np.log(float(hi))
    if not L <= U:
        raise ValueError(f'empty interval [{lo}, {hi}]')
    if np.any((x >= L) & (x <= U)):
        return 'inside'
    if L > 0:
        return 'flatter' if np.all(x < L) else 'steeper' if np.all(x > U) else 'mixed'
    if U < 0:
        return 'flatter' if np.all(x > U) else 'steeper' if np.all(x < L) else 'mixed'
    return 'steeper' if (np.all(x > U) or np.all(x < L)) else 'mixed'


def comparison_verdict(classes: Sequence[str], min_wirings: int = VERDICT_MIN_WIRINGS) -> str:
    """'inside', 'flatter' or 'steeper' if at least `min_wirings` of the wirings are; otherwise 'not decided'."""
    c = list(classes)
    for v in ('inside', 'flatter', 'steeper'):
        if sum(x == v for x in c) >= min_wirings:
            return v
    return 'not decided'


def hull(intervals: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    """The smallest interval holding all of `intervals` (V1: the two readings of paa's concentration)."""
    iv = [(float(a), float(b)) for a, b in intervals]
    if not iv or any(not (np.isfinite(a) and np.isfinite(b) and a <= b) for a, b in iv):
        raise ValueError(f'bad intervals {intervals}')
    return min(a for a, _ in iv), max(b for _, b in iv)


def v1_verdict(per_estimate: Dict[str, str]) -> str:
    """V1 from its comparison verdict under each estimate of ms's input: 'inside' if inside under any; 'flatter' or
    'steeper' only if that under all; otherwise 'not decided'."""
    v = list(per_estimate.values())
    if not v:
        raise ValueError('no estimate of ms\'s input')
    if 'inside' in v:
        return 'inside'
    for d in ('flatter', 'steeper'):
        if all(x == d for x in v):
            return d
    return 'not decided'


def wiring_classes(apl: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring the model's ratio of odor `num` over odor `den` (the panel's runs) through the indicator family, its
    range, the linear reading (n 1, evoked) and its class against [lo, hi]."""
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        a0 = float(g['a0'].iloc[0])
        if not np.allclose(g['a0'], a0):
            raise RuntimeError(f'{w}: the pre-odor level differs between odors')
        rs = model_ratio_set(g.loc[num, 'a'], g.loc[den, 'a'], a0)
        rows.append(dict(wiring=w, ratio_linear=float(indicator_ratio(g.loc[num, 'a'], g.loc[den, 'a'], a0, 1.0, 'evoked')),
                         ratio_min=float(rs['ratio'].min()), ratio_max=float(rs['ratio'].max()),
                         ratio_n2p25_baseline=float(rs[(rs.n == 2.25) & (rs.form == 'baseline')]['ratio'].iloc[0]),
                         ratio_n2p25_evoked=float(rs[(rs.n == 2.25) & (rs.form == 'evoked')]['ratio'].iloc[0]),
                         cls=classify(rs['ratio'], lo, hi), cls_linear=classify([rs['ratio'].iloc[0]], lo, hi)))
    return pd.DataFrame(rows)


def kc_rates(criterion: str = 'z1') -> pd.DataFrame:
    """The KC data's per-odor response rates (columns 'ctrl' and 'tnt', rows the 16 odors and 'pfo' where available):
    'z1' the mean over recordings of the fraction above Yang's 1 z; 'noise' the registry's rates in baseline-noise
    units (noise_ctrl_frac_<odor>, noise_tnt_frac_<odor>; the bands' criterion)."""
    if criterion == 'z1':
        return pd.DataFrame({k: H.yang_response_fractions(g, thresh=1.0).mean() for k, g in (('ctrl', 'TNT3c'), ('tnt', 'TNT'))})
    if criterion == 'noise':
        return pd.DataFrame({k: pd.Series({o: H.data_target(f'noise_{k}_frac_{o}') for o in H.APL16_ODORS})
                             for k in ('ctrl', 'tnt')})
    raise ValueError(f'unknown criterion {criterion!r}')


ACROSS_ODORS = (('paa_over_ms', 'paa', 'ms'), ('pa_over_ms', 'pa', 'ms'), ('pa_over_paa', 'pa', 'paa'))
E37_TABLE = ('E37_spike_only_apl', 'apl_activity.csv')


def alternative_apls(apl: pd.DataFrame) -> pd.DataFrame:
    """A5: the ratios of the across-odor comparisons (ACROSS_ODORS, the panel's concentrations) by other APLs: the
    family's (evoked, and the window mean; median and range over the wirings of `apl`); in proportion to the KCs'
    response rates in each genotype (Yang's criterion of 1 z, and the noise units of the bands beside); in proportion to
    the PN dendrites' totals and to the pebbled ORNs' totals (the megamat odors, the 28 glomeruli they share with the
    model; E90's `al_tables` / `al_statistics`; paa is not megamat); E37's saved table (window mean and peak of every
    variant, no pre-odor level)."""
    rows = []

    def add(source: str, kind: str, values: Dict[str, float], note: str = '', lo=None, hi=None):
        for rid, a, b in ACROSS_ODORS:
            va, vb = values.get(a, np.nan), values.get(b, np.nan)
            r = va / vb if np.isfinite(va) and np.isfinite(vb) and vb > 0 else np.nan
            rows.append(dict(source=source, kind=kind, ratio_id=rid, ratio=float(r),
                             ratio_min=float(lo[rid]) if lo is not None else np.nan,
                             ratio_max=float(hi[rid]) if hi is not None else np.nan, note=note))

    for col, kind in (('evoked', 'family APL, evoked (a - a0), median over wirings'),
                      ('a', 'family APL, window mean (a), median over wirings')):
        per = {}
        for w, g in apl.groupby('wiring', sort=False):
            s = g.set_index('odor')[col]
            per[w] = {rid: s[a] / s[b] for rid, a, b in ACROSS_ODORS}
        per = pd.DataFrame(per).T
        med = per.median()
        rows += [dict(source='family', kind=kind, ratio_id=rid, ratio=float(med[rid]), ratio_min=float(per[rid].min()),
                      ratio_max=float(per[rid].max()), note='8 wirings, E90 stage runs') for rid, _, _ in ACROSS_ODORS]
    for crit in ('z1', 'noise'):
        k = kc_rates(crit)
        for g, lab in (('ctrl', 'control flies'), ('tnt', 'APL output blocked (TNT)')):
            add(f'kc_rate_{crit}', f'in proportion to the KCs\' response rates, {lab} ({crit})', k[g].to_dict(),
                note='KC panel concentrations (pa, ms, paa at 10^-3)')
    from experiments import E90_model_stages as E90
    st = E90.al_statistics(**{k: v for k, v in E90.al_tables().items()})
    add('pn_dendrites', 'in proportion to the PN dendrites\' summed response (GH146, 28 glomeruli)', st['total_pn'].to_dict(),
        note='megamat odors only (paa not imaged)')
    add('orn_pebbled', 'in proportion to the ORNs\' summed response (pebbled, 28 glomeruli)', st['total_orn'].to_dict(),
        note='megamat odors only')
    tot = input_totals()
    add('input', "in proportion to the model's input total (36 glomeruli)", tot['input_total'].to_dict(),
        note='the lab input; paa at 10^-4')
    e37 = pd.read_csv(H.shared_dir() / 'results' / E37_TABLE[0] / E37_TABLE[1], index_col=0)
    for c in e37.columns:
        add(f'E37:{c}', 'E37 saved APL (time >= 0; mean or peak, no pre-odor level)', e37[c].to_dict(),
            note='older bases fitted at 5 % under the structure loss')
    return pd.DataFrame(rows)


def power_sds(a3_main: pd.DataFrame) -> dict:
    """A6's between-fly SDs of a log2 ratio: the median and the largest over A3's across-odor ratios (numerator and
    denominator different odors) whose two concentrations are both POWER_MIN_CONC or above (the panel's odors are at
    about 10^-3; ratios of weak responses spread more), under one reading and subset (`a3_main`), with the ratios used."""
    ok = []
    for rid, num, den, _ in RATIOS:
        if num[0] != den[0] and min(num[1], den[1]) >= POWER_MIN_CONC - CONC_TOL:
            ok.append(rid)
    x = a3_main.set_index('ratio_id').loc[ok, 'sd_log2_flies'].astype(float)
    return dict(median=float(x.median()), max=float(x.max()), ratios={k: float(v) for k, v in x.items()})


def flies_needed(delta_log2: float, sd_log2: float, z: Tuple[float, float] = POWER_Z) -> float:
    """A6: the number of flies whose geometric mean's 90 % interval excludes a value `delta_log2` (log2) away with 80 %
    power, at a between-fly SD `sd_log2` (log2): ceil(((z_a + z_b) sd / delta)^2), at least 2; inf when delta is 0."""
    d = abs(float(delta_log2))
    if d == 0:
        return float('inf')
    return float(max(2, int(np.ceil(((z[0] + z[1]) * float(sd_log2) / d) ** 2 - 1e-9))))


SUBSTITUTED_INPUT = ('paa', 'euc', '-aPine')    # the model's ORN input at another concentration than the KC panel's


def power_table(apl: pd.DataFrame, sd_log2: Dict[str, float]) -> pd.DataFrame:
    """A6: for every pair of panel odors (oriented so that the proportional APL's ratio is at least 1) the family's
    evoked ratio (median and range over the wirings), the proportional APL's (the KCs' response rates with APL output
    blocked (TNT), control flies beside; Yang's 1 z), their log2 difference (and that from the wiring nearest the
    proportional APL) and the flies needed at each between-fly SD of `sd_log2`; sorted by the size of the difference,
    so the pairs that need the fewest flies come first. The ratios are linear readings (n 1, evoked): through an
    indicator of n above 1 the family's stays near 1 and the proportional APL's moves further away, so the counts are
    upper bounds in that respect."""
    from itertools import combinations
    k = kc_rates('z1')
    ev = apl.pivot(index='odor', columns='wiring', values='evoked')
    rows = []
    for a, b in combinations(list(H.APL16_ODORS), 2):
        num, den = (a, b) if k.loc[a, 'tnt'] >= k.loc[b, 'tnt'] else (b, a)
        fam = ev.loc[num] / ev.loc[den]
        prop = float(k.loc[num, 'tnt'] / k.loc[den, 'tnt'])
        fmed = float(fam.median())
        d = float(np.log2(prop) - np.log2(fmed))
        row = dict(numerator=num, denominator=den, substituted_input=bool({num, den} & set(SUBSTITUTED_INPUT)), family_ratio=fmed, family_min=float(fam.min()), family_max=float(fam.max()),
                   proportional_tnt=prop, proportional_ctrl=float(k.loc[num, 'ctrl'] / k.loc[den, 'ctrl']),
                   delta_log2=d, delta_log2_nearest_wiring=float(np.min(np.abs(np.log2(prop) - np.log2(fam.to_numpy(float))))))
        for name, sd in sd_log2.items():
            row[f'flies_needed_sd_{name}'] = flies_needed(d, sd)
            row[f'flies_needed_nearest_wiring_sd_{name}'] = flies_needed(row['delta_log2_nearest_wiring'], sd)
        rows.append(row)
    df = pd.DataFrame(rows)
    return df.reindex(df['delta_log2'].abs().sort_values(ascending=False, kind='stable').index).reset_index(drop=True)


# ---------------------------------------------------------------------------
# part B: the extended input and its runs (phase 2)

B_ARMS = ('B1', 'B2', 'B3', 'B4')
CONDS = (('on', 'on'), ('off', 'acute_off'))    # (condition, its file tag; E90's)
SESSIONS = ('kiwi-control', 'megamat', 'validation2')    # the lab's three ORN sessions (E100's SESSIONS)
LOW_SESSION = 'kiwi-control'     # the session whose diagnostics carry B1's low-concentration stimuli (E49's: ma, ea)
PAA5_SESSION = 'validation2'     # paa at 10^-5: the diagnostics of the flies of paa at 10^-4 (the entry)
DIAG = 'glomeruli_diagnostics'
ADDED_PANEL = 'E98'
BLANK_LABEL = 'blank @ 0'
HALLEM_NAMES = {'pa': 'pentyl acetate', 'ms': 'methyl salicylate', 'paa': 'phenylacetaldehyde'}
HALLEM_STIMULI: Tuple[Stim, ...] = (('pa', -2.0), ('pa', -4.0), ('pa', -6.0), ('pa', -8.0), ('ms', -2.0), ('ms', -4.0),
                                    ('ms', -6.0), ('ms', -8.0), ('paa', -2.0))
APL_TOL = 1e-9                   # part 0: APL's means against E90's, and the blank's evoked activity (absolute)
SETTLED_PRE_WINDOW = (-0.25, 0.0)    # APL's pre-odor level after the KC layer's start-up (it starts at time_start, -0.5 s)
A0_DEFINITIONS = ('e90', 'settled', 'blank')   # part C's pre-odor level: E90's pre-odor mean, 'pre_settled', the blank's window


def session_label(odor: str, conc: float, session: str, panel: str) -> str:
    """A session column's B1 label, E100's `estimate_label` form: '<odor> [<session>/<panel>] @ <conc>' ('diag' for the
    diagnostics), one ' @ ' as fit_mb_model's odor-name parsing wants."""
    return f"{odor} [{session}/{'diag' if panel == DIAG else panel}] @ {float(conc):g}"


def hallem_label(odor: str, conc: float) -> str:
    return f'{odor} [hallem] @ {float(conc):g}'


def b1_specs() -> List[dict]:
    """B1's 52 stimuli in their fixed order: label, odor, log10 concentration, source ('panel', '<session>/<panel>',
    'hallem', 'blank'), role, and for session columns the session, panel and column. Roles: 'panel' (the 16 odors of the
    lab's input), 'V1' (the estimates of ms's input at 10^-3 besides the panel's), 'V2' (paa at 10^-5 from PAA5_SESSION),
    'beside' (the other estimates of paa at 10^-5; 2h, t2h, 2-but and 1-6ol at 10^-6 and 2h at 10^-5; ma at 10^-7 and ea
    at 10^-8, from LOW_SESSION), 'E100' (4h-ol at 10^-6 and farn at 10^-2 from all three sessions, and the other two
    sessions' 2h and 2-but at 10^-6 and ma at 10^-7: E100's model side, read by none of E98's parts), 'hallem' (V3 to V5,
    descriptive), 'blank'."""
    specs = [dict(label=o, odor=o, conc=float(H.ORN_INPUT_CONC[o]), source='panel', role='panel') for o in H.APL16_ODORS]

    def ses(odor, conc, session, panel, role):
        specs.append(dict(label=session_label(odor, conc, session, panel), odor=odor, conc=float(conc),
                          source=f"{session}/{'diag' if panel == DIAG else panel}", role=role, session=session, panel=panel,
                          column=f'{odor} @ {conc:g}'))

    ses('paa', -5, PAA5_SESSION, DIAG, 'V2')
    for sn in SESSIONS:
        if sn != PAA5_SESSION:
            ses('paa', -5, sn, DIAG, 'beside')
    for sn in SESSIONS:
        ses('ms', -3, sn, DIAG, 'V1')
    ses('ms', -3, 'kiwi-control', 'control', 'V1')
    for o in ('2h', 't2h', '2-but', '1-6ol'):
        ses(o, -6, LOW_SESSION, DIAG, 'beside')
    ses('2h', -5, 'kiwi-control', 'control', 'beside')
    ses('ma', -7, LOW_SESSION, DIAG, 'beside')
    ses('ea', -8, LOW_SESSION, DIAG, 'beside')
    for o, c in (('4h-ol', -6), ('farn', -2)):
        ses(o, c, LOW_SESSION, DIAG, 'E100')
        for sn in SESSIONS:
            if sn != LOW_SESSION:
                ses(o, c, sn, DIAG, 'E100')
    for o, c in (('2h', -6), ('2-but', -6), ('ma', -7)):
        for sn in SESSIONS:
            if sn != LOW_SESSION:
                ses(o, c, sn, DIAG, 'E100')
    specs += [dict(label=hallem_label(o, c), odor=o, conc=float(c), source='hallem', role='hallem') for o, c in HALLEM_STIMULI]
    specs.append(dict(label=BLANK_LABEL, odor='blank', conc=0.0, source='blank', role='blank'))
    labels = [x['label'] for x in specs]
    assert len(set(labels)) == len(labels) == 52, labels
    return specs


def session_inputs() -> Dict[str, pd.DataFrame]:
    """The lab's three ORN sessions' spike deltas (E100's `session_inputs`: olfsysm's
    pyolfsysm/data/internal/<session>_mean_est_spike_deltas.csv, glomerulus x (panel, 'odor @ conc'))."""
    from experiments import E100_apl_weak_odors as E100
    return E100.session_inputs(SESSIONS)


HALLEM_FILE = 'hc_data.csv'          # olfsysm's tabulation of Hallem & Carlson 2006 (pyolfsysm/data; tracked)
OLFDB_HALLEM_COPY = Path('outputs/_planning/atlas_2026-11-15/hallem_carlson_2006_olfdb.json')   # olfdb's, fetched 16:30
HALLEM_DM3_RECEPTOR = '47a'           # the file's header gives dm3 to both 33b and 47a; DM3 is Or47a (olfdb, E38's drosolf)


def hallem_glomerulus(name: str) -> str:
    """hc_data.csv's lower-case glomerulus name in the form drosolf and the model use ('da4m' -> 'DA4m', 'vc3' -> 'VC3')."""
    import re
    m = re.fullmatch(r'([a-z]+)(\d+)([a-z]*)', str(name).strip())
    if m is None:
        raise ValueError(f'unexpected glomerulus name {name!r}')
    return m.group(1).upper() + m.group(2) + m.group(3)


def hallem_hc_table(path: Optional[Path] = None) -> Tuple[pd.DataFrame, dict]:
    """olfsysm's Hallem & Carlson 2006 table (pyolfsysm/data/hc_data.csv): odorant rows (the class-11 concentration
    series labelled '<odorant> -<k>') x glomeruli (the first header row, put in drosolf's form; dm3 taken from
    HALLEM_DM3_RECEPTOR, the other dm3 column (Or33b) dropped as E38's drosolf table drops the multiglomerular one),
    spike-rate changes; and its record (path, sha256, the olfsysm commit)."""
    import hashlib
    path = H.pyolfsysm_data_dir() / HALLEM_FILE if path is None else Path(path)
    raw = pd.read_csv(path, header=None)
    gl, rec = list(raw.iloc[0, 2:]), list(raw.iloc[1, 2:])
    keep = [i for i, (g, r) in enumerate(zip(gl, rec)) if not (str(g).strip() == 'dm3' and str(r).strip() != HALLEM_DM3_RECEPTOR)]
    t = raw.iloc[2:, :]
    vals = t.iloc[:, 2:].iloc[:, keep].astype(float)
    vals.columns = [hallem_glomerulus(gl[i]) for i in keep]
    if vals.columns.duplicated().any():
        raise RuntimeError(f'{path}: repeated glomeruli {list(vals.columns[vals.columns.duplicated()])}')
    vals.index = [str(x).strip() for x in t.iloc[:, 1]]
    vals.index.name = 'odorant'
    vals.columns.name = 'glomerulus'
    info = dict(path=str(path), sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(), dm3_receptor=HALLEM_DM3_RECEPTOR,
                olfsysm_commit=H.code_commits().get('pyolfsysm_commit'))
    return vals, info


def hallem_rows(path: Optional[Path] = None, lower: bool = True) -> pd.DataFrame:
    """Hallem's spike deltas (glomerulus x (odor, log10 concentration)) for HALLEM_STIMULI, from olfsysm's hc_data.csv
    (`hallem_hc_table`): the 10^-2 rows of pa and ms from its concentration series ('pentyl acetate -2', ...), paa's from
    its 10^-2 table ('phenylacetaldehyde'), the lower rows from the series (`lower` False: the 10^-2 rows only). Part 0
    checks that the 10^-2 rows equal E38's (drosolf's table, `hallem_identity`)."""
    t, _ = hallem_hc_table(path)
    cols = {}
    for o, c in HALLEM_STIMULI:
        if not lower and abs(c + 2.0) > CONC_TOL:
            continue
        name = HALLEM_NAMES[o]
        row = f'{name} {int(round(c))}' if f'{name} {int(round(c))}' in t.index else (name if abs(c + 2.0) < CONC_TOL else None)
        if row is None or row not in t.index:
            continue
        cols[(o, c)] = t.loc[row]
    out = pd.DataFrame(cols)
    out.index.name = 'glomerulus'
    return out


def hallem_identity() -> dict:
    """Part 0's identity of Hallem's arm: through this path (`hallem_rows`, olfsysm's hc_data.csv) Hallem's pa, ms and paa
    at 10^-2 equal E38's inputs (drosolf's 10^-2 table as E38's `hallem_deltas` reads it; pa through the same call) on
    every glomerulus the input shares with it, and the file's concentration series at 10^-2 equals its 10^-2 table."""
    from drosolf import orns
    from experiments import E38_hallem_input as E38
    gl = list(H.apl16_orn_deltas().index)
    ours = hallem_rows(lower=False)
    e38 = E38.hallem_deltas(gl)
    d = orns.orns(add_sfr=False, drop_sfr=False, columns='glomerulus')
    d = d.loc[:, [g for g in d.columns if '+' not in g]]
    shared = [g for g in gl if g in d.columns]
    ref = {'pa': d.loc[HALLEM_NAMES['pa'], shared].astype(float), 'ms': e38['ms'], 'paa': e38['paa']}
    t, info = hallem_hc_table()
    res = dict(file=info, shared_glomeruli=shared)
    for o in ('pa', 'ms', 'paa'):
        mine = ours[(o, -2.0)].reindex(shared)
        r = ref[o].reindex(shared).astype(float)
        res[o] = dict(max_abs_diff=float((mine - r).abs().max()), n=int(len(shared)), missing=int(mine.isna().sum()))
    for o in ('pa', 'ms'):
        res[f'{o}_series_vs_table'] = float((t.loc[f'{HALLEM_NAMES[o]} -2'] - t.loc[HALLEM_NAMES[o]]).abs().max())
    res['exact'] = bool(all(res[o]['max_abs_diff'] == 0 and res[o]['missing'] == 0 for o in ('pa', 'ms', 'paa'))
                        and res['pa_series_vs_table'] == 0 and res['ms_series_vs_table'] == 0)
    return res


def olfdb_hallem_check(path: Path) -> dict:
    """A one-time cross-check of `hallem_rows` against olfdb's reference dataset hallem-carlson-2006 (a JSON fetched from
    its API: rows glomeruli, cols '<odorant> @ <k>', float32 values_b64): the largest difference over Hallem's stimuli
    and the glomeruli both have."""
    import base64
    d = json.loads(Path(path).read_text())
    v = np.frombuffer(base64.b64decode(d['values_b64']), dtype='<f4').reshape(d['shape'])
    lab = lambda x: str(x['label'] if isinstance(x, dict) else x)
    ref = pd.DataFrame(v, index=[lab(r) for r in d['rows']], columns=[lab(c) for c in d['cols']])
    ref = ref.loc[:, ~ref.columns.duplicated()]
    ours = hallem_rows()
    out = {}
    for (o, c), col in ours.items():
        name = f'{HALLEM_NAMES[o]} @ {int(round(c))}'
        if name not in ref.columns:
            out[f'{o} {c:g}'] = 'not in olfdb'
            continue
        g = [x for x in col.index if x in ref.index]
        out[f'{o} {c:g}'] = dict(n=len(g), max_abs_diff=float((col.loc[g].astype(float) - ref.loc[g, name].astype(float)).abs().max()))
    return out


def b1_input(include_hallem_lower: bool = True, hallem_path: Optional[Path] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """B1's extended input: the lab's 16 odors (H.apl16_orn_deltas, unchanged and first), then the added stimuli of
    `b1_specs` (panel ADDED_PANEL, odor the stimulus label), on the input's 36 glomeruli; a session's glomeruli outside
    them dropped and the ones it lacks filled with 0; Hallem's stimuli on the glomeruli Hallem and the input share, 0
    elsewhere (E38's 'hallem_only'); the blank all 0. Returns (input, table): the table has one row per B1 stimulus
    (`b1_specs`' fields, 'included', why not, the glomeruli filled with 0 and dropped, the input total and its positive
    part) and the input vector of every included stimulus (columns orn_<glomerulus>). Hallem's rows come from olfsysm's
    hc_data.csv (`hallem_rows`); without `include_hallem_lower` the rows below 10^-2 are left out and listed ('included'
    False, with the reason)."""
    panel = H.apl16_orn_deltas()
    gl = list(panel.index)
    sess = session_inputs()
    hal = hallem_rows(hallem_path, lower=include_hallem_lower)
    cols, rows = {}, []
    for sp in b1_specs():
        row = dict(sp, included=True, reason='', glomeruli_filled_0='', glomeruli_dropped='')
        if sp['source'] == 'panel':
            v = panel[[c for c in panel.columns if H.odor_abbrev(c) == sp['odor']][0]].astype(float)
        elif sp['source'] == 'blank':
            v = pd.Series(0.0, index=gl)
        elif sp['source'] == 'hallem':
            key = (sp['odor'], sp['conc'])
            if key not in hal.columns:
                row.update(included=False, reason=("left out (include_hallem_lower False)" if not include_hallem_lower
                                                   else f"not in {HALLEM_FILE}"))
                rows.append(row)
                continue
            h = hal[key]
            shared = [g for g in gl if g in h.index]
            v = pd.Series(0.0, index=gl)
            v.loc[shared] = h.loc[shared].to_numpy(float)
            row.update(glomeruli_filled_0=';'.join(g for g in gl if g not in shared),
                       glomeruli_dropped=';'.join(sorted(set(h.index) - set(gl))))
        else:
            src = sess[sp['session']]
            col = (sp['panel'], sp['column'])
            if col not in src.columns:
                raise RuntimeError(f"{sp['label']}: {col} not in the {sp['session']} session's table")
            raw = src[col]
            v = raw.reindex(gl)
            row.update(glomeruli_filled_0=';'.join(g for g in gl if g not in raw.index or pd.isna(raw.get(g))),
                       glomeruli_dropped=';'.join(sorted(set(raw.index) - set(gl))))
            v = v.fillna(0.0).astype(float)
        cols[sp['label']] = v.to_numpy(float)
        row.update(input_total=float(v.sum()), input_total_pos=float(v.clip(lower=0).sum()),
                   **{f'orn_{g}': float(x) for g, x in zip(gl, v.to_numpy(float))})
        rows.append(row)
    extra_labels = [k for k in cols if k not in H.APL16_ODORS]
    extra = pd.DataFrame({k: cols[k] for k in extra_labels}, index=panel.index)
    extra.columns = pd.MultiIndex.from_tuples([(ADDED_PANEL, k) for k in extra_labels], names=panel.columns.names)
    orn = pd.concat([panel.astype(float), extra], axis='columns')
    orn.index.name = panel.index.name
    for o, c in zip(H.APL16_ODORS, panel.columns):
        if not np.array_equal(orn[c].to_numpy(float), panel[c].to_numpy(float)):
            raise RuntimeError(f'{o}: the panel column changed')
    return orn, pd.DataFrame(rows)


def input_labels(orn: pd.DataFrame) -> List[str]:
    """B1 labels of an extended input's columns: the panel odors' abbreviations, the added stimuli's labels."""
    return [c[-1] if c[0] == ADDED_PANEL else H.odor_abbrev(c) for c in orn.columns]


class ExtendedRun:
    """One wiring at its fit (E83's `build_simulator`, breadth held at the lab input's, so the APL->KC weights are the
    fit's) on an extended input (default B1's), APL's activity saved (`mp.kc.save_inh_sims`); `run(cond)` simulates APL
    on ('on') or acutely off ('off': the APL->KC and APL->bouton scales 0, as E90's StageRun) at the fit's parameters
    and returns the responses and spike counts (columns the stimulus labels, in the input's order), APL's window and
    pre-odor means per stimulus ('apl', E90's MODEL_WINDOW and PRE_WINDOW) and the PN stage's window and pre-odor means
    ('pn_win', 'pn_pre', glomeruli x stimuli)."""

    def __init__(self, label: str, orn: Optional[pd.DataFrame] = None, params: Optional[Dict[str, float]] = None):
        from experiments import E83_cross_wiring_fit as E83
        P = E83.wiring_parts(label)
        self.label, self.W, self.saved = label, P['W'], P['saved']
        self.p = dict(P['p'] if params is None else params)
        self.orn = b1_input()[0] if orn is None else orn
        self.labels = input_labels(self.orn)
        self.sim, self.info = E83.build_simulator(label, orn=self.orn, breadth=P['breadth'], params=self.p)
        cells = self.info['cells']
        self.mp, self.rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        self.order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        self.mp.kc.save_inh_sims = True

    def params(self, cond: str) -> Dict[str, float]:
        return dict(self.p) if cond == 'on' else dict(self.p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)

    def run(self, cond: str) -> dict:
        r, s = self.sim(self.params(cond))
        gl = pd.Index(list(dict.fromkeys(self.W.columns.get_level_values('glomerulus'))), name='glomerulus')
        return dict(read_apl_pn(self.mp, self.rv, self.order, self.labels, gl), **relabel(r, s, self.labels))


def relabel(r: pd.DataFrame, s: pd.DataFrame, labels: Sequence[str]) -> dict:
    """The simulator's responses and spike counts with the stimulus labels as columns (the outputs' columns must be the
    input's, in its order)."""
    if r.shape[1] != len(labels) or s.shape[1] != len(labels):
        raise RuntimeError(f'{r.shape[1]} output columns for {len(labels)} stimuli')
    raw = [str(c[-1]) if isinstance(c, tuple) else str(c) for c in s.columns]
    out_lab = [x if x == lab else H.odor_abbrev(x) for x, lab in zip(raw, labels)]
    if out_lab != list(labels):
        raise RuntimeError(f"the outputs' columns are not the input's: {out_lab} / {list(labels)}")
    r, s = r.copy(), s.copy()
    r.columns = s.columns = pd.Index(list(labels), name='stimulus')
    return dict(responses=r, spike_counts=s)


def read_apl_pn(mp, rv, order: np.ndarray, labels: Sequence[str], glomeruli: pd.Index,
                pulse_window: Optional[Tuple[float, float]] = None) -> dict:
    """APL's activity and the PN stage of the last simulation (olfsysm's variables; the simulator's output column j is
    olfsysm's odor order[j]): 'apl' (index the stimulus labels; 'win' the mean over E90's MODEL_WINDOW, 0-0.5 s, 'pre'
    over E90's PRE_WINDOW, 'pre_settled' over SETTLED_PRE_WINDOW, and 'win_pulse' over `pulse_window` when given, e.g.
    0-1 s for a 1 s pulse), 'apl_trace' (APL's activity at every step from olfsysm's time_start, steps x stimuli) and
    'pn_win' / 'pn_pre' (glomeruli x stimuli, the PN rates' window and pre-odor means).

    olfsysm's KC layer (and so APL) starts from zero at time_start (-0.5 s) and reaches its spontaneous level about 0.2 s
    later, so E90's pre-odor mean over -0.5 to 0 s includes that start-up (on the real wiring 8.27 against a settled
    8.79); 'pre_settled' leaves it out."""
    from experiments import E90_model_stages as E90
    from pyolfsysm.mb_model import get_sim_var
    t = E90.step_times(mp, rv)
    w = E90._windows(t, E90.MODEL_WINDOW)
    wp = E90._windows(t, pulse_window)['win'] if pulse_window is not None else None
    ws = np.flatnonzero((t >= SETTLED_PRE_WINDOW[0] - 1e-9) & (t < SETTLED_PRE_WINDOW[1] - 1e-9))
    kept = np.flatnonzero(t >= float(mp.time_start) - 1e-9)
    inh = get_sim_var(rv.kc, 'inh_sims')
    pn = get_sim_var(rv.pn, 'pn_sims')
    apl, pw, pp, trace = {}, {}, {}, {}
    for j, lab in enumerate(labels):
        i = int(order[j])
        a = np.asarray(inh[i], dtype=float).ravel()
        apl[lab] = dict(win=float(a[w['win']].mean()), pre=float(a[w['pre']].mean()), pre_settled=float(a[ws].mean()))
        if wp is not None:
            apl[lab]['win_pulse'] = float(a[wp].mean())
        trace[lab] = a[kept]
        x = np.asarray(pn[i], dtype=float)
        pw[lab], pp[lab] = x[:, w['win']].mean(axis=1), x[:, w['pre']].mean(axis=1)
    if len(glomeruli) != len(next(iter(pw.values()))):
        raise RuntimeError("PN rows are not the wiring's glomeruli")
    apl = pd.DataFrame(apl).T
    apl.index.name = 'stimulus'
    tr = pd.DataFrame(trace, index=pd.Index(np.round(t[kept], 6), name='time_s'))
    return dict(apl=apl, apl_trace=tr, pn_win=pd.DataFrame(pw, index=glomeruli), pn_pre=pd.DataFrame(pp, index=glomeruli),
                n_steps=dict(win=int(len(w['win'])), pre=int(len(w['pre'])), dt=float(t[1] - t[0])))


def b_root(root: Optional[Path] = None) -> Path:
    """The folder of part B's runs: `root`/B (default outputs/E98_apl_response_strength/B)."""
    d = (Path(root) if root is not None else H.outputs_dir(EXPERIMENT)) / 'B'
    d.mkdir(parents=True, exist_ok=True)
    return d


def b1_roots() -> List[Path]:
    """Where B1's runs may be, searched in this order: the local outputs, the model-output archive, and hal's runs in
    the archive (from_hal/)."""
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [H.HERE / 'outputs' / EXPERIMENT / 'B' / 'B1', arc / 'B' / 'B1', arc / 'from_hal' / 'B' / 'B1']


B1_LAYOUT = """B1's output layout (what E100's `read_e98_b1` and this module's `load_b1` read). One folder per root
(`b1_roots`: outputs/E98_apl_response_strength/B/B1/, data/model_outputs/E98_apl_response_strength/B/B1/,
data/model_outputs/E98_apl_response_strength/from_hal/B/B1/), per wiring w and condition c ('on', 'acute_off'):
- <w>_<c>_apl.csv: index 'stimulus' (the B1 label), columns win and pre (APL's mean over 0-0.5 s and over the
  0.5 s before the odor, E90's windows) and pre_settled (over SETTLED_PRE_WINDOW, after the KC layer's start-up);
- <w>_<c>_apl_trace.parquet: APL's activity at every step from olfsysm's time_start (index time_s) x stimuli;
- <w>_<c>_responses.parquet, <w>_<c>_spike_counts.parquet: KCs x stimuli, columns the B1 labels;
- <w>_<c>_pn_win.parquet, <w>_<c>_pn_pre.parquet: glomeruli x stimuli, the PN stage's window and pre-odor means;
- <w>_identity.json: part 0's run identities of that wiring's run;
and once per folder b1_stimuli.csv (`b1_input`'s table: the labels, odors, concentrations, sources, roles, the input
vectors) and runs.csv (one row per wiring and condition: the absolute parameters, the commits, host, time).
Labels: a panel odor by its abbreviation (its input at the lab's concentration, H.ORN_INPUT_CONC); a session's column
'<odor> [<session>/<panel>] @ <conc>' (sessions 'kiwi-control', 'megamat', 'validation2'; panel 'diag' for the
diagnostics, else the panel's name), as E100's `estimate_label`; Hallem's '<odor> [hallem] @ <conc>'; the blank
'blank @ 0'. The labels and their order are fixed (`b1_specs`); a Hallem stimulus without a source is absent."""


def save_b1(d: Path, label: str, cond_tag: str, o: dict) -> None:
    for k in ('responses', 'spike_counts', 'pn_win', 'pn_pre', 'apl_trace'):
        if k not in o:
            continue
        x = o[k].copy()
        x.columns = [str(c) for c in x.columns]
        x.to_parquet(d / f'{label}_{cond_tag}_{k}.parquet')
    o['apl'].to_csv(d / f'{label}_{cond_tag}_apl.csv')


def load_b1_run(label: str, roots: Optional[Sequence[Path]] = None) -> Optional[Dict[str, dict]]:
    """One wiring's saved B1 run ({cond_tag: {'apl', 'responses', 'spike_counts', 'pn_win', 'pn_pre', 'dir'}}) from the
    first root holding both conditions, or None."""
    for d in (roots or b1_roots()):
        d = Path(d)
        if all((d / f'{label}_{ct}_apl.csv').exists() for _, ct in CONDS):
            out = {}
            for _, ct in CONDS:
                o = {k: pd.read_parquet(d / f'{label}_{ct}_{k}.parquet') for k in ('responses', 'spike_counts', 'pn_win', 'pn_pre')}
                o['apl'] = pd.read_csv(d / f'{label}_{ct}_apl.csv', index_col=0)
                o['dir'] = str(d)
                out[ct] = o
            return out
    return None


def load_b1(wirings: Sequence[str] = ALL_WIRINGS, roots: Optional[Sequence[Path]] = None) -> pd.DataFrame:
    """B1's saved runs as E100's `read_e98_b1` gives them: one row per wiring, condition ('on', 'acute_off') and stimulus,
    with APL's window mean ('win'), pre-odor mean ('pre'), evoked activity ('evoked'), the KCs' response rate
    ('response_rate', two or more spikes), the folder ('file') and the label parsed (b1_odor, b1_conc, b1_source; E100's
    `parse_b1_label`); attrs['missing'] the (wiring, condition) pairs not found. See B1_LAYOUT."""
    from experiments import E100_apl_weak_odors as E100
    rows, found = [], set()
    for w in wirings:
        run = load_b1_run(w, roots)
        if run is None:
            continue
        for ct, o in run.items():
            rate = o['responses'].astype(float).mean(axis=0)
            for stim, a in o['apl'].iterrows():
                lab = E100.parse_b1_label(stim)
                rows.append(dict(wiring=w, cond=ct, stimulus=str(stim), win=float(a['win']), pre=float(a['pre']),
                                 evoked=float(a['win'] - a['pre']), response_rate=float(rate.get(str(stim), np.nan)),
                                 file=o['dir'], b1_odor=lab['odor'], b1_conc=lab['conc'], b1_source=lab['source'],
                                 pre_settled=float(a['pre_settled']) if 'pre_settled' in a.index else np.nan))
            found.add((w, ct))
    if not rows:
        raise NotImplementedError("E98's B1 runs are not there yet (" + ', '.join(str(d) for d in (roots or b1_roots())) + ')')
    df = pd.DataFrame(rows)
    df.attrs['missing'] = sorted(set((w, ct) for w in wirings for _, ct in CONDS) - found)
    return df


E101_A3_TABLE = ('results', 'E101_apl_vs_kc_claws', 'a3_apl.csv')
E101_A3_PROVENANCE = ('results', 'E101_apl_vs_kc_claws', 'provenance_A3.json')
B1_COMPARE_COLS = ('win', 'pre', 'pre_settled', 'response_rate')


E101_A3_ARCHIVE = ('data', 'model_outputs', 'E101_apl_vs_kc_claws', 'A3')


def load_b1_e101(path: Optional[Path] = None) -> Tuple[pd.DataFrame, dict]:
    """B1's stimuli as E101's part A3 ran them on the eight wirings (hal, 2026-11-15; E98's `ExtendedRun` and `save_b1`,
    E98's B1 layout plus E101's stages): `load_b1`'s rows from its archived runs
    (data/model_outputs/E101_apl_vs_kc_claws/A3/), or from its table (results/E101_apl_vs_kc_claws/a3_apl.csv, E98's
    `load_b1` written out) when `path` names a CSV or the archive is absent; and its record (where, the sha256 of every
    file read, E101's A3 provenance: the commits of the runs and of the table)."""
    import hashlib
    prov = _read_json(H.shared_dir().joinpath(*E101_A3_PROVENANCE)) or {}
    rec = dict(runs=[dict(commits=r.get('commits'), host=r.get('host'), time=r.get('time'), runs_only=r.get('runs_only'))
                     for r in prov.get('runs', [])])
    arc = H.shared_dir().joinpath(*E101_A3_ARCHIVE)
    if path is None and arc.is_dir():
        b = load_b1(ALL_WIRINGS, roots=[arc])
        files = sorted(f for f in arc.iterdir() if f.name.endswith(('_apl.csv', '_responses.parquet')))
        rec.update(path=str(arc), sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
        b['file'] = 'E101:' + str(arc)
        return b, rec
    path = H.shared_dir().joinpath(*E101_A3_TABLE) if path is None else Path(path)
    b = pd.read_csv(path, float_precision='round_trip')
    rec.update(path=str(path), sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
    b['file'] = f'E101:{Path(path).name}'
    return b, rec


def load_b1_all(wirings: Sequence[str] = ALL_WIRINGS, e101: bool = True, roots: Optional[Sequence[Path]] = None,
                e101_path: Optional[Path] = None) -> pd.DataFrame:
    """B1's rows (`load_b1`'s columns) on `wirings`: E98's own runs where they exist (`b1_roots`), the other wirings from
    E101's A3 (`load_b1_e101`) when `e101`. On every wiring both have, the two must agree exactly (B1_COMPARE_COLS on
    every condition and stimulus), else it raises. attrs: 'missing', 'sources' ({wiring: 'E98' / 'E101'}), 'e101' (the
    record and the overlap's check)."""
    try:
        own = load_b1(wirings, roots)
    except NotImplementedError:
        own = pd.DataFrame(columns=['wiring'])
    have = set(own['wiring']) if len(own) else set()
    parts, sources, check = ([own] if len(own) else []), {w: 'E98' for w in have}, {}
    rec = None
    if e101:
        eb, rec = load_b1_e101(e101_path)
        for w in sorted(have & set(eb['wiring'])):
            a = own[own['wiring'] == w].set_index(['cond', 'stimulus'])[list(B1_COMPARE_COLS)].sort_index()
            b = eb[eb['wiring'] == w].set_index(['cond', 'stimulus'])[list(B1_COMPARE_COLS)].sort_index()
            same_idx = a.index.equals(b.index)
            d = float((a - b).abs().to_numpy().max()) if same_idx else float('nan')
            check[w] = dict(index_equal=bool(same_idx), max_abs_diff=d, n=int(len(a)))
            if not same_idx or d != 0.0:
                raise RuntimeError(f"{w}: E101's A3 rows differ from E98's B1 run: {check[w]}")
        add = [w for w in wirings if w not in have and w in set(eb['wiring'])]
        if add:
            parts.append(eb[eb['wiring'].isin(add)])
            sources.update({w: 'E101' for w in add})
    if not parts:
        raise NotImplementedError("no B1 runs (E98's or E101's)")
    out = pd.concat(parts, ignore_index=True)
    out = out[out['wiring'].isin(wirings)]
    out.attrs['missing'] = sorted(set((w, ct) for w in wirings for _, ct in CONDS)
                                  - set(zip(out['wiring'], out['cond'])))
    out.attrs['sources'] = sources
    out.attrs['e101'] = dict(record=rec, overlap_check=check)
    return out


def run_b1(label: str, root: Optional[Path] = None, force: bool = False) -> Dict[str, dict]:
    """B1 on one wiring: build the extended input (`b1_input`), run APL on and acutely off (`ExtendedRun`), save the
    outputs under `b_root(root)`/B1/ (B1_LAYOUT) with b1_stimuli.csv and a runs.csv row per condition, and the run's
    identities (`b1_identity`, <wiring>_identity.json). Skipped when the wiring's run exists in that folder (unless
    `force`). Returns the identities."""
    import gc
    d = b_root(root) / 'B1'
    d.mkdir(parents=True, exist_ok=True)
    if not force and all((d / f'{label}_{ct}_apl.csv').exists() for _, ct in CONDS):
        return _read_json(d / f'{label}_identity.json') or b1_identity(label, [d])
    t0 = time.time()
    orn, table = b1_input()
    hinfo = hallem_hc_table()[1]
    table['hallem_file_sha256'] = np.where(table['source'] == 'hallem', hinfo['sha256'], '')
    table.to_csv(d / 'b1_stimuli.csv', index=False)
    er = ExtendedRun(label, orn)
    commits = process_commits()
    for c, ct in CONDS:
        t1 = time.time()
        o = er.run(c)
        save_b1(d, label, ct, o)
        row = dict(wiring=label, cond=ct, n_stimuli=len(er.labels), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                   runtime_s=time.time() - t1, threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=_frame_sha(orn), hallem_file_sha256=hinfo['sha256'],
                   **{f'param_{k}': v for k, v in er.params(c).items()},
                   **{k: v for k, v in commits.items() if not isinstance(v, (dict, list))})
        _append_csv(d / 'runs.csv', row)
        print(f'{label} {ct}: {time.time() - t1:.0f} s', flush=True)
        del o
    del er
    gc.collect()
    ident = b1_identity(label, [d])
    _write_json(d / f'{label}_identity.json', ident)
    print(f'{label}: B1 done in {time.time() - t0:.0f} s; reproduces the saved fit and E90: {ident["reproduces"]}; '
          f'blank evoked {ident["on"]["blank_evoked"]:.4f} (settled pre {ident["on"].get("blank_evoked_settled")})', flush=True)
    return ident


def _frame_sha(df: pd.DataFrame) -> str:
    from experiments.E77_input_screen import frame_sha256
    return frame_sha256(df)


def _append_csv(path: Path, row: dict) -> None:
    df = pd.DataFrame([row])
    df.to_csv(path, mode='a', header=not path.exists(), index=False)


def b1_identity(label: str, roots: Optional[Sequence[Path]] = None) -> dict:
    """Part 0's run identities of one wiring's saved B1 run: on the 16 panel odors the spike counts and responses equal
    the saved fit's (E83's `saved_fit`, APL on and off; 0 differences), APL's window and pre-odor means equal E90's saved
    stage runs' (within APL_TOL), and the blank's evoked activity (window minus pre-odor mean) is zero within APL_TOL;
    beside, the blank's KC response rates with APL on and off."""
    from experiments import E83_cross_wiring_fit as E83
    run = load_b1_run(label, roots)
    if run is None:
        raise RuntimeError(f'{label}: no saved B1 run')
    saved = E83.saved_fit(label)['saved']
    res = dict(wiring=label, dir=run['on']['dir'])
    for c, ct in CONDS:
        o = run[ct]
        rs, ss = (H.abbrev_columns(x.astype(float)) for x in saved[c])
        s16 = o['spike_counts'][list(H.APL16_ODORS)].astype(float)
        r16 = o['responses'][list(H.APL16_ODORS)].astype(float)
        idx_ok = bool(np.array_equal(np.asarray(s16.index), np.asarray(ss.index)) if len(s16) == len(ss) else False)
        nd = int((s16.to_numpy() != ss[list(H.APL16_ODORS)].to_numpy()).sum()) if idx_ok else -1
        nr = int((r16.to_numpy() != rs[list(H.APL16_ODORS)].to_numpy()).sum()) if idx_ok else -1
        e90 = pd.read_csv(e90_dir(f'{label}_{ct}') / f'{label}_{ct}_apl.csv', index_col=0)
        a = o['apl'].loc[list(H.APL16_ODORS)]
        dwin = float((a['win'] - e90.loc[list(H.APL16_ODORS), 'win']).abs().max())
        dpre = float((a['pre'] - e90.loc[list(H.APL16_ODORS), 'pre']).abs().max())
        b = o['apl'].loc[BLANK_LABEL]
        res[c] = dict(index_equal=idx_ok, counts_differ=nd, responses_differ=nr, apl_win_max_abs_diff=dwin,
                      apl_pre_max_abs_diff=dpre, apl_mean_win=float(a['win'].mean()),
                      apl_pre_share=float(a['pre'].mean() / a['win'].mean()),
                      blank_evoked=float(b['win'] - b['pre']), blank_pre=float(b['pre']),
                      blank_evoked_settled=float(b['win'] - b['pre_settled']) if 'pre_settled' in b.index else None,
                      blank_pre_settled=float(b['pre_settled']) if 'pre_settled' in b.index else None,
                      blank_response_rate=float(o['responses'][BLANK_LABEL].astype(float).mean()),
                      n_stimuli=int(o['apl'].shape[0]))
    res['reproduces'] = bool(all(res[c]['counts_differ'] == 0 and res[c]['responses_differ'] == 0
                                 and res[c]['apl_win_max_abs_diff'] <= APL_TOL and res[c]['apl_pre_max_abs_diff'] <= APL_TOL
                                 for c, _ in CONDS))
    res['blank_zero'] = bool(abs(res['on']['blank_evoked']) <= APL_TOL)
    res['exact'] = bool(res['reproduces'] and res['blank_zero'])
    return res


def b1_identity_e101(label: str) -> dict:
    """Part 0's run identities of a wiring that E98 has not run, from E101's part A3 (`load_b1_e101`; E101's
    a3_identity.json for the 16 odors' counts and responses against the saved fit): APL's window and pre-odor means on
    the 16 panel odors against E90's saved ones, and the blank's evoked APL, as `b1_identity`."""
    eb, rec = load_b1_e101()
    ident = (_read_json(H.shared_dir() / 'results' / 'E101_apl_vs_kc_claws' / 'a3_identity.json') or {}).get(label, {})
    res = dict(wiring=label, source='E101 A3', e101_path=rec['path'],
               e101_sha256={k: v for k, v in rec['sha256'].items() if k.startswith(f'{label}_')}
               if isinstance(rec['sha256'], dict) else rec['sha256'])
    for c, ct in CONDS:
        g = eb[(eb['wiring'] == label) & (eb['cond'] == ct)].set_index('stimulus')
        e90 = pd.read_csv(e90_dir(f'{label}_{ct}') / f'{label}_{ct}_apl.csv', index_col=0)
        a = g.loc[list(H.APL16_ODORS)]
        b = g.loc[BLANK_LABEL]
        ie = ident.get(c, {})
        res[c] = dict(counts_differ=ie.get('counts_differ'), responses_differ=ie.get('responses_differ'),
                      apl_win_max_abs_diff=float((a['win'] - e90.loc[list(H.APL16_ODORS), 'win']).abs().max()),
                      apl_pre_max_abs_diff=float((a['pre'] - e90.loc[list(H.APL16_ODORS), 'pre']).abs().max()),
                      apl_mean_win=float(a['win'].mean()), apl_pre_share=float(a['pre'].mean() / a['win'].mean()),
                      blank_evoked=float(b['win'] - b['pre']), blank_evoked_settled=float(b['win'] - b['pre_settled']),
                      blank_response_rate=float(b['response_rate']), n_stimuli=int(len(g)))
    res['reproduces'] = bool(all(res[c]['counts_differ'] == 0 and res[c]['responses_differ'] == 0
                                 and res[c]['apl_win_max_abs_diff'] <= APL_TOL and res[c]['apl_pre_max_abs_diff'] <= APL_TOL
                                 for c, _ in CONDS))
    res['blank_zero'] = bool(abs(res['on']['blank_evoked']) <= APL_TOL)
    res['exact'] = bool(res['reproduces'] and res['blank_zero'])
    return res


def run_identity(label: str, run: bool = True) -> dict:
    """Part 0's run identities on one wiring: from E98's saved B1 run (`b1_identity`); else from E101's A3
    (`b1_identity_e101`); else run here first (`run_b1`, the local outputs) when `run`."""
    if load_b1_run(label) is not None:
        return b1_identity(label)
    eb, _ = load_b1_e101()
    if label in set(eb['wiring']):
        return b1_identity_e101(label)
    if not run:
        raise RuntimeError(f'{label}: no B1 run')
    run_b1(label)
    return b1_identity(label)


def part_b(arms: Sequence[str], wirings: Sequence[str], root: Optional[Path] = None) -> None:
    """Part B's runs (resumable), saved under `b_root(root)` (default outputs/E98_apl_response_strength/B/): B1 on
    `wirings` (each with its identities); B2 to B4 on the real wiring (`run_b2`, `run_b3`, `run_b4`)."""
    for arm in arms:
        if arm == 'B1':
            for w in wirings:
                ident = run_b1(w, root)
                if not ident.get('reproduces'):
                    raise RuntimeError(f'{w}: B1\'s run does not reproduce the saved fit and E90\'s APL means: {ident}')
                if not ident.get('blank_zero'):
                    print(f"{w}: the blank's evoked APL is {ident['on']['blank_evoked']:.4f} with E90's pre-odor window "
                          f"(KC layer start-up), {ident['on'].get('blank_evoked_settled')} with SETTLED_PRE_WINDOW; "
                          'recorded, part C waits on the definition', flush=True)
        elif arm in ('B2', 'B3', 'B4'):
            {'B2': run_b2, 'B3': run_b3, 'B4': run_b4}[arm](root)
        else:
            raise ValueError(f'unknown arm {arm!r}')


# B2 to B4: the real wiring's other models at fixed parameters, on B1's stimuli
ARCHIVE = ('data', 'model_outputs')
PULSE_B2 = 1.0
B_POINTS = {'B2': ('pulse1', 'E93_F', 'E95_f1'), 'B3': ('E36_spk40', 'E37_spko40'), 'B4': ('E77_inhadd1000', 'E77_inhadd3141')}
POINT_NAMES = {
    'pulse1': "C2 with the flies' 1 s pulse (E80's C2_pulse1_fixed: C2's parameters and APL->KC weights, the simulation "
              "ending 0.25 s after the pulse)",
    'E93_F': "C2 with E93's fuller input F (46 glomeruli) at E93's held point B_F_real (its breadth recomputed for F)",
    'E95_f1': "E95's refit of the family at the connectome's bouton share with the unweighted APL->KC component "
              "(rate0.068_aplboth_aplpnhold_flatu_k2_f1)",
    'E36_spk40': "E36's KC-spike-driven fit at 40 Hz with the gain (spk_thr40_gain1.5_k2_cap3_s2; older base, 5 %, "
                 "structure loss)",
    'E37_spko40': "E37's spike-only APL at that fit (the bouton->APL input off, threshold re-tuned; "
                  "spk_thr40_gain1.5_k2_cap3_s2_nopnapl)",
    'E77_inhadd1000': "C2 with E77's weaker PN normalization pn_inhadd 1000 at its held reading (S1_inhadd1000)",
    'E77_inhadd3141': "C2 with E77's weaker PN normalization pn_inhadd 3141 at its held reading (S1_inhadd3141)"}


def _archive(*parts) -> Path:
    return H.shared_dir().joinpath(*ARCHIVE, *parts)


def _saved_pair(d: Path, tag: str) -> Optional[Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]]:
    """{'on', 'off'}: (responses, spike counts) saved as <tag>_<on|acute_off>_<kind>.parquet in `d`, or None."""
    out = {}
    for c, ct in CONDS:
        fr, fs = d / f'{tag}_{ct}_responses.parquet', d / f'{tag}_{ct}_spike_counts.parquet'
        if not (fr.exists() and fs.exists()):
            return None
        out[c] = (pd.read_parquet(fr), pd.read_parquet(fs))
    return out


def _sim_cells(sim) -> Tuple[object, object, np.ndarray]:
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    mp.kc.save_inh_sims = True
    return mp, rv, np.asarray(cells['odor_order'].cell_contents, dtype=int)


def _off(p: Dict[str, float]) -> Dict[str, float]:
    return dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)


def _json_params(x) -> Dict[str, float]:
    return {k: float(v) for k, v in (json.loads(x) if isinstance(x, str) else x).items()}


def point_spec(name: str, orn: Optional[pd.DataFrame] = None) -> dict:
    """One B2 to B4 point (POINT_NAMES): 'build' (a function of nothing returning the simulator on B1's stimuli, or on
    E93's 46-row version of them for E93_F), 'params' ({'on', 'off'}: absolute parameters, as the source ran them),
    'saved' (the source's saved outputs for the identity on the 16 panel odors, or None), 'pulse_window', 'orn', 'record'
    (the source)."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E77_input_screen as S
    orn = b1_input()[0] if orn is None else orn
    if name == 'pulse1':
        from experiments import E80_pulse as E80
        P = S.c2_parts()
        p = dict(P['p'])
        kws = E80.timing_kws(PULSE_B2)
        return dict(build=lambda: S.build_c2_simulator(kws, orn=orn, breadth=P['breadth'], params=p)[0],
                    params={'on': p, 'off': _off(p)}, saved=_saved_pair(_archive('E80_remaining_candidates', 'pulse'),
                    'C2_pulse1_fixed'), pulse_window=(0.0, PULSE_B2), orn=orn, record=dict(source='E80 C2_pulse1_fixed', kws=kws))
    if name == 'E93_F':
        return e93_spec(orn)
    if name == 'E95_f1':
        r = pd.read_csv(H.shared_dir() / 'results' / 'E95_static_apl' / 'fits.csv', float_precision='round_trip')
        r = r.set_index('tag').loc['rate0.068_aplboth_aplpnhold_flatu_k2_f1']
        p = {k: float(r[k]) for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')}
        p['wKCKC_scale'] = float(r['kc_kc_inh'])
        size, unit = float(r['aplkc_flat_size']), float(r['aplkc_flat_unit'])
        if str(r['aplkc_flat']) != 'uniform' or not np.isclose(float(r['k']), 2.0):
            raise RuntimeError(f"E95's f1 refit is not the uniform component at breadth 2: {r['aplkc_flat']}, k {r['k']}")

        def build():
            from experiments import E10_joint_fit as J
            P = E83.wiring_parts('real')
            sim, info = E83.build_simulator('real', orn=orn, breadth=P['breadth'], params=p)
            w = info['cells']['wAPLKC_unscaled'].cell_contents
            comp = unit * J.aplkc_flat_component(P['wnorm']['wAPLKC'], 'uniform')
            w[...] = (info['wAPLKC'].to_numpy(dtype=float) + size * comp).reshape(w.shape)
            return sim
        return dict(build=build, params={'on': p, 'off': _off(p)}, saved=_saved_pair(_archive('E95_static_apl'),
                    'rate0.068_aplboth_aplpnhold_flatu_k2_f1'), pulse_window=None, orn=orn,
                    record=dict(source='E95 fits.csv rate0.068_aplboth_aplpnhold_flatu_k2_f1', aplkc_flat='uniform',
                                aplkc_flat_size=size, aplkc_flat_unit=unit, aplpn_hold=float(r['aplpn_hold'])))
    if name in ('E36_spk40', 'E37_spko40'):
        from experiments import E37_spike_only_apl as E37
        from pyolfsysm.mb_model import make_fixed_param_kc_simulator
        fit = 'spk_thr40_gain1.5_k2_cap3_s2'
        row = E37.fit_row(fit)
        setup = E37.setup_fit(H.apl16_orn_deltas(), local_dir('cache'), row)
        if name == 'E36_spk40':
            p = dict(setup['params'])
            saved = _saved_pair(_archive('E36_spike_apl_jointfit'), fit)
        else:
            tp = json.loads((H.shared_dir() / 'results' / 'E37_spike_only_apl' / 'tuned_params.json').read_text())
            p = {k: float(tp[f'{fit}_nopnapl'][k]) for k in E37.PARAMS}
            if p.get('wPNAPL_scale', 1.0) != 0.0:
                raise RuntimeError(f"E37's spike-only point has wPNAPL_scale {p.get('wPNAPL_scale')}")
            saved = _saved_pair(_archive('E37_spike_only_apl'), f'{fit}_nopnapl')
        kws = setup['kws']
        return dict(build=lambda: make_fixed_param_kc_simulator(orn, initial_params=p, **kws)[0],
                    params={'on': p, 'off': dict(p, **E37.OFF)}, saved=saved, pulse_window=None, orn=orn,
                    record=dict(source=f"{'E36_spike_apl_jointfit' if name == 'E36_spk40' else 'E37 tuned_params'} {fit}",
                                claw_thr=setup['claw_thr'], gain=setup['gain'], gloms=setup['gloms'], k=setup['k']))
    if name in ('E77_inhadd1000', 'E77_inhadd3141'):
        pid = 'S1_' + name.split('_')[1]
        sp = json.loads((H.shared_dir() / 'results' / 'E77_input_side_common_error' / 'screen_params.json').read_text())[pid]
        kws = {k: float(v) for k, v in sp['kws'].items()}
        breadth = S.handle_breadth(kws)
        bsha = _breadth_sha(breadth)
        if bsha != sp['breadth_sha256']:
            raise RuntimeError(f"{pid}: the breadth's sha256 {bsha} is not E77's {sp['breadth_sha256']}")
        p, poff = _json_params(sp['params']), _json_params(sp['params_off'])
        return dict(build=lambda: S.build_c2_simulator(kws, orn=orn, breadth=breadth, params=p)[0],
                    params={'on': p, 'off': poff}, saved=_saved_pair(_archive('E77_input_side_common_error', 'screen'), f'C2_{pid}'),
                    pulse_window=None, orn=orn, record=dict(source=f'E77 screen_params.json {pid}', kws=kws, breadth_sha256=bsha))
    raise ValueError(f'unknown point {name!r}')


def _breadth_sha(b: pd.Series) -> str:
    """The breadth's sha256 as E77's `hold_rates` records it (the values as little-endian float64)."""
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(b.to_numpy(dtype=float), dtype='<f8').tobytes()).hexdigest()


E93_TAG = 'E93_F_fe8e7367feec'


def e93_f_frame() -> pd.DataFrame:
    """E93's F input as its B_F_real ran it: read from results/E93_missing_glomeruli_input/inputs/F.csv (its file's sha256
    checked against runs_B.csv) with the lab input's columns, its frame sha256 checked too. (E93's `input_of('F')`
    rebuilds it within 1.4e-14 of these values, not exactly, so the saved file is read.)"""
    import hashlib
    from experiments.E77_input_screen import frame_sha256
    rb = pd.read_csv(H.shared_dir() / 'results' / 'E93_missing_glomeruli_input' / 'runs_B.csv').set_index('tag').loc['B_F_real']
    f = H.shared_dir() / 'results' / 'E93_missing_glomeruli_input' / 'inputs' / 'F.csv'
    if hashlib.sha256(f.read_bytes()).hexdigest() != rb['input_file_sha256']:
        raise RuntimeError(f'{f}: not the file B_F_real ran')
    t = pd.read_csv(f, index_col=0, float_precision='round_trip')
    lab = H.apl16_orn_deltas()
    if list(t.columns) != [H.odor_abbrev(c) for c in lab.columns]:
        raise RuntimeError(f"{f}: columns are not the lab input's odors")
    F = pd.DataFrame(t.to_numpy(dtype=float), index=pd.Index([str(g) for g in t.index], name=lab.index.name), columns=lab.columns)
    if frame_sha256(F) != rb['input_frame_sha256']:
        raise RuntimeError("E93's F input differs from the one its B_F_real ran")
    return F


def e93_input(orn: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """B1's stimuli on E93's 46 glomeruli: the 16 panel odors E93's F input (its file, sha256 checked against E93's
    runs_B.csv); the added stimuli B1's 36 rows and, on the ten added glomeruli, the session's own value where its table
    has that glomerulus (as F's lab rows), Hallem's where Hallem has it (unscaled: the stimulus is Hallem's), else 0.
    Returns the input and the added glomeruli."""
    F = e93_f_frame()
    added = [g for g in F.index if g not in orn.index]
    labels = input_labels(orn)
    sess = session_inputs()
    hal = hallem_rows(lower=False)
    specs = {x['label']: x for x in b1_specs()}
    extra = {}
    for c, lab in zip(orn.columns[16:], labels[16:]):
        sp = specs[lab]
        v = pd.Series(0.0, index=F.index)
        v.loc[orn.index] = orn[c].to_numpy(dtype=float)
        for g in added:
            if sp['source'] not in ('panel', 'blank', 'hallem'):
                t = sess[sp['session']]
                col = (sp['panel'], sp['column'])
                if g in t.index and pd.notna(t.loc[g, col]):
                    v[g] = float(t.loc[g, col])
            elif sp['source'] == 'hallem' and (sp['odor'], sp['conc']) in hal.columns and g in hal.index:
                v[g] = float(hal.loc[g, (sp['odor'], sp['conc'])])
        extra[c] = v
    out = pd.concat([F.astype(float), pd.DataFrame(extra, index=F.index)], axis='columns')
    out.columns = pd.MultiIndex.from_tuples(list(out.columns), names=orn.columns.names)
    out.index.name = F.index.name
    if input_labels(out) != input_labels(orn):
        raise RuntimeError("E93's extended input's stimuli are not B1's")
    return out, added


def e93_spec(orn: pd.DataFrame) -> dict:
    """E93_F's point: C2 (the real wiring's fit) on `e93_input`, the breadth of F's base run (E93's archived run copied
    into this experiment's cache, so nothing is simulated for it; its sha256 checked against runs_B.csv), E93's held
    parameters for B_F_real (APL on and off as recorded)."""
    import shutil
    from experiments import E10_joint_fit as J
    from experiments import E83_cross_wiring_fit as E83
    rb = pd.read_csv(H.shared_dir() / 'results' / 'E93_missing_glomeruli_input' / 'runs_B.csv').set_index('tag').loc['B_F_real']
    cache = local_dir('cache')
    for kind in ('responses', 'spike_counts'):
        f = f'input-{E93_TAG}_base_acute_off_{kind}.parquet'
        if not (cache / f).exists() and (_archive('E93_missing_glomeruli_input') / f).exists():
            shutil.copy2(_archive('E93_missing_glomeruli_input') / f, cache / f)
    F = e93_f_frame()
    _, wnorm, W, breadth = J.bouton_base(F, cache, None, input_tag=E93_TAG)
    P = E83.wiring_parts('real')
    if not W.equals(P['W']):
        raise RuntimeError("E93's F wiring is not the lab input's")
    bsha = _breadth_sha(breadth)
    if bsha != rb['breadth_sha256']:
        raise RuntimeError(f"F's breadth sha256 {bsha} is not E93's {rb['breadth_sha256']}")
    ext, added = e93_input(orn)
    p, poff = _json_params(rb['params_on']), _json_params(rb['params_off'])
    return dict(build=lambda: E83.build_simulator('real', orn=ext, breadth=breadth, params=p)[0],
                params={'on': p, 'off': poff}, saved=_saved_pair(_archive('E93_missing_glomeruli_input', 'B'), 'B_F_real'),
                pulse_window=None, orn=ext, record=dict(source='E93 runs_B.csv B_F_real', added_glomeruli=added,
                                                        breadth_sha256=bsha))


def run_point_b(arm: str, name: str, root: Optional[Path] = None, force: bool = False) -> dict:
    """One B2 to B4 point on the real wiring (`point_spec`): APL on and off at the source's parameters, the outputs and
    APL's means saved under `b_root(root)`/<arm>/ (B1_LAYOUT's files, tag the point's name), a runs.csv row per
    condition and <name>_identity.json (the 16 panel odors' counts and responses against the source's saved run, 0
    differences required; the run stops otherwise). Skipped when saved (unless `force`)."""
    import gc
    d = b_root(root) / arm
    d.mkdir(parents=True, exist_ok=True)
    if not force and all((d / f'{name}_{ct}_apl.csv').exists() for _, ct in CONDS):
        return _read_json(d / f'{name}_identity.json') or {}
    t0 = time.time()
    spec = point_spec(name)
    orn = spec['orn']
    labels = input_labels(orn)
    gl = pd.Index([str(g) for g in orn.index], name='glomerulus')
    sim = spec['build']()
    mp, rv, order = _sim_cells(sim)
    commits = process_commits()
    ident = dict(arm=arm, point=name, name=POINT_NAMES[name], record=spec['record'])
    for c, ct in CONDS:
        t1 = time.time()
        r, s = sim(spec['params'][c])
        o = relabel(r, s, labels)
        pnr = read_apl_pn(mp, rv, order, labels, _pn_index(rv, gl), spec['pulse_window'])
        o.update(pnr)
        save_b1(d, name, ct, o)
        if spec['saved'] is not None:
            rs, ss = (H.abbrev_columns(x.astype(float)) for x in spec['saved'][c])
            s16, r16 = o['spike_counts'][list(H.APL16_ODORS)], o['responses'][list(H.APL16_ODORS)]
            ok = len(s16) == len(ss) and np.array_equal(np.asarray(s16.index), np.asarray(ss.index))
            ident[c] = dict(index_equal=bool(ok),
                            counts_differ=int((s16.to_numpy(float) != ss[list(H.APL16_ODORS)].to_numpy(float)).sum()) if ok else -1,
                            responses_differ=int((r16.to_numpy(float) != rs[list(H.APL16_ODORS)].to_numpy(float)).sum()) if ok else -1)
        else:
            ident[c] = dict(index_equal=None, counts_differ=None, responses_differ=None, note='no saved run found')
        ident[c].update(n_stimuli=len(labels), blank_evoked=float(o['apl'].loc[BLANK_LABEL, 'win'] - o['apl'].loc[BLANK_LABEL, 'pre']))
        _append_csv(d / 'runs.csv', dict(arm=arm, point=name, cond=ct, n_stimuli=len(labels), host=_host(),
                                         time=time.strftime('%Y-%m-%d %H:%M'), runtime_s=time.time() - t1,
                                         threads=os.environ.get('OMP_NUM_THREADS'), input_sha256=_frame_sha(orn),
                                         **{f'param_{k}': v for k, v in spec['params'][c].items()},
                                         **{k: v for k, v in commits.items() if not isinstance(v, (dict, list))}))
        print(f'{arm} {name} {ct}: {time.time() - t1:.0f} s, identity {ident[c]}', flush=True)
        del o, r, s
    ident['exact'] = bool(all(ident[c].get('counts_differ') == 0 and ident[c].get('responses_differ') == 0 for c, _ in CONDS))
    _write_json(d / f'{name}_identity.json', ident)
    del sim, mp, rv
    gc.collect()
    if spec['saved'] is not None and not ident['exact']:
        raise RuntimeError(f'{arm} {name}: the run differs from its source on the panel: {ident}')
    print(f'{arm} {name}: done in {time.time() - t0:.0f} s', flush=True)
    return ident


def _pn_index(rv, gl: pd.Index) -> pd.Index:
    """The PN stage's glomeruli: the input's when olfsysm's PN rows are as many, else the wiring's (E90's)."""
    from pyolfsysm.mb_model import get_sim_var
    n = np.asarray(get_sim_var(rv.pn, 'pn_sims')[0]).shape[0]
    return gl if n == len(gl) else pd.RangeIndex(n, name='pn_row')


def run_b2(root: Optional[Path] = None) -> None:
    for name in B_POINTS['B2']:
        run_point_b('B2', name, root)


def run_b3(root: Optional[Path] = None) -> None:
    for name in B_POINTS['B3']:
        run_point_b('B3', name, root)


def run_b4(root: Optional[Path] = None) -> None:
    for name in B_POINTS['B4']:
        run_point_b('B4', name, root)


# ---------------------------------------------------------------------------
# the parts

E79_RESULTS = ('results', 'E79_apl_calyx_constraints')
IDENTITY_TOL = 1e-12


def _max_diff(a: pd.Series, b: pd.Series) -> Tuple[float, int]:
    """The largest |a - b| over their common index and the number of index entries in only one of them."""
    common = a.index.intersection(b.index)
    unmatched = int(len(a.index.symmetric_difference(b.index)))
    d = (a.loc[common].astype(float) - b.loc[common].astype(float)).abs()
    return (float(d.max()) if len(d) else np.nan), unmatched


def reader_identity(recs=None, resp: Optional[pd.DataFrame] = None) -> dict:
    """Part 0's reader identities: E79's per-stimulus calyx means under its three baselines
    (E79_maps_across_baselines.csv: mean_order_matched, mean_first_in_block, mean_all_pfo) and its blank table
    (E79_baseline_noise.csv: the blank's and the paraffin oil presentations' calyx means per recording;
    E79_blank_vs_odor_maps.csv: the blank and each stimulus's first-in-block mean per recording and stimulus)."""
    from experiments import E79_apl_calyx_constraints as E79
    recs = load_recordings() if recs is None else recs
    resp = stimulus_responses(recs) if resp is None else resp
    e79 = H.shared_dir().joinpath(*E79_RESULTS)
    out: Dict[str, object] = {}
    maps = pd.read_csv(e79 / 'E79_maps_across_baselines.csv').set_index(['recording', 'stimulus'])
    raw = resp[~resp['blank_subtracted'].astype(bool)]
    for bl in BASELINES:
        ours = raw[raw['baseline'] == bl].set_index(['recording', 'stimulus'])['response']
        d, u = _max_diff(ours, maps[f'mean_{bl}'])
        out[f'maps_{bl}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(ours)))
    noise = pd.read_csv(e79 / 'E79_baseline_noise.csv').set_index('recording')
    rows = {}
    for rec in recs:
        bf = E79.blank_frames(rec.frames)
        m, st = rec.mask, rec.stack
        rows[rec.name] = dict(blank_mean=float(E79.position_matched_blank(st, rec.frames)[m].mean()),
                              pfo_first_in_block_mean=float(st[bf['first_in_block']].mean(axis=0)[m].mean()),
                              pfo_block_pos2_mean=float(st[bf['pfo_pos2'][0]][m].mean()),
                              pfo_block_pos3_mean=float(st[bf['pfo_pos3'][0]][m].mean()))
    mine = pd.DataFrame(rows).T
    for col in mine.columns:
        d, u = _max_diff(mine[col], noise[col])
        out[f'blank_table_{col}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(mine)))
    bvo = pd.read_csv(e79 / 'E79_blank_vs_odor_maps.csv').set_index(['recording', 'stimulus'])
    fib = raw[raw['baseline'] == 'first_in_block'].drop_duplicates(['recording', 'stimulus']).set_index(['recording', 'stimulus'])
    for ours, theirs in (('blank', 'blank_mean'), ('response', 'odor_mean_first_in_block')):
        d, u = _max_diff(fib[ours], bvo[theirs])
        out[f'blank_vs_odor_{theirs}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(fib)))
    out['exact'] = bool(all(v['unmatched'] == 0 and np.isfinite(v['max_abs_diff']) and v['max_abs_diff'] <= IDENTITY_TOL
                            for v in out.values() if isinstance(v, dict)))
    out['tolerance'] = IDENTITY_TOL
    return out


def identity_passed(which: str = 'reader') -> bool:
    """Whether identity.json records the reader's ('reader') or every ('all') identity as passed."""
    d = _read_json(out_dir() / 'identity.json') or {}
    if which == 'reader':
        return bool(d.get('reader', {}).get('exact'))
    if which == 'all':
        return bool(d.get('all_exact'))
    raise ValueError(f'unknown identity {which!r}')


def part_0(wirings: Sequence[str] = ALL_WIRINGS, runs: bool = True) -> dict:
    """Part 0 (module docstring): identity.json. The reader's identities always; the run identities on `wirings` when
    `runs` (part B's builder)."""
    path = out_dir() / 'identity.json'
    res = _read_json(path) or {}
    t0 = time.time()
    res['reader'] = dict(reader_identity(), time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(path, res)
    print(f"reader identity: exact {res['reader']['exact']} ({time.time() - t0:.0f} s)", flush=True)
    if not res['reader']['exact']:
        raise RuntimeError(f"the reader does not reproduce E79: {res['reader']}")
    res['hallem'] = hallem_identity()
    olfdb = H.HERE / OLFDB_HALLEM_COPY
    if olfdb.exists():
        res['hallem']['olfdb_check'] = dict(file=str(OLFDB_HALLEM_COPY), result=olfdb_hallem_check(olfdb))
    _write_json(path, res)
    print(f"Hallem identity: exact {res['hallem']['exact']}", flush=True)
    if not res['hallem']['exact']:
        raise RuntimeError(f"Hallem's rows through hc_data.csv are not E38's: {res['hallem']}")
    if runs:
        for w in wirings:
            res[f'run_{w}'] = run_identity(w)
            _write_json(path, res)
    res['all_exact'] = bool(res['reader']['exact'] and res['hallem']['exact']
                            and all(res.get(f'run_{w}', {}).get('exact') for w in ALL_WIRINGS))
    _write_json(path, res)
    return res


def a1_tables(recs, resp: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """A1: the presentations' agreement per recording (under every reading: Pearson r over stimuli of the two
    presentations' calyx means, the median |rep1 - rep2| / 2, the median ratio rep2 / rep1 of responses above
    MIN_RESPONSE) and the drift of the first-in-block paraffin oil response (per recording and block its calyx mean, the
    block's stimulus; per recording the slope over blocks, the first and last; the Spearman ρ of block with concentration
    within each odor, the order's confound with concentration)."""
    from scipy.stats import spearmanr
    agree = []
    for (rec, bl, sub), g in resp.groupby(['recording', 'baseline', 'blank_subtracted']):
        ok = (g['rep1'] > MIN_RESPONSE) & (g['rep2'] > MIN_RESPONSE)
        agree.append(dict(recording=rec, fly=g['fly'].iloc[0], pulse_s=float(g['pulse_s'].iloc[0]),
                          reading=reading_label((bl, bool(sub))), n_stimuli=int(len(g)),
                          r_rep1_rep2=float(np.corrcoef(g['rep1'], g['rep2'])[0, 1]), presentation_se=presentation_se(g),
                          median_rep2_over_rep1=float(np.median(g.loc[ok, 'rep2'] / g.loc[ok, 'rep1'])) if ok.any() else np.nan,
                          mean_rep2_minus_rep1=float((g['rep2'] - g['rep1']).mean()), blank=float(g['blank'].iloc[0])))
    drift_rows, drift_sum = [], []
    for rec in recs:
        fr = rec.frames
        pfo1 = fr[(fr['odor'] == 'pfo') & (fr['position'] == 1)]
        vals = np.array([float(rec.stack[i][rec.mask].mean()) for i in pfo1.index])
        blocks = pfo1['block'].to_numpy(int)
        stim = {int(b): '; '.join(sorted({f"{o} {c:g}" for o, c in zip(g['odor'], g['log10_conc']) if o != 'pfo'}) or ['pfo'])
                for b, g in fr.groupby('block')}
        for b, v in zip(blocks, vals):
            drift_rows.append(dict(recording=rec.name, block=int(b), stimulus=stim[int(b)], pfo_first_in_block=v))
        r = resp[(resp['recording'] == rec.name) & (resp['baseline'] == 'order_matched') & ~resp['blank_subtracted'].astype(bool)]
        order = {f'spearman_block_conc_{o}': float(spearmanr(g['block'], g['conc'])[0]) for o, g in r.groupby('odor')}
        drift_sum.append(dict(recording=rec.name, n_blocks=int(len(blocks)), pfo_first_mean=float(vals.mean()),
                              pfo_first_sd=float(vals.std(ddof=1)), pfo_first_first=float(vals[0]), pfo_first_last=float(vals[-1]),
                              slope_per_block=float(np.polyfit(blocks, vals, 1)[0]),
                              spearman_block=float(spearmanr(blocks, vals)[0]), **order))
    return dict(agreement=pd.DataFrame(agree), drift=pd.DataFrame(drift_rows), drift_summary=pd.DataFrame(drift_sum))


def a2_medians(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING) -> pd.DataFrame:
    """Per odor and imaged concentration the median, minimum and maximum over recordings of the response, and the
    geometric mean over flies (responses above MIN_RESPONSE) with its 90 % interval."""
    rows = []
    for (o, c), g in _reading(resp, reading).groupby(['odor', 'conc']):
        gm = fly_geomean(g, 'response')
        rows.append(dict(odor=o, conc=float(c), reading=reading_label(reading), n_recordings=int(len(g)),
                         median=float(g['response'].median()), min=float(g['response'].min()), max=float(g['response'].max()),
                         n_above=int((g['response'] > MIN_RESPONSE).sum()), fly_geomean=gm['ratio'], fly_lo90=gm['lo90'],
                         fly_hi90=gm['hi90'], n_flies=gm['n_flies']))
    return pd.DataFrame(rows)


def a2_fly_summary(dose: pd.DataFrame) -> pd.DataFrame:
    """Per odor the fold change per log unit over flies (geometric mean of the recordings' folds within flies, then over
    flies, with the 90 % interval) and the lowest concentration above the blank (median and range over recordings)."""
    rows = []
    for (o, rd), g in dose.groupby(['odor', 'reading']):
        gm = fly_geomean(g, 'fold_per_log_unit')
        rows.append(dict(odor=o, reading=rd, fold_per_log_unit=gm['ratio'], lo90=gm['lo90'], hi90=gm['hi90'],
                         n_flies=gm['n_flies'], fold_min=float(g['fold_per_log_unit'].min()),
                         fold_max=float(g['fold_per_log_unit'].max()),
                         lowest_above_blank_median=float(g['lowest_above_blank'].median()),
                         lowest_above_blank_min=float(g['lowest_above_blank'].min()),
                         lowest_above_blank_max=float(g['lowest_above_blank'].max()),
                         n_recordings_none_above=int(g['lowest_above_blank'].isna().sum())))
    return pd.DataFrame(rows)


def a4_model_ratios(apl: pd.DataFrame) -> pd.DataFrame:
    """The model's ratios that the panel's runs give, per wiring, n and form: paa over ms (the lab's input, paa at
    10^-4 and ms at 10^-3: V1 under the megamat estimate of ms's input) and pa over ms (pa at 10^-3; beside the flies'
    pa at 10^-4.386, outside the one-log rule)."""
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        a0 = float(g['a0'].iloc[0])
        for rid, a, b in ACROSS_ODORS:
            rs = model_ratio_set(g.loc[a, 'a'], g.loc[b, 'a'], a0)
            for _, r in rs.iterrows():
                rows.append(dict(wiring=w, ratio_id=rid, n=r['n'], form=r['form'], ratio=r['ratio']))
    return pd.DataFrame(rows)


PREVIEW = (('V1 under the megamat estimate of ms\'s input', 'paa', 'ms', ('V1_paa-4', 'V1_paa-2.954')),
           ('pa over ms (the model\'s pa at 10^-3, the flies\' at 10^-4.386, outside the one-log rule)', 'pa', 'ms',
            ('probe_pa-4.386',)))


PREVIEW_SUBSETS = ('all', 'all_recording_ranges')


def a4_preview(apl: pd.DataFrame, a3: pd.DataFrame) -> pd.DataFrame:
    """The comparisons that the panel's runs already allow, under every reading of the flies (subsets PREVIEW_SUBSETS): per wiring
    the class (`classify`, n and both forms) and the linear reading's, and the comparison's verdict. V1 here is under one
    of ms's estimates only: inside here means V1 is inside; outside here leaves V1 to part C."""
    rows = []
    a3 = a3[a3['subset'].isin(PREVIEW_SUBSETS)]
    for label, num, den, ids in PREVIEW:
        for (sub, rd), g in a3.groupby(['subset', 'reading'], sort=False):
            g = g.set_index('ratio_id')
            iv = [(g.loc[i, 'lo90'], g.loc[i, 'hi90']) for i in ids]
            if any(not np.isfinite(x) for v in iv for x in v):
                continue
            lo, hi = hull(iv)
            wc = wiring_classes(apl, num, den, lo, hi)
            rows.append(dict(comparison=label, reading=rd, subset=sub, flies_lo90=lo, flies_hi90=hi,
                             verdict=comparison_verdict(wc['cls']), verdict_linear=comparison_verdict(wc['cls_linear']),
                             **{f'n_{c}': int((wc['cls'] == c).sum()) for c in ('inside', 'flatter', 'steeper', 'mixed')},
                             model_linear_min=float(wc['ratio_linear'].min()), model_linear_max=float(wc['ratio_linear'].max()),
                             model_min=float(wc['ratio_min'].min()), model_max=float(wc['ratio_max'].max()),
                             per_wiring='; '.join(f"{r.wiring} {r.cls} ({r.ratio_linear:.3f})" for r in wc.itertuples())))
    return pd.DataFrame(rows)


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def input_files(wirings: Sequence[str] = ALL_WIRINGS) -> Dict[str, str]:
    """Part A's inputs (path relative to the shared directory or the repository root -> sha256): the calyx maps (TIFFs,
    ROIs, recordings.csv, frames.csv), E79's tables the reader is checked against, E90's stage files read, E37's APL
    table and the data-target registry."""
    from experiments import E79_apl_calyx_constraints as E79
    sh = H.shared_dir()
    files = sorted(E79.DATA.glob('*.tif')) + sorted(E79.DATA.glob('*.roi')) + [E79.DATA / 'recordings.csv', E79.DATA / 'frames.csv']
    e79 = sh.joinpath(*E79_RESULTS)
    files += [e79 / f for f in ('E79_maps_across_baselines.csv', 'E79_baseline_noise.csv', 'E79_blank_vs_odor_maps.csv')]
    for w in wirings:
        files += [e90_dir(f'{w}_on') / f'{w}_on_{k}' for k in ('apl.csv', 'pn_win.parquet', 'pn_pre.parquet', 'responses.parquet')]
        files.append(e90_dir(f'{w}_acute_off') / f'{w}_acute_off_responses.parquet')
    files += [sh / 'results' / E37_TABLE[0] / E37_TABLE[1], H.data_targets_path()]
    out = {}
    for f in files:
        f = Path(f)
        key = str(f.relative_to(sh)) if f.is_relative_to(sh) else str(f.relative_to(H.HERE)) if f.is_relative_to(H.HERE) else str(f)
        out[key] = _sha256(f)
    return out


def _records(df: pd.DataFrame, cols: Sequence[str]) -> List[dict]:
    return [{c: (r[c] if not isinstance(r[c], (np.floating, np.integer)) else r[c].item()) for c in cols} for _, r in df.iterrows()]


def part_a() -> dict:
    """Part A (module docstring): A1 to A6 and a_summary.json."""
    if not identity_passed('reader'):
        raise RuntimeError("part A needs part 0's reader identities (run --part 0 --no-runs first)")
    out = out_dir()
    t0 = time.time()
    recs = load_recordings()
    resp = stimulus_responses(recs)
    resp.to_csv(out / 'a1_responses.csv', index=False)
    a1 = a1_tables(recs, resp)
    a1['agreement'].to_csv(out / 'a1_agreement.csv', index=False)
    a1['drift'].to_csv(out / 'a1_pfo_drift.csv', index=False)
    a1['drift_summary'].to_csv(out / 'a1_pfo_drift_summary.csv', index=False)
    print(f'A1 done ({time.time() - t0:.0f} s)', flush=True)
    dose, at, med = [], [], []
    for rd in READINGS:
        d, a = dose_response(resp, rd)
        dose.append(d)
        at.append(a)
        med.append(a2_medians(resp, rd))
    dose, at, med = (pd.concat(x, ignore_index=True) for x in (dose, at, med))
    dose.to_csv(out / 'a2_dose_response.csv', index=False)
    at.to_csv(out / 'a2_at_model_concs.csv', index=False)
    med.to_csv(out / 'a2_medians.csv', index=False)
    a2s = a2_fly_summary(dose)
    a2s.to_csv(out / 'a2_fly_summary.csv', index=False)
    print(f'A2 done ({time.time() - t0:.0f} s)', flush=True)
    a3, a3r = [], []
    for rd in READINGS:
        for sub in SUBSETS:
            s, r = flies_ratios(resp, rd, sub)
            a3.append(s)
            a3r.append(r)
    a3, a3r = pd.concat(a3, ignore_index=True), pd.concat(a3r, ignore_index=True)
    a3.to_csv(out / 'a3_ratios.csv', index=False)
    a3r.to_csv(out / 'a3_per_recording.csv', index=False)
    print(f'A3 done ({time.time() - t0:.0f} s)', flush=True)
    apl = family_apl()
    apl.to_csv(out / 'a4_apl.csv', index=False)
    prof = profile_summary(apl)
    prof.to_csv(out / 'a4_profile.csv', index=False)
    mr = a4_model_ratios(apl)
    mr.to_csv(out / 'a4_model_ratios.csv', index=False)
    prev = a4_preview(apl, a3)
    prev.to_csv(out / 'a4_preview.csv', index=False)
    print(f'A4 done ({time.time() - t0:.0f} s)', flush=True)
    alt = alternative_apls(apl)
    alt.to_csv(out / 'a5_alternatives.csv', index=False)
    print(f'A5 done ({time.time() - t0:.0f} s)', flush=True)
    main = a3[(a3['reading'] == reading_label(MAIN_READING)) & (a3['subset'] == 'all')]
    literal = a3[(a3['reading'] == reading_label(('order_matched', True))) & (a3['subset'] == 'all')]
    sds = power_sds(main)
    power = power_table(apl, {k: v for k, v in sds.items() if k in ('median', 'max')})
    power.to_csv(out / 'a6_power.csv', index=False)
    print(f'A6 done ({time.time() - t0:.0f} s)', flush=True)
    ratio_cols = ['ratio_id', 'numerator', 'denominator', 'ratio', 'lo90', 'hi90', 'sd_log2_flies', 'n_flies',
                  'n_recordings', 'n_recordings_left_out', 'num_how', 'den_how', 'probe_ratio', 'probe_lo90', 'probe_hi90']
    summ = dict(
        reader_identity=(_read_json(out / 'identity.json') or {}).get('reader', {}).get('exact'),
        blank_identity_max_abs_diff=blank_identity(resp),
        large_blank_recordings=sorted(resp.loc[resp['blank'] > LARGE_BLANK, 'recording'].unique()),
        a3_main=_records(main if len(main) else a3[a3['reading'] == reading_label(MAIN_READING)], ratio_cols),
        a3_literal_blank_twice=_records(literal, ratio_cols),
        a4_profile=_records(prof, ['wiring', 'a_mean', 'a_cv', 'a0_share', 'evoked_mean', 'evoked_cv', 'evoked_min_odor',
                                   'evoked_max_odor', 'ms_over_mean', 'spearman_evoked_input_total',
                                   'spearman_evoked_pn_total', 'spearman_evoked_kc_rate_off']),
        a4_preview=_records(prev, ['comparison', 'reading', 'subset', 'flies_lo90', 'flies_hi90', 'verdict', 'verdict_linear',
                                   'n_inside', 'n_flatter', 'n_steeper', 'model_linear_min', 'model_linear_max']),
        a5=_records(alt[~alt['source'].str.startswith('E37:') | alt['source'].str.contains('_mean')],
                    ['source', 'kind', 'ratio_id', 'ratio', 'ratio_min', 'ratio_max']),
        a6_sd_log2=sds, a6_best5=_records(power.head(5), list(power.columns)),
        a6_best5_without_substituted_input=_records(power[~power['substituted_input']].head(5), list(power.columns)),
        a6_n_pairs_flies_needed_at_most_3=dict(
            median=int((power['flies_needed_sd_median'] <= 3).sum()), max=int((power['flies_needed_sd_max'] <= 3).sum())),
        runtime_s=time.time() - t0)
    _write_json(out / 'a_summary.json', summ)
    _write_json(out / 'provenance.json', dict(experiment=EXPERIMENT, script=H.script_path(__file__), parts=['0', 'A'],
        host=_host(), time=time.strftime('%Y-%m-%d %H:%M'), commits=dict(H.code_commits()), inputs=input_files(),
        settings=dict(MIN_RESPONSE=MIN_RESPONSE, MATCH_LOG_UNITS=MATCH_LOG_UNITS, N_BOOT=N_BOOT, SEED=SEED, CI=CI,
                      INDICATOR_NS=INDICATOR_NS, FORMS=FORMS, VERDICT_MIN_WIRINGS=VERDICT_MIN_WIRINGS,
                      MAIN_READING=MAIN_READING, LARGE_BLANK=LARGE_BLANK, POWER_Z=POWER_Z, POWER_MIN_CONC=POWER_MIN_CONC,
                      RATIOS=[list(r) for r in RATIOS])))
    return summ


# the verdict's comparisons on B1's stimuli (the model's labels)
V1_NUM = 'paa'                                   # the lab's input, paa at 10^-4 (validation2's panel)
V1_ESTIMATES: Dict[str, str] = {                 # the estimates of ms's input at 10^-3: name -> B1 label
    'megamat/megamat (the panel)': 'ms', 'megamat/diag': session_label('ms', -3, 'megamat', DIAG),
    'validation2/diag': session_label('ms', -3, 'validation2', DIAG),
    'kiwi-control/diag': session_label('ms', -3, 'kiwi-control', DIAG),
    'kiwi-control/control': session_label('ms', -3, 'kiwi-control', 'control')}
V2_NUM, V2_DEN = session_label('paa', -5, PAA5_SESSION, DIAG), 'paa'
V2_BESIDE = {sn: session_label('paa', -5, sn, DIAG) for sn in SESSIONS if sn != PAA5_SESSION}
# V3 to V5 with Hallem's input (described, decide nothing): ratio id of RATIOS -> (numerator, denominator) B1 labels
HALLEM_COMPARISONS: Dict[str, Tuple[str, str]] = {
    'V3_ms-4': (hallem_label('ms', -4), hallem_label('ms', -2)), 'V3_ms-6': (hallem_label('ms', -6), hallem_label('ms', -2)),
    'V3_ms-8': (hallem_label('ms', -8), hallem_label('ms', -2)), 'V4': (hallem_label('pa', -6), hallem_label('pa', -4)),
    'V5_-4': (hallem_label('pa', -4), hallem_label('ms', -4)), 'V5_-6': (hallem_label('pa', -6), hallem_label('ms', -6)),
    'V5_paa-2': (hallem_label('paa', -2), hallem_label('ms', -2))}
FLIES_SOLVENT_RATE = {'on': ('ctrl', 'control flies'), 'off': ('tnt', 'APL output blocked (TNT)')}
OUTCOMES = {'a': 'V1 and V2 inside', 'b': 'V2 flatter (in ascending series)', 'c': 'V2 steeper (in ascending series)',
            'd': "V1 outside under every estimate of ms's input",
            'e': 'V1 or V2 not decided (reading-dependent, short of 6 of 8, or the 0.5 s and 1 s recordings disagree)',
            'h': "the blank's response rate with APL off above the flies' solvent rate with APL output blocked (TNT)"}
# the flies' readings (the amendment of 2026-11-15 15:24): the carrying reading, its partners, the entry's literal reading
CARRYING = (('order_matched', False), 'all')
PARTNERS: Dict[str, Tuple[Tuple[str, bool], str]] = {
    'first_in_block': (('first_in_block', False), 'all'),              # the other baseline, no blank removed
    'own_range': (('order_matched', False), 'all_recording_ranges'),     # ms at 10^-3.431 in the two recordings; V1 only
    'without_three': (('order_matched', False), 'without_large_blank_and_3FRfun1')}
V1_ONLY_PARTNERS = ('own_range',)
LITERAL = (('order_matched', True), 'all')                               # the blank removed twice: decides nothing
READING_SETS = dict(carrying=CARRYING, **PARTNERS, literal=LITERAL)
CATEGORIES = ('inside', 'flatter', 'steeper')
PULSE_SHARE = 0.1                    # the pulse sets disagree when fewer than this share of the splits are as large
REPEAT_RECORDINGS = ('20210805_3FRfun1', '20210823_4FRfun2')     # pa at 10^-6.771 delivered twice
SHUFFLED_RECORDING = '20210827_3FRfun4'                          # the one recording with shuffled order


def stimulus_classes(on: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring the model's ratio of stimulus `num` over stimulus `den` (rows of `on`: wiring, stimulus, win, pre; APL
    on) through the indicator family (INDICATOR_NS x FORMS), its range, the linear reading (n 1, evoked), the class of
    the set and of the linear reading against [lo, hi] (`classify`), and the distance from the nearer end of the
    interval: of the linear reading (signed, the linear ratio minus that end) and of the set (0 if any of it is
    inside). The pre-odor level must be the same for both stimuli."""
    rows = []
    for w, g in on.groupby('wiring', sort=False):
        g = g.set_index('stimulus')
        if num not in g.index or den not in g.index:
            continue
        a1, a2 = float(g.loc[num, 'win']), float(g.loc[den, 'win'])
        p1, p2 = float(g.loc[num, 'pre']), float(g.loc[den, 'pre'])
        if not np.isclose(p1, p2, rtol=1e-9, atol=1e-9):
            raise RuntimeError(f'{w}: the pre-odor levels of {num} and {den} differ ({p1}, {p2})')
        rs = model_ratio_set(a1, a2, p2)
        lin = float(indicator_ratio(a1, a2, p2, 1.0, 'evoked'))
        end = lo if abs(lin - lo) <= abs(lin - hi) else hi
        x = rs['ratio'].to_numpy(float)
        inside = (x >= lo) & (x <= hi)
        rows.append(dict(wiring=w, numerator=num, denominator=den, a_num=a1, a_den=a2, a0=p2, ratio_linear=lin,
                         ratio_min=float(x.min()), ratio_max=float(x.max()), flies_lo90=lo, flies_hi90=hi,
                         nearer_end='lo' if end == lo else 'hi', linear_minus_nearer_end=lin - end,
                         set_distance=0.0 if inside.any() else float(np.min(np.minimum(np.abs(x - lo), np.abs(x - hi)))),
                         nearest_ratio=float(x[np.argmin(np.minimum(np.abs(np.log(x / lo)), np.abs(np.log(x / hi))))]),
                         cls=classify(x, lo, hi), cls_linear=classify([lin], lo, hi)))
    return pd.DataFrame(rows)


def held_category(categories: Dict[str, str], partners: Sequence[str]) -> str:
    """'inside', 'flatter' or 'steeper' when the carrying reading ('carrying') and every partner in `partners` give that
    category (each by the 6-of-8 rule); otherwise 'reading-dependent'."""
    c = categories['carrying']
    if c in CATEGORIES and all(categories[k] == c for k in partners):
        return c
    return 'reading-dependent'


def v1_from_estimates(held: Dict[str, str]) -> str:
    """V1 from the held category under each estimate of ms's input: 'inside' if inside is held under at least one;
    'flatter' / 'steeper' if that is held under every one ('outside, mixed' if each is held flatter or steeper but not
    all the same); otherwise 'not decided'."""
    v = list(held.values())
    if not v:
        raise ValueError("no estimate of ms's input")
    if 'inside' in v:
        return 'inside'
    if all(x in ('flatter', 'steeper') for x in v):
        return v[0] if len(set(v)) == 1 else 'outside, mixed'
    return 'not decided'


def pulse_permutation(per_fly: pd.Series, short: Sequence[str]) -> dict:
    """Whether the 0.5 s and the 1 s recordings disagree on a ratio: `per_fly` the flies' mean log ratios, `short` the
    flies of the 0.5 s set; over every split of the flies into sets of the same sizes (56 for 3 and 5), the share whose
    difference of mean log ratios is at least as large (absolute) as the observed one; they disagree when that share is
    below PULSE_SHARE."""
    from itertools import combinations
    x = per_fly.dropna()
    flies = list(x.index)
    k = sum(f in set(short) for f in flies)
    if k == 0 or k == len(flies):
        return dict(n_flies=len(flies), n_short=k, n_splits=0, n_as_large=0, share=np.nan, observed=np.nan, disagree=False)
    v = x.to_numpy(float)
    obs = abs(np.mean([x[f] for f in flies if f in set(short)]) - np.mean([x[f] for f in flies if f not in set(short)]))
    n, n_as = 0, 0
    for idx in combinations(range(len(flies)), k):
        m = np.zeros(len(flies), bool)
        m[list(idx)] = True
        d = abs(v[m].mean() - v[~m].mean())
        n += 1
        n_as += int(d >= obs - 1e-12)
    return dict(n_flies=len(flies), n_short=k, n_splits=n, n_as_large=n_as, share=n_as / n, observed=float(obs),
                ratio_short=float(np.exp(np.mean([x[f] for f in flies if f in set(short)]))),
                ratio_long=float(np.exp(np.mean([x[f] for f in flies if f not in set(short)]))),
                disagree=bool(n_as / n < PULSE_SHARE))


def _intervals(flies: pd.DataFrame, ids: Sequence[str]) -> Optional[Tuple[float, float]]:
    f = flies.set_index('ratio_id') if 'ratio_id' in flies.columns else flies
    if any(i not in f.index for i in ids):
        return None
    iv = [(float(f.loc[i, 'lo90']), float(f.loc[i, 'hi90'])) for i in ids]
    if any(not (np.isfinite(a) and np.isfinite(b)) for a, b in iv):
        return None
    return hull(iv)


def verdict_core(flies: Dict[str, pd.DataFrame], b1: pd.DataFrame, pulse: Optional[Dict[str, dict]] = None,
                 solvent: Optional[Dict[str, float]] = None) -> Tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Part C's verdict (the entry's (C) as amended 2026-11-15 15:24). `flies`: the flies' ratios (rows of `flies_ratios`)
    under every reading of READING_SETS ('carrying', the partners, 'literal'); `b1`: `load_b1` rows (APL on is read);
    `pulse`: `pulse_permutation` results for 'V1_paa-4', 'V1_paa-2.954' and 'V2' (carrying reading); `solvent`: the
    flies' solvent response rates ({'ctrl', 'tnt'}) for outcome (h).
    - Each comparison is read under every reading (the 6-of-8 rule on the same model rows); its category is held when
      the carrying reading and every partner give it (`held_category`; V2 without the V1-only partner), otherwise
      reading-dependent (not decided). The literal reading is shown and decides nothing.
    - V1: paa (the lab's input at 10^-4) over each estimate of ms's input at 10^-3 (V1_ESTIMATES, the panel's own first)
      against the hull of the flies' two readings of paa's concentration; `v1_from_estimates` over the held categories.
    - V2: paa at 10^-5 (PAA5_SESSION's diagnostics) over the panel's paa; the other sessions' estimates beside.
    - A comparison on which the pulse sets disagree is not decided (V1 if either of its ratios disagrees).
    - V3 to V5 on Hallem's stimuli where B1 has them, under the carrying reading (described).
    Returns (verdict, one row per comparison and reading, the per-wiring rows of every comparison and reading)."""
    on = b1[b1['cond'] == 'on'][['wiring', 'stimulus', 'win', 'pre']]
    comps, rows = [], []

    def read(name, num, den, ids, readings, decides):
        out = {}
        for rk in readings:
            iv = _intervals(flies[rk], ids)
            if iv is None:
                continue
            wc = stimulus_classes(on, num, den, *iv)
            if wc.empty:
                continue
            wc.insert(0, 'reading_set', rk)
            wc.insert(0, 'comparison', name)
            rows.append(wc)
            c = dict(comparison=name, reading_set=rk, reading=reading_label(READING_SETS[rk][0]), subset=READING_SETS[rk][1],
                     numerator=num, denominator=den, flies_lo90=iv[0], flies_hi90=iv[1], n_wirings=int(len(wc)),
                     category=comparison_verdict(wc['cls']), category_linear=comparison_verdict(wc['cls_linear']),
                     **{f'n_{k}': int((wc['cls'] == k).sum()) for k in CATEGORIES + ('mixed',)},
                     model_linear_min=float(wc['ratio_linear'].min()), model_linear_max=float(wc['ratio_linear'].max()),
                     model_min=float(wc['ratio_min'].min()), model_max=float(wc['ratio_max'].max()), decides=decides)
            comps.append(c)
            out[rk] = c
        return out

    v1_partners = tuple(PARTNERS)
    v2_partners = tuple(k for k in PARTNERS if k not in V1_ONLY_PARTNERS)
    v1_ids = ('V1_paa-4', 'V1_paa-2.954')
    v1_held, v1_lines = {}, {}
    for est, lab in V1_ESTIMATES.items():
        r = read(f"V1, ms's input from {est}", V1_NUM, lab, v1_ids, ('carrying',) + v1_partners + ('literal',), True)
        if 'carrying' not in r or any(k not in r for k in v1_partners):
            continue
        cats = {k: r[k]['category'] for k in r}
        v1_held[est] = held_category(cats, v1_partners)
        v1_lines[est] = _line(r, v1_partners, v1_held[est])
    r2 = read(f'V2, paa at 10^-5 from {PAA5_SESSION}/diag', V2_NUM, V2_DEN, ('V2',), ('carrying',) + v2_partners + ('literal',), True)
    for sn, lab in V2_BESIDE.items():
        read(f'V2 beside, paa at 10^-5 from {sn}/diag', lab, V2_DEN, ('V2',), ('carrying',) + v2_partners + ('literal',), False)
    for rid, (num, den) in HALLEM_COMPARISONS.items():
        read(f"{rid} (Hallem's input)", num, den, (rid,), ('carrying',), False)
    pulse = pulse or {}
    v1_pulse = any(pulse.get(i, {}).get('disagree') for i in v1_ids)
    v2_pulse = bool(pulse.get('V2', {}).get('disagree'))
    v1 = v1_from_estimates(v1_held) if v1_held else 'not run'
    if v1_pulse and v1 != 'not run':
        v1 = 'not decided'
    if 'carrying' in r2 and all(k in r2 for k in v2_partners):
        v2_held = held_category({k: r2[k]['category'] for k in r2}, v2_partners)
        v2 = 'not decided' if (v2_pulse or v2_held == 'reading-dependent') else v2_held
        v2_line = _line(r2, v2_partners, v2_held)
    else:
        v2_held, v2, v2_line = 'not run', 'not run', {}
    out = dict(V1=v1, V1_held_per_estimate=v1_held, V1_lines=v1_lines, V1_pulse_disagree=v1_pulse, V2=v2, V2_held=v2_held,
               V2_line=v2_line, V2_pulse_disagree=v2_pulse, pulse=pulse,
               V2_note='V2 is read in the series as recorded, in which concentration rises with block (11 of 12 recordings)')
    oc = []
    if v1 == 'inside' and v2 == 'inside':
        oc.append('a')
    if v2 == 'flatter':
        oc.append('b')
    if v2 == 'steeper':
        oc.append('c')
    if v1 in ('flatter', 'steeper', 'outside, mixed'):
        oc.append('d')
    if v1 in ('not decided',) or v2 in ('not decided',) or v1_pulse or v2_pulse:
        oc.append('e')
    blank = b1[b1['stimulus'] == BLANK_LABEL]
    if len(blank):
        rates = blank.groupby('cond')['response_rate'].agg(['min', 'median', 'max'])
        out['blank_response_rate'] = {c: rates.loc[c].to_dict() for c in rates.index}
        if solvent is not None:
            out['flies_solvent_rate'] = dict(solvent)
            if 'acute_off' in rates.index and rates.loc['acute_off', 'median'] > solvent['tnt']:
                oc.append('h')
    out['outcomes'] = oc
    out['outcome_names'] = {k: OUTCOMES[k] for k in oc}
    allrows = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(allrows):
        c2 = allrows[(allrows['comparison'].str.startswith('V2, ')) & (allrows['reading_set'] == 'carrying')]
        if len(c2):
            out['V2_order_factor'] = order_factor(c2)
    return out, pd.DataFrame(comps), allrows


def _line(r: Dict[str, dict], partners: Sequence[str], held: str) -> dict:
    """A verdict line: the carrying reading's category with its wiring count, each partner's, the literal reading's, and
    the held category."""
    def cell(c):
        return f"{c['category']} ({c.get('n_' + c['category'], 0) if c['category'] in CATEGORIES else max(c[f'n_{k}'] for k in CATEGORIES)} of {c['n_wirings']})"
    line = dict(carrying=cell(r['carrying']), **{k: cell(r[k]) for k in partners if k in r})
    if 'literal' in r:
        line['literal (decides nothing)'] = cell(r['literal'])
    line['held'] = held
    return line


def order_factor(rows: pd.DataFrame) -> dict:
    """V2's order description: the factor between the nearer end of the flies' interval and the model's nearest ratio
    over the wirings and the indicator family (1 when a wiring is inside): the order effect that would be needed."""
    lo, hi = float(rows['flies_lo90'].iloc[0]), float(rows['flies_hi90'].iloc[0])
    if (rows['set_distance'] == 0).any():
        return dict(factor=1.0, note='a wiring is inside')
    x = rows['nearest_ratio'].to_numpy(float)
    i = int(np.argmin(np.minimum(np.abs(np.log(x / lo)), np.abs(np.log(x / hi)))))
    end = lo if abs(np.log(x[i] / lo)) <= abs(np.log(x[i] / hi)) else hi
    return dict(factor=float(end / x[i]), nearest_model_ratio=float(x[i]), nearer_end=float(end),
                wiring=str(rows['wiring'].iloc[i]))


def pa_repeat(recs) -> pd.DataFrame:
    """V2's description: pa at 10^-6.771 delivered twice (REPEAT_RECORDINGS), the later block's response over the
    earlier's under the carrying baseline (each presentation minus the paraffin oil block's at the same position, the
    calyx mean, the mean of the two presentations)."""
    from experiments import E79_apl_calyx_constraints as E79
    rows = []
    for rec in recs:
        if rec.name not in REPEAT_RECORDINGS:
            continue
        fr, st, mask = rec.frames, rec.stack, rec.mask
        b = E79.blank_frames(fr)
        by_pos = {2: st[b['pfo_pos2'][0]], 3: st[b['pfo_pos3'][0]]}
        g = fr[(fr['odor'] == 'pa') & np.isclose(fr['log10_conc'].astype(float), -6.771)]
        resp = {int(blk): float(np.mean([(st[i] - by_pos[int(pos)])[mask].mean() for i, pos in zip(gg.index, gg['position'])]))
                for blk, gg in g.groupby('block')}
        if len(resp) != 2:
            continue
        (b0, r0), (b1_, r1) = sorted(resp.items())
        rows.append(dict(recording=rec.name, earlier_block=b0, later_block=b1_, earlier=r0, later=r1, later_over_earlier=r1 / r0))
    return pd.DataFrame(rows)


def v2_descriptions(per_recording: pd.DataFrame, recs, drift: pd.DataFrame) -> dict:
    """V2's descriptions beside the verdict (no gate): SHUFFLED_RECORDING's V2 ratio alone against the geometric mean
    and range of the other recordings' (carrying reading, `per_recording`: `flies_ratios`' rows of V2); the pa 10^-6.771
    repeat (`pa_repeat`); the first-in-block paraffin oil drift (A1's summary: recordings up and down, the median
    slope). The order factor is `order_factor`, from the model's rows."""
    v = per_recording[(per_recording['ratio_id'] == 'V2') & per_recording['entered']].set_index('recording')['log_ratio']
    out = {}
    if SHUFFLED_RECORDING in v.index:
        others = v.drop(SHUFFLED_RECORDING)
        out['shuffled_recording'] = dict(recording=SHUFFLED_RECORDING, ratio=float(np.exp(v[SHUFFLED_RECORDING])),
                                         others_geomean=float(np.exp(others.mean())), others_min=float(np.exp(others.min())),
                                         others_max=float(np.exp(others.max())), n_others=int(len(others)))
    rep = pa_repeat(recs)
    out['pa_repeat'] = rep.to_dict('records')
    out['pfo_drift'] = dict(n_up=int((drift['slope_per_block'] > 0).sum()), n_down=int((drift['slope_per_block'] < 0).sum()),
                            median_slope_per_block=float(drift['slope_per_block'].median()))
    return out


def pulse_matched(b1: pd.DataFrame, flies_short: pd.DataFrame, flies_long: pd.DataFrame, root: Optional[Path] = None) -> dict:
    """The pulse-matched readings (described): V1 under the megamat estimate of ms's input, the 0.5 s recordings' hull
    against the family's eight wirings, and the 1 s recordings' hull against B2's 1 s run on the real wiring (APL's
    mean over 0-0.5 s and over 0-1 s) when it exists."""
    on = b1[b1['cond'] == 'on'][['wiring', 'stimulus', 'win', 'pre']]
    out = {}
    iv = _intervals(flies_short, ('V1_paa-4', 'V1_paa-2.954'))
    if iv is not None:
        wc = stimulus_classes(on, V1_NUM, 'ms', *iv)
        out['short_vs_family'] = dict(flies=list(iv), n_inside=int((wc['cls'] == 'inside').sum()), n_wirings=int(len(wc)),
                                      category=comparison_verdict(wc['cls']))
    iv = _intervals(flies_long, ('V1_paa-4', 'V1_paa-2.954'))
    d = b_root(root) / 'B2'
    if iv is not None and (d / 'pulse1_on_apl.csv').exists():
        a = pd.read_csv(d / 'pulse1_on_apl.csv', index_col=0)
        res = dict(flies=list(iv))
        for col in ('win', 'win_pulse'):
            rs = model_ratio_set(a.loc[V1_NUM, col], a.loc['ms', col], a.loc['ms', 'pre'])
            res[col] = dict(ratio_linear=float(rs['ratio'].iloc[0]), ratio_min=float(rs['ratio'].min()),
                            ratio_max=float(rs['ratio'].max()), cls=classify(rs['ratio'], *iv))
        out['long_vs_b2_pulse1'] = res
    return out


def with_a0(b1: pd.DataFrame, how: str) -> pd.DataFrame:
    """B1's rows with the pre-odor level ('pre', and 'evoked') set by `how` (A0_DEFINITIONS): 'e90' E90's pre-odor mean
    over -0.5 to 0 s (as saved; it includes the KC layer's start-up); 'settled' the mean over SETTLED_PRE_WINDOW; 'blank'
    the blank stimulus's window mean in the same wiring and condition (the model's paraffin oil)."""
    if how not in A0_DEFINITIONS:
        raise ValueError(f'unknown pre-odor level {how!r}')
    b = b1.copy()
    if how == 'settled':
        if b['pre_settled'].isna().any():
            raise RuntimeError('runs without pre_settled (saved before it was added): re-run them')
        b['pre'] = b['pre_settled']
    elif how == 'blank':
        bl = b[b['stimulus'] == BLANK_LABEL].set_index(['wiring', 'cond'])['win']
        b['pre'] = [bl[(w, c)] for w, c in zip(b['wiring'], b['cond'])]
    b['evoked'] = b['win'] - b['pre']
    return b


def part_c(wirings: Sequence[str] = ALL_WIRINGS, allow_missing: bool = False, a0: Optional[str] = None) -> dict:
    """Part C (module docstring; the entry's (C) as amended 2026-11-15 15:24): reads B1's runs (`load_b1`; all eight
    wirings unless `allow_missing`), the flies' ratios under every reading of READING_SETS (`flies_ratios`), the pulse
    sets' permutation test on V1's two ratios and V2 (carrying reading), V2's descriptions and the pulse-matched
    readings; writes c_comparisons.csv, c_wirings.csv and verdict.json."""
    if a0 is None:
        raise RuntimeError("part C needs the definition of the model's pre-odor level (--a0 e90|settled|blank): the blank "
                           "identity fails with E90's window (the KC layer's start-up), so the definition waits on a decision")
    b1 = load_b1_all(wirings)
    missing = b1.attrs.get('missing') or []
    if missing and not allow_missing:
        raise RuntimeError(f'B1 runs missing: {missing}')
    b1_sources = dict(sources=b1.attrs.get('sources'), e101=b1.attrs.get('e101'))
    b1 = with_a0(b1, a0)
    recs = load_recordings()
    resp = stimulus_responses(recs)
    flies, per = {}, {}
    for k, (rd, sub) in READING_SETS.items():
        flies[k], per[k] = flies_ratios(resp, rd, sub)
    short = sorted(resp.loc[np.isclose(resp['pulse_s'], 0.5), 'fly'].unique())
    pulse = {}
    pc = per['carrying']
    for rid in ('V1_paa-4', 'V1_paa-2.954', 'V2'):
        d = pc[(pc['ratio_id'] == rid) & pc['entered']]
        pulse[rid] = pulse_permutation(d.groupby('fly')['log_ratio'].mean(), short)
    k = kc_rates('z1')
    verdict, comps, rows = verdict_core(flies, b1, pulse, solvent=dict(ctrl=float(k.loc['pfo', 'ctrl']), tnt=float(k.loc['pfo', 'tnt'])))
    drift = a1_tables(recs, resp)['drift_summary']
    verdict['V2_descriptions'] = v2_descriptions(pc, recs, drift)
    rd = CARRYING[0]
    verdict['pulse_matched'] = pulse_matched(b1, flies_ratios(resp, rd, 'pulse_0.5s')[0], flies_ratios(resp, rd, 'pulse_1s')[0])
    out = out_dir()
    comps.to_csv(out / 'c_comparisons.csv', index=False)
    rows.to_csv(out / 'c_wirings.csv', index=False)
    verdict.update(a0=a0, b1_sources=b1_sources, reading_sets={k: dict(reading=reading_label(v[0]), subset=v[1]) for k, v in READING_SETS.items()},
                   wirings=sorted(set(b1['wiring'])), missing_runs=missing, time=time.strftime('%Y-%m-%d %H:%M'),
                   host=_host(), commits=dict(H.code_commits()))
    _write_json(out / 'verdict.json', verdict)
    return verdict


# ---------------------------------------------------------------------------
# part D: Yang's 2021 recordings outside the 12 (the entry's (C) as amended 2026-11-15 17:50)

D_FOLDER = ('data', 'yang_apl_recordings', '2021_jul_sep')
D_SET = 'A2021'
D1_SHUFFLED = ('20210714_2FRfun2', '20210714_2FRfun5', '20210716_1FLfun2', '20210716_3FLfun4', '20210823_4FLfun5',
               '20210826_1FLfun5', '20210903_1FLfun6')
D2_PLANES: Dict[int, Tuple[str, ...]] = {
    15: ('20210714_2FLfun1', '20210714_2FRfun4', '20210805_3FLfun3', '20210823_4FLfun2', '20210823_4FRfun1',
         '20210826_1FLfun3', '20210827_3FRfun1', '20210905_1FLfun1', '20210905_1FRfun1'),
    18: ('20210716_3FLfun3', '20210804_2FLfun3', '20210823_4FLfun4', '20210827_3FRfun2'),
    21: ('20210823_4FLfun3', '20210826_1FLfun4', '20210903_1FLfun2', '20210903_1FLfun5', '20210903_1FLfun7',
         '20210903_1FRfun2'),
    24: ('20210903_1FLfun3', '20210903_1FRfun3')}
D2_ASCENDING = tuple(x for v in D2_PLANES.values() for x in v)
D3_REPEATS_ONLY = ('20210805_3FLfun2', '20210903_1FLfun1')
D_RECORDINGS = D1_SHUFFLED + D2_ASCENDING + D3_REPEATS_ONLY
YANG_LISTS: Dict[int, Tuple[str, ...]] = {
    15: ('20210905_1FLfun1', '20210905_1FRfun1', '20210826_1FLfun3', '20210823_4FRfun1', '20210714_2FRfun4'),
    21: ('20210903_1FLfun2', '20210903_1FRfun2', '20210826_1FLfun4', '20210823_4FLfun3', '20210716_3FLfun4',
         '20210714_2FRfun2')}
D2_NOT_INTERLEAVED = ('20210714_2FLfun1', '20210714_2FRfun4', '20210716_3FLfun3')   # one odor's series in consecutive blocks
D2_V2_ONLY = ('20210823_4FLfun2', '20210823_4FLfun4')                             # no ms at 10^-2
D_EXTRA_OUT: Dict[str, Tuple[int, ...]] = {'20210903_1FLfun6': (1,)}             # Yang's note: not vented (blocks from 1)
D_BLOCK_OUT: Dict[str, Tuple[int, ...]] = {                                      # the yaml's block_out, blocks from 1
    '20210714_2FLfun1': (1, 12, 13), '20210805_3FLfun2': (8, 9, 10), '20210716_3FLfun3': (1,), '20210903_1FRfun3': (1,),
    '20210903_1FLfun5': (6,), '20210903_1FLfun7': (1,)}
D1_LAGS: Dict[str, float] = {'20210714_2FRfun5': -0.55, '20210716_1FLfun2': -0.55, '20210823_4FLfun5': -0.55,
                             '20210714_2FRfun2': -0.97, '20210716_3FLfun4': 6.91, '20210826_1FLfun5': 6.91,
                             '20210903_1FLfun6': 6.91}
LAG_TOL = 0.005
ORDER_RHO_MIN = 0.9      # the 12's ascending design swaps one pa pair in three recordings (rho 0.943)
D_OUTSIDE_THE_8 = ('2021-07-16_fly1', '2021-07-16_fly3')                       # the two flies outside the 12's 8
D3_SAME_PLANE = (('20210714_2FRfun5', '20210714_2FRfun1'), ('20210823_4FLfun5', '20210823_4FLfun3'),
                 ('20210903_1FLfun6', '20210903_1FLfun3'), ('20210827_3FRfun4', '20210827_3FRfun2'))
D3_ADJACENT = (('20210826_1FLfun5', '20210826_1FLfun4'), ('20210716_3FLfun4', '20210716_3FLfun3'),
               ('20210714_2FRfun2', '20210714_2FRfun1'))
D3_REPEATS: Tuple[Tuple[str, Stim], ...] = (
    ('20210805_3FLfun2', ('pa', -6.771)), ('20210823_4FLfun2', ('pa', -6.771)), ('20210823_4FLfun4', ('pa', -6.771)),
    ('20210805_3FRfun1', ('pa', -6.771)), ('20210823_4FRfun2', ('pa', -6.771)),
    ('20210903_1FLfun1', ('paa', -5.34)), ('20210903_1FLfun1', ('ms', -6.294)))
D3_RERUNS = (('20210804_2FLfun3', '20210804_2FLfun2'), ('20210823_4FLfun4', '20210823_4FLfun1'),
             ('20210903_1FLfun5', '20210903_1FLfun2'), ('20210903_1FLfun7', '20210903_1FLfun2'))
D_V_RATIOS = tuple(r for r in RATIOS if r[0] in ('V1_paa-4', 'V1_paa-2.954', 'V2'))
D_BASELINES = ('order_matched', 'first_in_block')
SOLVENT = 'pfo'


def part_c_done() -> bool:
    """Whether part C's verdict is written (part D runs only after it)."""
    return (out_dir() / 'verdict.json').exists()


def twelve() -> List[str]:
    """The 12 recordings' stems (E79's recordings.csv)."""
    return list(recording_info().index)


def fly_id(date, fly) -> str:
    """E79's fly id, '<YYYY-MM-DD>_fly<n>', from a date written YYYYMMDD or YYYY-MM-DD."""
    d = str(date).replace('-', '')
    return f'{d[:4]}-{d[4:6]}-{d[6:8]}_fly{int(float(fly))}'


def d_rows(stems: Sequence[str]) -> pd.DataFrame:
    """viki's recordings.csv rows of `stems` (E100's `read_recordings`, every set; the row whose date is the stem's, so the
    one stem written twice resolves), indexed by stem, with 'fly' (E79's form), 'side' and 'role' (D1, D2, D3 or the 12).
    Raises on a stem outside part D's lists and the 12, on a missing or repeated row and on a row that is not processed
    (`use_until_further_notice`)."""
    from experiments import E100_apl_weak_odors as E100
    allowed = set(D_RECORDINGS) | set(twelve())
    bad = [x for x in stems if x not in allowed]
    if bad:
        raise PermissionError(f'part D reads only its lists and the 12: {bad}')
    r = E100.read_recordings(sets=None)
    out = []
    for st in stems:
        g = r[(r['stem'] == st) & (r['date'].astype(str).str.replace('-', '') == st[:8])]
        if len(g) != 1:
            raise RuntimeError(f'{st}: {len(g)} rows in recordings.csv')
        row = g.iloc[0].to_dict()
        if not str(row.get('use_until_further_notice', '')).startswith('processed'):
            raise PermissionError(f"{st}: not processed ({row.get('use_until_further_notice')})")
        row['fly_id'] = fly_id(row['date'], row['fly'])
        row['side'] = str(row['hemisphere'])
        row['role'] = 'D1' if st in D1_SHUFFLED else 'D2' if st in D2_ASCENDING else 'D3' if st in D3_REPEATS_ONLY else '12'
        out.append(row)
    return pd.DataFrame(out).set_index('stem', drop=False)


def d_yaml_path(row) -> Path:
    from experiments import E100_apl_weak_odors as E100
    p = E100.recording_paths(row)['yaml']
    if p is None:
        raise FileNotFoundError(f"{row['stem']}: no yaml")
    return Path(p)


def yaml_block_out(path: Path) -> Tuple[int, ...]:
    """The yaml's `block_out` (blocks numbered from 1), empty if none."""
    import yaml
    y = yaml.safe_load(Path(path).read_text())
    return tuple(int(x) for x in (y.get('block_out') or []))


def d_stim_table(path: Path, stem: str, extra_out: Sequence[int] = ()) -> pd.DataFrame:
    """One recording's presentations from its yaml (E100's `yaml_stimuli`; blocks from 0 there, 'block1' from 1 here):
    'odor' with the solvent as SOLVENT, 'out' (the block in the yaml's `block_out` or in `extra_out`, from 1),
    'pfo_block' (every presentation of the block is the solvent), 'first' (the stimulus's first block among the kept
    ones)."""
    from experiments import E100_apl_weak_odors as E100
    t = E100.yaml_stimuli(path, stem)
    t['odor'] = [SOLVENT if o == E100.SOLVENT else o for o in t['odor']]
    t['block1'] = t['block'] + 1
    out = set(yaml_block_out(path)) | set(int(x) for x in extra_out)
    t['out'] = t['block1'].isin(out)
    t['pfo_block'] = t.groupby('block')['odor'].transform(lambda x: bool((x == SOLVENT).all())).astype(bool)
    kept = t[~t['out'] & (t['odor'] != SOLVENT)]
    first = kept.groupby(['odor', 'log10_conc'])['block1'].min()
    t['first'] = [bool(o != SOLVENT and not ou and first.get((o, c)) == b)
                  for o, c, b, ou in zip(t['odor'], t['log10_conc'], t['block1'], t['out'])]
    return t


def d_stimulus_blocks(table: pd.DataFrame, which: str = 'first') -> pd.DataFrame:
    """The (odor, conc, block1) of every stimulus read: its first kept block ('first'), or every kept block ('all')."""
    k = table[~table['out'] & (table['odor'] != SOLVENT)]
    if which == 'first':
        k = k[k['first']]
    return k.drop_duplicates(['odor', 'log10_conc', 'block1'])[['odor', 'log10_conc', 'block1']]


def d_responses_from_means(means: Sequence[float], table: pd.DataFrame, baseline: str, which: str = 'first') -> pd.DataFrame:
    """Each stimulus's response from per-presentation calyx means (`means`, in the yaml's presentation order) and the
    stimulus table (`d_stim_table`): 'order_matched' each presentation minus the mean over the kept solvent blocks of
    their presentation at the same position; 'first_in_block' each presentation minus the mean of the kept blocks'
    position-1 solvent presentations; the response the mean over the stimulus's presentations in its block (two in a
    block of three, three in a block of five), nothing further subtracted. `which` 'first' (the first kept block) or 'all'
    (every kept block, part D3's repeats). Rows: odor, conc, block1, response, rep1 to rep3, n_reps."""
    m = np.asarray(means, dtype=float)
    if len(m) != len(table):
        raise ValueError(f'{len(m)} presentation values for {len(table)} presentations')
    t = table.reset_index(drop=True)
    kept = ~t['out']
    if baseline == 'order_matched':
        pb = t[kept & t['pfo_block']]
        if pb.empty:
            raise RuntimeError('no kept solvent block for the order-matched baseline')
        base = pb.groupby('position').apply(lambda g: float(m[g.index].mean())).to_dict()
    elif baseline == 'first_in_block':
        f = t[kept & (t['position'] == 1) & (t['odor'] == SOLVENT)]
        base = float(m[f.index].mean())
    else:
        raise ValueError(f'unknown baseline {baseline!r}')
    rows = []
    for _, sb in d_stimulus_blocks(t, which).iterrows():
        g = t[(t['block1'] == sb['block1']) & (t['odor'] == sb['odor']) & (t['log10_conc'] == sb['log10_conc'])]
        reps = [float(m[i] - (base[int(pos)] if isinstance(base, dict) else base)) for i, pos in zip(g.index, g['position'])]
        rows.append(dict(odor=sb['odor'], conc=float(sb['log10_conc']), block1=int(sb['block1']), response=float(np.mean(reps)),
                         n_reps=len(reps), **{f'rep{j + 1}': (reps[j] if j < len(reps) else np.nan) for j in range(3)}))
    return pd.DataFrame(rows)


def d_frame_means(row, source: str = 'maps') -> np.ndarray:
    """One recording's per-presentation calyx means: 'maps' its peak dF/F maps inside the calyx outline (E79's mask,
    E100's `calyx_responses`), 'yang' Yang's ROI table (roi_peak_dff in presentation order, its pins checked against the
    yaml's). Refuses a recording of part D's lists before part C's verdict exists."""
    from experiments import E100_apl_weak_odors as E100
    from data_prep import apl_calyx_maps as ACM
    import tifffile
    if row['stem'] in D_RECORDINGS and not part_c_done():
        raise PermissionError(f"{row['stem']}: part D opens no map or ROI table before part C's verdict exists")
    paths = E100.recording_paths(row)
    if source == 'maps':
        stack = tifffile.imread(paths['peak_maps']).astype(float)
        mask = ACM.roi_mask(ACM.read_imagej_polygon(paths['roi']), stack.shape[1:])
        return E100.calyx_responses(stack, mask)
    if source == 'yang':
        t = pd.read_csv(paths['roi_table'])
        y = E100.yaml_stimuli(d_yaml_path(row), row['stem'])
        if list(t['pin_n'].astype(int)) != list(y['pin'].astype(int)):
            raise RuntimeError(f"{row['stem']}: Yang's ROI table's pins are not the yaml's")
        return t['roi_peak_dff'].to_numpy(float)
    raise ValueError(f'unknown source {source!r}')


def d_resp(stems: Sequence[str], source: str = 'maps') -> pd.DataFrame:
    """Part D's responses in E98's `stimulus_responses` format (baselines D_BASELINES, nothing subtracted beyond the
    baseline, 'blank' NaN), for `stems` (part D's lists and the 12), from the maps or Yang's tables."""
    rows = d_rows(stems)
    out = []
    for st, row in rows.iterrows():
        tab = d_stim_table(d_yaml_path(row), st, D_EXTRA_OUT.get(st, ()))
        m = d_frame_means(row, source)
        for bl in D_BASELINES:
            r = d_responses_from_means(m, tab, bl)
            for x in r.to_dict('records'):
                out.append(dict(recording=st, fly=row['fly_id'], side=row['side'], pulse_s=float(row['pulse_s']),
                                plane_um=row['plane_depth_um'], role=row['role'], odor=x['odor'], conc=x['conc'],
                                stimulus=f"{x['odor']} {x['conc']:g}", block=x['block1'] - 1, baseline=bl,
                                blank_subtracted=False, blank=np.nan, response=x['response'], rep1=x['rep1'], rep2=x['rep2'],
                                rep3=x['rep3'], n_reps=x['n_reps'], source=source))
    return pd.DataFrame(out)


def interp_weights(concs: Sequence[float], c: float, design: Optional[Sequence[float]] = None) -> Dict[float, float]:
    """The weights `value_at` puts on a recording's imaged concentrations to read concentration `c` (1 on an imaged or
    substituted one, the two neighbours' linear weights on log10 concentration when interpolated); empty when it cannot
    be read; `design` the series' concentrations, whose range and members `value_at` applies)."""
    x = sorted(float(v) for v in concs)
    rng = (min(design), max(design)) if design is not None else None
    _, how = value_at(x, [0.0] * len(x), c, rng, design)
    if how == 'imaged':
        return {next(v for v in x if abs(v - c) <= CONC_TOL): 1.0}
    if how == 'interpolated':
        hi = next(v for v in x if v > c)
        lo = max(v for v in x if v < c)
        w = (hi - c) / (hi - lo)
        return {lo: w, hi: 1.0 - w}
    if how == 'substituted':
        return {concentration_match(c, [x[0], x[-1]]): 1.0}
    return {}


def block_lag(table: pd.DataFrame, num: Stim, den: Stim, design: Optional[Dict[str, Sequence[float]]] = None) -> float:
    """A recording's block lag for a ratio: the interpolation-weighted block (from 1) of the numerator's concentration
    minus the denominator's, each over the first kept blocks of its imaged neighbours (NaN if either cannot be read)."""
    sb = d_stimulus_blocks(table)
    out = []
    for o, c in (num, den):
        g = sb[sb['odor'] == o].set_index('log10_conc')['block1']
        w = interp_weights(g.index, c, (design or {}).get(o))
        if not w:
            return float('nan')
        out.append(sum(wt * float(g[k]) for k, wt in w.items()))
    return float(out[0] - out[1])


def order_class(table: pd.DataFrame, odors: Sequence[str] = ('pa', 'paa', 'ms')) -> str:
    """'ascending' if every odor's concentration rises with block over its stimuli's first kept blocks (Spearman's rho
    of concentration on block at least ORDER_RHO_MIN for each odor), else 'shuffled'."""
    from scipy.stats import spearmanr
    sb = d_stimulus_blocks(table)
    for o in odors:
        g = sb[sb['odor'] == o].sort_values('block1')
        if len(g) > 2 and not spearmanr(g['block1'], g['log10_conc'])[0] >= ORDER_RHO_MIN - 1e-12:
            return 'shuffled'
    return 'ascending'


def twelve_design() -> Dict[str, List[float]]:
    """The 12's design: each odor's imaged concentrations over the 12 recordings (part A's responses); its range is the
    design range part D passes to `value_at`."""
    a1 = pd.read_csv(out_dir() / 'a1_responses.csv')
    return {o: sorted(float(c) for c in g.unique()) for o, g in a1.groupby('odor')['conc']}


def eligibility(table: pd.DataFrame, design: Dict[str, Sequence[float]]) -> dict:
    """From the stimulus table alone: the order class, V1's and V2's availability (every stimulus they read can be read
    from the kept blocks, `value_at` with the 12's design ranges) and V2's block lag."""
    sb = d_stimulus_blocks(table)

    def ok(st: Stim) -> bool:
        g = sb[sb['odor'] == st[0]]
        dc = design.get(st[0])
        return bool(np.isfinite(value_at(g['log10_conc'], [1.0] * len(g), st[1], (min(dc), max(dc)), dc)[0]))
    v1 = all(ok(x) for x in (('paa', -4.0), ('paa', -2.954), ('ms', -3.0)))
    v2 = all(ok(x) for x in (('paa', -5.0), ('paa', -4.0)))
    return dict(order=order_class(table), v1=v1, v2=v2,
                v2_lag=block_lag(table, ('paa', -5.0), ('paa', -4.0), design) if v2 else np.nan)


def d0_expected(stem: str) -> dict:
    """The eligibility the amendment's lists give a recording (D3's: order and ratios not fixed, NaN)."""
    if stem in D1_SHUFFLED:
        return dict(order='shuffled', v1=True, v2=True, v2_lag=D1_LAGS[stem])
    if stem in D2_ASCENDING:
        return dict(order='ascending', v1=stem not in D2_V2_ONLY, v2=True, v2_lag=np.nan)
    if stem == '20210805_3FLfun2':
        return dict(order=None, v1=None, v2=False, v2_lag=np.nan)
    return dict(order=None, v1=None, v2=None, v2_lag=np.nan)


def d0_eligibility() -> pd.DataFrame:
    """D0's eligibility table, from the yaml files alone (no map or ROI table opened): per recording of part D's lists and
    of the 12, its role, plane, pulse, fly, the yaml's `block_out` and part D's extra blocks out, the number of
    presentations (and viki's frames_equal_presentations), the order class, V1's and V2's availability and V2's block
    lag, beside the amendment's (`d0_expected`) and whether they match."""
    ranges = twelve_design()
    stems = list(D_RECORDINGS) + twelve()
    rows = d_rows(stems)
    out = []
    for st, row in rows.iterrows():
        path = d_yaml_path(row)
        tab = d_stim_table(path, st, D_EXTRA_OUT.get(st, ()))
        el = eligibility(tab, ranges)
        ex = d0_expected(st) if st in D_RECORDINGS else dict(order=None, v1=None, v2=None, v2_lag=np.nan)
        bo = yaml_block_out(path)
        rec = dict(stem=st, role=row['role'], fly=row['fly_id'], side=row['side'], plane_um=row['plane_depth_um'],
                   pulse_s=float(row['pulse_s']), n_presentations=int(len(tab)),
                   frames_equal_presentations=row.get('frames_equal_presentations'), nstim=int(tab['position'].max()),
                   pfo_blocks=';'.join(str(b) for b in sorted(tab.loc[tab['pfo_block'], 'block1'].unique())),
                   block_out=';'.join(str(b) for b in bo), extra_out=';'.join(str(b) for b in D_EXTRA_OUT.get(st, ())),
                   block_out_expected=';'.join(str(b) for b in D_BLOCK_OUT.get(st, ())),
                   outside_the_8=row['fly_id'] in D_OUTSIDE_THE_8, **el,
                   **{f'{k}_expected': v for k, v in ex.items()})
        m = [rec['block_out'] == rec['block_out_expected'] if st in D_RECORDINGS else True]
        for k in ('order', 'v1', 'v2'):
            if ex[k] is not None:
                m.append(rec[k] == ex[k])
        if np.isfinite(ex['v2_lag']):
            m.append(bool(abs(rec['v2_lag'] - ex['v2_lag']) <= LAG_TOL))
        rec['matches'] = bool(all(m))
        out.append(rec)
    return pd.DataFrame(out)


def d0_reproduce_12() -> dict:
    """D0: the generalized reader pointed at the 12 through their yaml (not frames.csv) reproduces part A's
    a1_responses.csv under both baselines (nothing subtracted), within IDENTITY_TOL, on the same recordings and stimuli."""
    a1 = pd.read_csv(out_dir() / 'a1_responses.csv')
    a1 = a1[~a1['blank_subtracted'].astype(bool) & a1['baseline'].isin(D_BASELINES)]
    mine = d_resp(twelve(), 'maps')
    key = ['recording', 'stimulus', 'baseline']
    a = a1.set_index(key)[['response', 'rep1', 'rep2']].sort_index()
    b = mine.set_index(key)[['response', 'rep1', 'rep2']].sort_index()
    same = a.index.equals(b.index)
    d = float((a - b).abs().to_numpy().max()) if same else float('nan')
    return dict(index_equal=bool(same), n=int(len(a)), max_abs_diff=d, exact=bool(same and d <= IDENTITY_TOL))


def part_d0() -> dict:
    """D0 (yaml only for part D's recordings; the 12's maps for the reader's check): d0_eligibility.csv, d0.json. Raises
    if the eligibility differs from the amendment's lists or the reader does not reproduce part A."""
    out = out_dir()
    el = d0_eligibility()
    el.to_csv(out / 'd0_eligibility.csv', index=False)
    rep = d0_reproduce_12()
    res = dict(eligibility_matches=bool(el['matches'].all()), mismatches=el.loc[~el['matches']].to_dict('records'),
               reproduce_12=rep, n_recordings=int(len(el)), time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(out / 'd0.json', res)
    if not res['eligibility_matches']:
        raise RuntimeError(f"D0: the eligibility from the yaml differs from the amendment's lists: {res['mismatches']}")
    if not rep['exact']:
        raise RuntimeError(f"D0: the generalized reader does not reproduce part A on the 12: {rep}")
    return res


def _cat(on: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> Tuple[str, pd.DataFrame]:
    wc = stimulus_classes(on, num, den, lo, hi)
    return comparison_verdict(wc['cls']), wc


def _ratio_row(name: str, fr: dict, cat: str, wc: pd.DataFrame, **extra) -> dict:
    return dict(reading=name, ratio=fr['ratio'], lo90=fr['lo90'], hi90=fr['hi90'], n_flies=fr['n_flies'],
                n_recordings=fr['n_recordings'], sd_log2_flies=fr['sd_log2_flies'], category=cat,
                **{f'n_{k}': int((wc['cls'] == k).sum()) for k in CATEGORIES + ('mixed',)}, **extra)


def d1_v2r(resp: pd.DataFrame, ranges: Dict[str, Tuple[float, float]], on: pd.DataFrame,
           extra_recordings: Sequence[str] = ('20210827_3FRfun4',)) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """D1: V2 in the shuffled-order recordings (V2R) under the carrying reading, the first-in-block baseline and each
    reading with one fly left out; its category is held when all of them give it and the 0.5 s and 1 s flies do not
    disagree (`pulse_permutation`, the 15 splits of 6 flies into 2 and 4); beside, V2R with `extra_recordings` (the 12's
    shuffled one) and V1's two ratios in the same recordings. Returns (one row per reading, per wiring classes under the
    carrying reading with each wiring's distance from the interval's nearer end, the result)."""
    rows, cats = [], {}
    rec1 = [r for r in D1_SHUFFLED if r in set(resp['recording'])]
    per = {}
    for name, rd in (('carrying', CARRYING[0]), ('first_in_block', PARTNERS['first_in_block'][0])):
        f, pr = flies_ratios(resp, rd, 'all', D_V_RATIOS, None, rec1, ranges)
        per[name] = pr
        fr = f.set_index('ratio_id').loc['V2']
        cat, wc = _cat(on, V2_NUM, V2_DEN, fr['lo90'], fr['hi90'])
        cats[name] = cat
        rows.append(_ratio_row(name, fr.to_dict(), cat, wc, decides=True))
        if name == 'carrying':
            classes = wc.assign(category=cat)
            for rid in ('V1_paa-4', 'V1_paa-2.954'):
                g = f.set_index('ratio_id').loc[rid]
                rows.append(dict(reading=f'carrying, {rid} (beside)', ratio=g['ratio'], lo90=g['lo90'], hi90=g['hi90'],
                                 n_flies=g['n_flies'], n_recordings=g['n_recordings'], decides=False))
    pc = per['carrying']
    v2 = pc[(pc['ratio_id'] == 'V2') & pc['entered']]
    for fly in sorted(v2['fly'].unique()):
        fr = fly_ratio(v2[v2['fly'] != fly])
        cat, wc = _cat(on, V2_NUM, V2_DEN, fr['lo90'], fr['hi90'])
        cats[f'without {fly}'] = cat
        rows.append(_ratio_row(f'carrying, without {fly}', fr, cat, wc, decides=True))
    flies = sorted(v2['fly'].unique())
    short = sorted(resp.loc[resp['recording'].isin(rec1) & np.isclose(resp['pulse_s'], 0.5), 'fly'].unique())
    pulse = pulse_permutation(v2.groupby('fly')['log_ratio'].mean(), short)
    c0 = cats['carrying']
    held = c0 if (c0 in CATEGORIES and all(v == c0 for v in cats.values()) and not pulse['disagree']) else 'not decided'
    extra = [r for r in extra_recordings if r in set(resp['recording'])]
    if extra:
        f, _ = flies_ratios(resp, CARRYING[0], 'all', D_V_RATIOS, None, rec1 + extra, ranges)
        fr = f.set_index('ratio_id').loc['V2']
        cat, wc = _cat(on, V2_NUM, V2_DEN, fr['lo90'], fr['hi90'])
        rows.append(_ratio_row(f"carrying, with {', '.join(extra)} (beside)", fr.to_dict(), cat, wc, decides=False))
    res = dict(V2R=held, categories=cats, pulse=pulse, flies=flies, short_flies=short,
               outside_the_8=[f for f in flies if f in D_OUTSIDE_THE_8])
    return pd.DataFrame(rows), classes, res


def d1_lag_slope(resp: pd.DataFrame, ranges: Dict[str, Tuple[float, float]], lags: pd.Series) -> pd.DataFrame:
    """D1 beside: over every recording (part D's and the 12) of the flies that have both orders, the slope of V2's log
    ratio (carrying reading, one recording) on its block lag: pooled, and within flies (each fly's mean removed)."""
    _, pr = flies_ratios(resp, CARRYING[0], 'all', D_V_RATIOS, None, sorted(set(resp['recording'])), ranges)
    v = pr[(pr['ratio_id'] == 'V2') & pr['entered']].copy()
    v['lag'] = v['recording'].map(lags)
    v['order'] = np.where(v['recording'].isin(D1_SHUFFLED) | (v['recording'] == '20210827_3FRfun4'), 'shuffled', 'ascending')
    both = v.groupby('fly')['order'].nunique()
    v = v[v['fly'].isin(both.index[both == 2]) & np.isfinite(v['lag'])]
    rows = []
    if len(v) >= 3:
        pooled = np.polyfit(v['lag'], v['log_ratio'], 1)[0]
        dv = v.assign(x=v['lag'] - v.groupby('fly')['lag'].transform('mean'),
                      y=v['log_ratio'] - v.groupby('fly')['log_ratio'].transform('mean'))
        within = float((dv['x'] * dv['y']).sum() / (dv['x'] ** 2).sum()) if (dv['x'] ** 2).sum() > 0 else np.nan
        rows.append(dict(n_recordings=int(len(v)), n_flies=int(v['fly'].nunique()), slope_pooled=float(pooled),
                         slope_within_flies=within, flies=';'.join(sorted(v['fly'].unique())),
                         note='log of V2 (natural log) per block of lag'))
    return pd.concat([pd.DataFrame(rows), v[['recording', 'fly', 'order', 'lag', 'log_ratio']]], ignore_index=True)


def d_twelve_categories(resp12: pd.DataFrame, ranges, on: pd.DataFrame) -> dict:
    """The 12's categories under the carrying reading (V1 under each estimate of ms's input on the hull; V2), as part C
    reads them."""
    f = flies_ratios(resp12, CARRYING[0], 'all', D_V_RATIOS, None, twelve(), ranges)[0]
    return _categories(f, on)


def _categories(f: pd.DataFrame, on: pd.DataFrame) -> dict:
    f = f.set_index('ratio_id')
    out = {}
    iv = _intervals(f.reset_index(), ('V1_paa-4', 'V1_paa-2.954'))
    for est, lab in V1_ESTIMATES.items():
        out[f'V1 {est}'] = _cat(on, V1_NUM, lab, *iv)[0] if iv is not None else 'not run'
    out['V2'] = _cat(on, V2_NUM, V2_DEN, f.loc['V2', 'lo90'], f.loc['V2', 'hi90'])[0] if np.isfinite(f.loc['V2', 'lo90']) \
        else 'not run'
    return out


def d2_sets(resp: pd.DataFrame) -> Dict[str, List[str]]:
    have = set(resp['recording'])
    sets = {'all21': [x for x in D2_ASCENDING if x in have]}
    for pl, xs in D2_PLANES.items():
        sets[f'{pl} um'] = [x for x in xs if x in have]
    yl = {x for v in YANG_LISTS.values() for x in v}
    sets['yang_lists'] = [x for x in D2_ASCENDING if x in yl and x in have]
    sets['interleaved'] = [x for x in D2_ASCENDING if x not in D2_NOT_INTERLEAVED and x in have]
    return sets


def d2_replication(resp: pd.DataFrame, ranges, on: pd.DataFrame, cats12: dict) -> pd.DataFrame:
    """D2: V1's two ratios and V2 over each of `d2_sets`, the model's category against each set's interval (V1 under
    each estimate of ms's input, on the hull; V2), and whether it replicates (equals the 12's category under the carrying
    reading: V1 with the megamat estimate, V2 with validation2's paa at 10^-5). Decides nothing."""
    rows = []
    for name, recs in d2_sets(resp).items():
        if not recs:
            continue
        f = flies_ratios(resp, CARRYING[0], 'all', D_V_RATIOS, None, recs, ranges)[0]
        cats = _categories(f, on)
        fi = f.set_index('ratio_id')
        iv = _intervals(f, ('V1_paa-4', 'V1_paa-2.954'))
        for k, c in cats.items():
            rid = 'V2' if k == 'V2' else 'V1'
            rows.append(dict(set=name, n_recordings=len(recs), comparison=k, category=c, category_12=cats12.get(k),
                             replicates=(bool(c == cats12.get(k)) if k in ('V2', 'V1 megamat/megamat (the panel)') else None),
                             flies_lo90=(fi.loc['V2', 'lo90'] if rid == 'V2' else (iv[0] if iv else np.nan)),
                             flies_hi90=(fi.loc['V2', 'hi90'] if rid == 'V2' else (iv[1] if iv else np.nan)),
                             V1_paa4=fi.loc['V1_paa-4', 'ratio'], V1_paa2954=fi.loc['V1_paa-2.954', 'ratio'],
                             V2=fi.loc['V2', 'ratio'], n_flies_V2=int(fi.loc['V2', 'n_flies']),
                             outside_the_8=';'.join(sorted({f for f in resp.loc[resp['recording'].isin(recs), 'fly']}
                                                           & set(D_OUTSIDE_THE_8)))))
    return pd.DataFrame(rows)


def _per_recording_logs(resp: pd.DataFrame, ranges, recs: Sequence[str], baseline=CARRYING[0]) -> pd.DataFrame:
    _, pr = flies_ratios(resp, baseline, 'all', D_V_RATIOS, None, list(recs), ranges)
    return pr.pivot_table(index='recording', columns='ratio_id', values='log_ratio')


def d2_plane_pairs(resp: pd.DataFrame, ranges) -> pd.DataFrame:
    """D2: per calyx (fly and side), the log ratio at 15, 21 or 24 um minus the log ratio at 18 um (the 12's recording of
    that calyx), per ratio; then per plane and ratio the mean over calyces and the count of positive differences."""
    meta = resp.drop_duplicates('recording').set_index('recording')
    logs = _per_recording_logs(resp, ranges, sorted(set(resp['recording'])))
    t12 = set(twelve())
    rows = []
    for rec in D2_ASCENDING:
        if rec not in logs.index:
            continue
        pl = next(k for k, v in D2_PLANES.items() if rec in v)
        if pl == 18:
            continue
        ref = [r for r in t12 if r in logs.index and meta.loc[r, 'fly'] == meta.loc[rec, 'fly']
               and meta.loc[r, 'side'] == meta.loc[rec, 'side']]
        if len(ref) != 1:
            continue
        for rid in ('V1_paa-4', 'V1_paa-2.954', 'V2'):
            d = logs.loc[rec, rid] - logs.loc[ref[0], rid] if rid in logs.columns else np.nan
            rows.append(dict(recording=rec, reference_18um=ref[0], fly=meta.loc[rec, 'fly'], side=meta.loc[rec, 'side'],
                             plane_um=pl, ratio_id=rid, log_diff=float(d)))
    per = pd.DataFrame(rows)
    summ = []
    if len(per):
        for (pl, rid), g in per.dropna(subset=['log_diff']).groupby(['plane_um', 'ratio_id']):
            c = g.groupby(['fly', 'side'])['log_diff'].mean()
            summ.append(dict(plane_um=pl, ratio_id=rid, n_calyces=int(len(c)), mean_log_diff=float(c.mean()),
                             n_positive=int((c > 0).sum()), recording='(summary)'))
    return pd.concat([per, pd.DataFrame(summ)], ignore_index=True)


def d3_tables(resp: pd.DataFrame, ranges, lags: pd.Series, means: Dict[str, np.ndarray],
              tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """D3's descriptions: shuffled against ascending in one calyx (d3_pairs: the difference of log V2 and of the block
    lags), stimuli repeated within a recording (d3_repeats: the later block's response over the earlier's, carrying
    baseline over every kept solvent block), repeat runs of one calyx and plane (d3_reruns: the later run's ratio over
    the earlier's)."""
    logs = _per_recording_logs(resp, ranges, sorted(set(resp['recording'])))
    pairs = []
    for kind, ps in (('same plane', D3_SAME_PLANE), ('adjacent planes', D3_ADJACENT)):
        for a, b in ps:
            pairs.append(dict(kind=kind, shuffled=a, ascending=b,
                              log_v2_diff=float(logs.loc[a, 'V2'] - logs.loc[b, 'V2']) if a in logs.index and b in logs.index
                              else np.nan, lag_diff=float(lags.get(a, np.nan) - lags.get(b, np.nan))))
    reps = []
    for rec, (o, c) in D3_REPEATS:
        if rec not in means:
            continue
        r = d_responses_from_means(means[rec], tables[rec], CARRYING[0][0], which='all')
        g = r[(r['odor'] == o) & np.isclose(r['conc'], c)].sort_values('block1')
        if len(g) == 2:
            reps.append(dict(recording=rec, stimulus=f'{o} {c:g}', earlier_block=int(g['block1'].iloc[0]),
                             later_block=int(g['block1'].iloc[1]), earlier=float(g['response'].iloc[0]),
                             later=float(g['response'].iloc[1]),
                             later_over_earlier=float(g['response'].iloc[1] / g['response'].iloc[0])))
    reruns = []
    for later, earlier in D3_RERUNS:
        for rid in ('V1_paa-4', 'V1_paa-2.954', 'V2'):
            ok = later in logs.index and earlier in logs.index and rid in logs.columns
            reruns.append(dict(later=later, earlier=earlier, ratio_id=rid,
                               later_over_earlier=float(np.exp(logs.loc[later, rid] - logs.loc[earlier, rid])) if ok else np.nan))
    return dict(d3_pairs=pd.DataFrame(pairs), d3_repeats=pd.DataFrame(reps), d3_reruns=pd.DataFrame(reruns))


def d_yang_table_ratios(resp_yang: pd.DataFrame, ranges) -> pd.DataFrame:
    """Every ratio of part D's sets also from Yang's ROI tables by the same arithmetic (a description): the 12 first, then
    V2R's recordings and D2's sets, under both baselines."""
    rows = []
    groups = [('the 12', twelve()), ('D1 (V2R)', list(D1_SHUFFLED))] + [(f'D2 {k}', v) for k, v in d2_sets(resp_yang).items()]
    for name, recs in groups:
        recs = [r for r in recs if r in set(resp_yang['recording'])]
        if not recs:
            continue
        for bl in D_BASELINES:
            f = flies_ratios(resp_yang, (bl, False), 'all', D_V_RATIOS, None, recs, ranges)[0]
            for _, r in f.iterrows():
                rows.append(dict(group=name, baseline=bl, ratio_id=r['ratio_id'], ratio=r['ratio'], lo90=r['lo90'],
                                 hi90=r['hi90'], n_flies=r['n_flies'], n_recordings=r['n_recordings']))
    return pd.DataFrame(rows)


def part_d() -> dict:
    """Part D (the amendment of 2026-11-15 17:50): refuses to run before part C's verdict exists; reads part C's pre-odor
    level and the model's rows (`load_b1_all`, `with_a0`), D0's checks, the responses of part D's recordings and the 12
    from the maps and from Yang's tables, and writes d_responses.csv, d1_v2r.csv, d1_classes.csv, d1_lag_slope.csv,
    d2_replication.csv, d2_plane_pairs.csv, d3_pairs.csv, d3_repeats.csv, d3_reruns.csv, d_yang_table_ratios.csv,
    d_summary.json and provenance_D.json (with every input file's sha256)."""
    import hashlib
    if not part_c_done():
        raise RuntimeError("part D runs after part C's verdict (results/E98_apl_response_strength/verdict.json)")
    out = out_dir()
    verdict = _read_json(out / 'verdict.json')
    a0 = verdict.get('a0')
    d0 = part_d0()
    b1 = with_a0(load_b1_all(ALL_WIRINGS), a0)
    on = b1[b1['cond'] == 'on'][['wiring', 'stimulus', 'win', 'pre']]
    ranges = twelve_design()
    stems = list(D_RECORDINGS) + twelve()
    rows = d_rows(stems)
    means, tables, files = {}, {}, {}
    from experiments import E100_apl_weak_odors as E100
    for st, row in rows.iterrows():
        tables[st] = d_stim_table(d_yaml_path(row), st, D_EXTRA_OUT.get(st, ()))
        means[st] = d_frame_means(row, 'maps')
        if len(means[st]) != len(tables[st]):
            raise RuntimeError(f'{st}: {len(means[st])} frames for {len(tables[st])} presentations')
        pth = E100.recording_paths(row)
        for k in ('yaml', 'peak_maps', 'roi', 'roi_table'):
            if pth.get(k) is not None:
                files[str(pth[k])] = hashlib.sha256(Path(pth[k]).read_bytes()).hexdigest()
    resp = d_resp(stems, 'maps')
    resp.to_csv(out / 'd_responses.csv', index=False)
    resp_yang = d_resp(stems, 'yang')
    el = pd.read_csv(out / 'd0_eligibility.csv').set_index('stem')
    lags = el['v2_lag']
    v2r, classes, d1 = d1_v2r(resp, ranges, on)
    v2r.to_csv(out / 'd1_v2r.csv', index=False)
    classes.to_csv(out / 'd1_classes.csv', index=False)
    d1_lag_slope(resp, ranges, lags).to_csv(out / 'd1_lag_slope.csv', index=False)
    cats12 = d_twelve_categories(resp, ranges, on)
    d2 = d2_replication(resp, ranges, on, cats12)
    d2.to_csv(out / 'd2_replication.csv', index=False)
    d2_plane_pairs(resp, ranges).to_csv(out / 'd2_plane_pairs.csv', index=False)
    for k, t in d3_tables(resp, ranges, lags, means, tables).items():
        t.to_csv(out / f'{k}.csv', index=False)
    d_yang_table_ratios(resp_yang, ranges).to_csv(out / 'd_yang_table_ratios.csv', index=False)
    summ = dict(a0=a0, d0=dict(eligibility_matches=d0['eligibility_matches'], reproduce_12=d0['reproduce_12']), D1=d1,
                categories_12=cats12, D2_replicates=d2[d2['replicates'].notna()][['set', 'comparison', 'category',
                                                                                      'category_12', 'replicates']].to_dict('records'),
                time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(out / 'd_summary.json', summ)
    for f in (H.shared_dir() / 'data' / 'yang_apl_recordings' / 'recordings.csv', out / 'verdict.json', out / 'a1_responses.csv'):
        files[str(f)] = hashlib.sha256(Path(f).read_bytes()).hexdigest()
    b1_rec = load_b1_all(ALL_WIRINGS).attrs
    _write_json(out / 'provenance_D.json', dict(part='D', script=H.script_path(__file__), commits=process_commits(),
                                                time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), a0=a0, files=files,
                                                b1_sources=dict(sources=b1_rec.get('sources'), e101=b1_rec.get('e101'))))
    return summ


def comparison_table() -> str:
    """comparison_table.md from the parts' tables."""
    raise NotImplementedError


def _wirings(arg: str) -> Tuple[str, ...]:
    if arg == 'all8':
        return ALL_WIRINGS
    ws = tuple(x.strip() for x in arg.split(',') if x.strip())
    bad = [w for w in ws if w not in ALL_WIRINGS]
    if bad:
        raise SystemExit(f'unknown wirings {bad}')
    return ws


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'C', 'D0', 'D', 'table'])
    ap.add_argument('--wirings', default='all8', help="all8 or a comma-separated list of wiring labels")
    ap.add_argument('--arm', default='B1', help="part B's arms, comma-separated (B1..B4)")
    ap.add_argument('--threads', type=int, default=4,
                    help='part B (and part 0 with runs): OMP threads (OMP_NUM_THREADS and OMP_THREAD_LIMIT must match)')
    ap.add_argument('--no-runs', action='store_true', help="part 0: the reader's identities only")
    ap.add_argument('--allow-missing', action='store_true', help='part C: read the B1 runs there are')
    ap.add_argument('--a0', default=None, choices=list(A0_DEFINITIONS), help="part C: the model's pre-odor level (with_a0)")
    ap.add_argument('--root', default=None,
                    help="part B: the folder its runs go under (<root>/B/B1/...; default outputs/E98_apl_response_strength), "
                         "e.g. the archive's from_hal/ for another host's lane")
    a = ap.parse_args()
    process_commits(strict=a.part == 'B' or (a.part == '0' and not a.no_runs))
    ws = _wirings(a.wirings)
    if a.part in ('B',) or (a.part == '0' and not a.no_runs):
        for v in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT'):
            if os.environ.get(v) != str(a.threads):
                ap.error(f'part {a.part}: run with OMP_NUM_THREADS={a.threads} OMP_THREAD_LIMIT={a.threads} in the '
                         'environment')
    if a.part == '0':
        provenance('0', wirings=list(ws), runs=not a.no_runs)
        print(json.dumps(part_0(ws, runs=not a.no_runs), indent=1, default=float), flush=True)
    elif a.part == 'A':
        provenance('A')
        print(json.dumps(part_a(), indent=1, default=float), flush=True)
    elif a.part == 'B':
        arms = tuple(x for x in a.arm.split(',') if x)
        if a.root is None:
            provenance('B', arms=list(arms), wirings=list(ws), threads=a.threads)
        part_b(arms, ws, None if a.root is None else Path(a.root))
    elif a.part == 'C':
        provenance('C', allow_missing=a.allow_missing, a0=a.a0)
        print(json.dumps(part_c(ws, allow_missing=a.allow_missing, a0=a.a0), indent=1, default=float), flush=True)
    elif a.part == 'D0':
        provenance('D0')
        print(json.dumps(part_d0(), indent=1, default=float), flush=True)
    elif a.part == 'D':
        print(json.dumps(part_d(), indent=1, default=float), flush=True)
    elif a.part == 'table':
        print(comparison_table(), flush=True)


if __name__ == '__main__':
    main()
