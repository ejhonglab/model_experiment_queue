#!/usr/bin/env python3
"""E102: the antennal lobe's normalization at the strength the flies' APL asks for. The flies' APL falls with concentration
where the best family's holds its strength; in the family APL follows its claws, and the claws follow PN totals that the
model's antennal lobe makes nearly equal. What a weaker PN normalization (olfsysm's pn_inhadd raised, pn_inhsc with it)
costs or gains on the per-pair APL effect, at fixed scales on the eight wirings and with the family fitted again, and what
it does to the KCs' responses at low concentrations.

The entry is item E102 of experiment_priorities.md (proposed by atlas 2026-11-15 19:50; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-15 19:33; claimed by hal 19:53). C2's family (the breadth-2 APL-effect fit
of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit` /
`build_simulator`; S: real, P1, P2, M1; T: P3, P4, M2, M3). Halves A and B of the recordings are E77's, C and D E86's.
The input is the lab's in every verdict run.

Reused, not duplicated:
- E77 (experiments/E77_input_screen.py): `inhsc_for` (pn_inhsc raised with pn_inhadd as E3 / E32 did), `handle_breadth`
  (the breadth behind the APL->KC weights for a handle, here on each wiring's claw x bouton matrix), `hold_rates` (both
  response rates re-held), and the S1 points of its screen (part 0's second identity, part A3).
- E83 (experiments/E83_cross_wiring_fit.py): the wirings (`saved_fit`, `wiring_parts`, `build_simulator`, which takes the
  PN kwargs as `extra_kws` since this experiment), the scores of a population (`score_runs`, `effect_scores`,
  `pair_distances`, `realism`).
- E86: the targets with halves C and D (`targets`).
- E93 / E95: E73's per-KC reading under the saturating and the linear map (`E95_static_apl.zcell_distances_map`).
- E98 (experiments/E98_apl_response_strength.py): the extended input (`b1_input`, 52 stimuli with the blank), its saved
  layout (B1_LAYOUT; `save_b1`, `load_b1`), APL's and the PN stage's reading (`read_apl_pn`, `relabel`), the pre-odor
  levels (`with_a0`: the blank's window carries, the settled level is the partner, E90's level is shown and decides
  nothing; E98's amendment of 2026-11-15 17:59), the amended identity (`a0_identity`), the indicator family, `classify`,
  `comparison_verdict` and the verdict itself (`verdict_core`, read with the flies' intervals of E98's saved
  a3_ratios.csv and its saved pulse test and solvent rates: no recording is read again).
- E100: W1's six weak-odor ratios (`read_ratio`, `family_on`, `match_b1`, `estimate_labels`, `held_class`,
  `w1_verdict`), against the flies' intervals of its saved b_ratios.csv; the flies' KC response rates at low
  concentration (`kc_data_rates`, E72 part F's data_fractions.csv).
- E101: the claw stage on the extended input (E90's `read_stages`, `family_table`, `unit_meta`), X2's model side
  (`model_x2`, `x2_model_classes`) against the flies' intervals of its saved c_x2_flies.csv (pa over ms near 10^-3, and
  APL's ratio over the claws', Q).
- E90: the PN stage against the PN dendrites (`pn_statistics`: the totals' CV over the 12 megamat odors, the rms
  difference from the dendrites' group distances).

Definitions (the entry's, fixed before any reading):
- The normalization: olfsysm's PN gain is pn_inhsc / (pn_inhadd + 0.25 inhA + 0.75 inhB) inside the PN tanh. The ladder
  (LADDER): pn_inhadd 100, 314, 1000, 3141, 10000 and the gain fixed (pn_inhadd FIXED_GAIN_INHADD, 10^6), pn_inhsc
  raised with it (`E77_input_screen.inhsc_for`) so that the gain at a typical inhibition is kept. The identity point
  (STEP_IDENTITY) passes the defaults (31.4088, 368.6631) explicitly.
- Rates re-held (E77's screen; `hold_point`): the breadth weights recomputed for the step on each wiring's own wiring
  (`point_breadth`), the somatic threshold bisected to the wiring's fit's count of responses with APL off, then the
  tied APL scale (APL->KC with KC->APL and bouton->APL) to its count with APL on; the KC-KC and APL->bouton scales as
  fitted.
- Each point is also run at its held parameters on E98's extended input (52 stimuli and the blank), APL on and acutely
  off (the APL-off run at hold_rates' params_off, as E98's B4), with PN rates, APL and the claws' summed drive saved
  (`run_extended`).
- The model's APL is read as E98 reads it: against the blank's window in the same run (the carrying level), with the
  settled pre-odor level as partner; E90's pre-odor mean decides nothing (unresolved_issues.md, the pre-odor start-up).
- Every pattern part is given with its amplitude part, its shape part, the spread ratio, ρ and Pearson r on the per-pair
  APL effect, the loss and the mean APL effect, under the calibrated readout (K set on each run's APL-on counts) and on
  raw counts. No level is chosen on the main target or on a correlation.

The levels (part C; `levels`, fixed here):
- L_apl: the smallest pn_inhadd at which V2 (paa at 10^-5 over 10^-4) is inside by E98's rule (`verdict_core`: at least 6
  of the 8 wirings, the carrying reading and level, held under its partners).
- L_lin: the smallest at which V2's linear reading (n 1, evoked; carrying reading and level) lies in the flies' interval
  on at least LEVEL_MIN_WIRINGS (6) of the 8 wirings.
- L_free: the largest that is an equal at fixed scales on S by E94's rule (`equal_rule`).
- A step with a flagged point (any of its eight runs outside the realism bands) is not a level. Any of the three may not
  exist. If L_apl does not exist, the refit runs at the last unflagged step (`refit_levels`).

The classes (part D; `refit_class`, fixed here), against the unchanged fits on S (BASE_PATTERN, recomputed in part 0):
- an equal (E94's rule): the pattern part within EQUAL_PATTERN (0.003) against all recordings and against halves A and B,
  the four fits inside the bands, the raw-count loss within EQUAL_RAW (0.005), the mean APL effect within EQUAL_MEAN
  (0.01);
- better: lower by BETTER_MARGIN (0.003) or more against all recordings and lower against each of halves A to D;
- a price: anything else.
- On an equal or better, T's four wirings are refitted at that level; better stands only if T's population is also lower
  than its fits by BETTER_MARGIN or more against all recordings and lower against each half, at least
  T_BETTER_MIN_WIRINGS (6) of the 8 single wirings are lower and the raw-count loss is not higher (`t_better_holds`).

The KCs at low concentration (outcome (e); `low_conc_class`, fixed here): the model's response rate at ma 10^-7 (the
kiwi-control session's estimate, E49's source, carrying; the other sessions beside), the median over the eight wirings,
against the flies' (E72 part F; control flies with APL on, APL output blocked (TNT) with APL off). The distance from
the flies is |log(model / flies)|; at a step it is 'further' if it exceeds C2's (the identity point's) by a factor of
LOW_CONC_FACTOR (1.25) or more, 'nearer' if C2's exceeds it by that factor, 'about the same' otherwise.

Parts (CLI `--part`):
- 0, the identities (identity.json; a failure stops the experiment, and every later part refuses to run past it,
  `require_identity`): the ladder's identity point reproduces each wiring's saved fit with 0 count differences (16-odor
  input, APL on and off; factors 1) and, on the real wiring, E98's saved B1 run on the extended input (counts and APL's
  means; `identity_point`); pn_inhadd 100 to 3141 on the real wiring reproduce E77's screen (thresholds, scales, counts and
  scores; `e77_identity`, read from part B's runs once they exist); the extended-input runs at 1000 and 3141 on the real
  wiring reproduce E98's B4 (counts and APL window means on their common stimuli; `b4_identity`), and every run passes
  E98's amended blank identity (`a0_identity`).
- A (no simulations): A1 the chain of totals from the ORN input to APL per stimulus on the real wiring at C2 (E98's B1
  run), 1000 and 3141 (E98's B4 runs): the input's total, its broad component (36 x the median glomerulus) and the rest,
  the PN total and APL's response at the carrying level, as ratios for E98's, E100's and E101's stimuli, with the exponent
  a power law from the input's total to APL would need, in the flies and in the model, and the log split of the model's
  flattening into the antennal lobe's part and the mushroom body's (a1_chain.csv, a1_ratios.csv); A2 the KCs at low
  concentration under E98's two B4 settings: response rates with APL on and off at ma 10^-7 and 2h 10^-5 (exact input),
  2h 10^-6 and farn 10^-2 (standing for 10^-7 and 10^-3), against the flies' and C2's (a2_low_conc.csv); A3 E77's S1
  points under the calibrated readout and on raw counts with the amplitude and shape parts, from its saved pair distances,
  the per-odor part beside as a lead (a3_s1_parts.csv, a3_s1_per_odor.csv); A4 the prediction for every ladder step on S
  (pattern part, mean APL effect on raw counts, the paa ratio), written before part B (a4_prediction.csv,
  a4_prediction.json). a_summary.json.
- B (simulations; resumable; `--wirings`, `--steps`): each step on each wiring, rates re-held (`hold_point`), its 16-odor
  outputs and parameters saved, then the extended run (`run_extended`), under outputs/E102_antennal_lobe_normalization/B/
  <step>/ in E98's B1 layout plus <wiring>_family.csv (E101's family table) and <wiring>_point.json.
- score (no simulations; from part B's saved runs): every step's readings (b_readings.csv per unit and readout and
  target; b_points.csv per wiring; b_apl.csv E98's, E100's and E101's comparisons per step; b_low_conc.csv; b_pn.csv) and
  part C's levels (levels.json).
- D (the refits): the commands of `E10_joint_fit.py` with C2's options at each level (`refit_commands`; printed, not run),
  and their classes once fitted (`--score-refits`).

Outputs: results/E102_antennal_lobe_normalization/ (the tables above, identity.json, provenance_<part>.json); counts,
responses, APL and PN means of part B under outputs/E102_antennal_lobe_normalization/, archived to
data/model_outputs/E102_antennal_lobe_normalization/ with pn_inhadd, pn_inhsc, the thresholds and every scale as absolute
values. Every table holding Yang's recordings' values (the flies' intervals and rates) opens with the quality caveat
(YANG_QUALITY_CAVEAT; her answers of 2026-11-15 18:08) as a '#' line.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E102_antennal_lobe_normalization --part A
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E102_antennal_lobe_normalization --part 0
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E102_antennal_lobe_normalization --part B --wirings S
    ./venv/bin/python -u -m experiments.E102_antennal_lobe_normalization --part score
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E102_antennal_lobe_normalization'
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
ALL_WIRINGS = S_WIRINGS + T_WIRINGS
WIRING_SETS = {'S': S_WIRINGS, 'T': T_WIRINGS, 'all8': ALL_WIRINGS}

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any reading)

PN_INHADD_DEFAULT, PN_INHSC_DEFAULT = 31.4088, 368.6631     # olfsysm's (E3's PN_INHADD_DEFAULT, PN_INHSC_DEFAULT)
FIXED_GAIN_INHADD = 1e6                                     # the ladder's last step: the gain fixed
LADDER = (100.0, 314.0, 1000.0, 3141.0, 10000.0, FIXED_GAIN_INHADD)
STEP_IDENTITY = 'identity'
E77_S1 = (100.0, 314.0, 1000.0, 3141.0)                     # E77's screen's S1 points (part 0's second identity)
E98_B4 = (1000.0, 3141.0)                                   # E98's B4 points (part 0's third identity)

# the levels (part C) and the classes (part D)
LEVEL_MIN_WIRINGS = 6                                       # of the 8 (E98's VERDICT_MIN_WIRINGS)
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01   # E94's equal rule
EQUAL_TARGETS = ('z8', 'halfA', 'halfB')                    # an equal: the pattern part against all recordings, A and B
BETTER_MARGIN = 0.003
BETTER_HALVES = ('halfA', 'halfB', 'halfC', 'halfD')
T_BETTER_MIN_WIRINGS = 6                                    # of the 8 single wirings (E83's gate)
# the unchanged fits' pattern parts under the calibrated readout against all recordings (the entry's; part 0 recomputes)
BASE_PATTERN = {'S': 0.0500, 'T': 0.0486, 'all8': 0.0482}
BASE_PATTERN_TOL = 5e-5                                     # part 0: the recomputed values agree to the entry's 4 decimals

# the KCs at low concentration (outcome (e))
LOW_CONC_FACTOR = 1.25
LOW_STIMULI = (                                             # (B1 label of the model's input, the flies' stimulus, role)
    ('ma [kiwi-control/diag] @ -7', ('ma', -7.0), 'exact'),
    ('2h [kiwi-control/control] @ -5', ('2h', -5.0), 'exact'),
    ('2h [kiwi-control/diag] @ -6', ('2h', -7.0), 'stands for 10^-7'),
    ('farn [kiwi-control/diag] @ -2', ('farn', -3.0), 'stands for 10^-3'))
LOW_CARRYING = 'ma [kiwi-control/diag] @ -7'                # outcome (e)'s stimulus (E49's source)
LOW_BESIDE = ('ma [megamat/diag] @ -7', 'ma [validation2/diag] @ -7')
FLIES_OFF_ON_MA7 = 2.29                                     # the flies' APL-off over APL-on response rate at ma 10^-7

# refits (part D): E10_joint_fit.py with C2's options
REFIT_ARGS = ('-k', '2', '--aplpn-max', '3', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both')
REFIT_MAXFEV = 200

YANG_QUALITY_CAVEAT = ("Yang has not checked these recordings' quality (her answers of 2026-11-15 18:08: all data can be "
                       "used, she is unsure of quality questions); the flies' values are read from E98's, E100's, E101's "
                       "and E72's saved tables")


# ---------------------------------------------------------------------------
# folders, records

def out_dir() -> Path:
    """results/E102_antennal_lobe_normalization/ in the shared directory."""
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    """outputs/E102_antennal_lobe_normalization/<sub>/ in this repository."""
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _is_identity(step) -> bool:
    return isinstance(step, str) and step == STEP_IDENTITY


def step_tag(step) -> str:
    """A step's folder and tag name: 'identity', 'inhadd100' ... 'inhadd10000', 'gain_fixed' (ASCII only)."""
    if _is_identity(step):
        return STEP_IDENTITY
    a = float(step)
    if a == FIXED_GAIN_INHADD:
        return 'gain_fixed'
    return f'inhadd{a:g}'


def step_value(step) -> float:
    """A step's pn_inhadd (the default at STEP_IDENTITY)."""
    return PN_INHADD_DEFAULT if _is_identity(step) else float(step)


def step_of_tag(tag: str):
    """The step a tag names (`step_tag`'s inverse)."""
    if tag == STEP_IDENTITY:
        return STEP_IDENTITY
    if tag == 'gain_fixed':
        return FIXED_GAIN_INHADD
    if not tag.startswith('inhadd'):
        raise ValueError(f'not a step tag: {tag!r}')
    return float(tag[len('inhadd'):])


def step_kws(step) -> Dict[str, float]:
    """A step's PN kwargs: pn_inhadd and pn_inhsc (`E77_input_screen.inhsc_for`; the defaults at STEP_IDENTITY)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT as A0, PN_INHSC_DEFAULT as S0
    if (A0, S0) != (PN_INHADD_DEFAULT, PN_INHSC_DEFAULT):
        raise RuntimeError(f"E3's defaults ({A0}, {S0}) are not this module's ({PN_INHADD_DEFAULT}, {PN_INHSC_DEFAULT})")
    if _is_identity(step):
        return dict(pn_inhadd=A0, pn_inhsc=S0)
    from experiments.E77_input_screen import inhsc_for
    a = float(step)
    return dict(pn_inhadd=a, pn_inhsc=float(inhsc_for(a)))


def write_table(df: pd.DataFrame, path: Path, caveat: bool = False) -> None:
    """A results table; with `caveat` its first line is '# ' + YANG_QUALITY_CAVEAT (report_values skips it)."""
    path = Path(path)
    with open(path, 'w') as f:
        if caveat:
            f.write(f'# {YANG_QUALITY_CAVEAT}\n')
        df.to_csv(f, index=False)


def read_table(path: Path) -> pd.DataFrame:
    """A results table written by `write_table` (leading '#' lines skipped)."""
    n = 0
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            n += 1
    return pd.read_csv(path, skiprows=n, float_precision='round_trip')


def _write_json(path: Path, obj) -> None:
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')


def _json_default(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (pd.Series, pd.Index, np.ndarray)):
        return list(np.asarray(x).tolist())
    return str(x)


def _read_json(path: Path) -> Optional[dict]:
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else None


def _host() -> str:
    return socket.gethostname().split('.')[0]


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                                         **extra))
    _write_json(path, d)


def dirty_code(repo: Path = H.HERE) -> List[str]:
    """The checkout's modified tracked files that are code (.py, .cpp, .hpp, .h, .toml, .cfg): what a run's outputs could
    depend on. Other modified files (a report template another session is editing) do not count."""
    import subprocess
    st = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=repo, capture_output=True, text=True,
                        check=True).stdout
    paths = [ln[3:].strip() for ln in st.splitlines() if ln.strip()]
    return [x for x in paths if x.endswith(('.py', '.cpp', '.hpp', '.h', '.toml', '.cfg'))]


_RUN_COMMITS: Dict[str, object] = {}


def run_commits() -> Dict[str, object]:
    """The commits of the code a run loads (read once per process), refusing modified code in this repository or in
    pyolfsysm (`dirty_code`); the modified non-code files of this repository are listed ('queue_dirty_files')."""
    if not _RUN_COMMITS:
        import subprocess
        c = {k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}
        bad = dirty_code()
        olf = Path(H.pyolfsysm_dir()) if hasattr(H, 'pyolfsysm_dir') else None
        bad_olf = dirty_code(olf) if olf is not None and olf.exists() else []
        if bad or bad_olf:
            raise RuntimeError(f'part B runs from committed code only: modified {bad + bad_olf}')
        st = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=H.HERE, capture_output=True,
                            text=True).stdout
        c['queue_dirty_files'] = ';'.join(ln[3:].strip() for ln in st.splitlines() if ln.strip())
        _RUN_COMMITS.update(c)
    return dict(_RUN_COMMITS)


def identity_path() -> Path:
    return out_dir() / 'identity.json'


class _Lock:
    """A lock beside a shared file (`mkdir` is atomic on NFS), for identity.json's read-modify-write when two lanes run
    part 0 at once; a lock older than LOCK_STALE_S is taken as stale."""

    def __init__(self, path: Path, timeout: float = 600.0):
        self.d, self.timeout = Path(str(path) + '.lock'), timeout

    def __enter__(self):
        t0 = time.time()
        while True:
            try:
                self.d.mkdir()
                return self
            except FileExistsError:
                try:
                    if time.time() - self.d.stat().st_mtime > LOCK_STALE_S:
                        self.d.rmdir()
                        continue
                except FileNotFoundError:
                    continue
                if time.time() - t0 > self.timeout:
                    raise TimeoutError(f'{self.d} held for more than {self.timeout} s')
                time.sleep(1.0)

    def __exit__(self, *exc):
        try:
            self.d.rmdir()
        except FileNotFoundError:
            pass


LOCK_STALE_S = 600.0


def update_identity(**entries) -> dict:
    """identity.json updated under its lock: each keyword replaces that entry; 'points' entries are merged by wiring."""
    path = identity_path()
    with _Lock(path):
        ident = _read_json(path) or {}
        for k, v in entries.items():
            if k == 'points':
                ident.setdefault('points', {}).update(v)
            else:
                ident[k] = v
        ident['time'], ident['host'] = time.strftime('%Y-%m-%d %H:%M'), _host()
        _write_json(path, ident)
    return ident


def require_identity(which: str = 'all', wirings: Sequence[str] = ALL_WIRINGS, path: Optional[Path] = None) -> None:
    """Refuses to go on past a failed (or missing) identity of part 0 (identity.json): 'point' needs the identity point of
    each of `wirings` passed; 'readers' the no-simulation identities (the unchanged fits' populations, the readings of
    E98's and E100's saved runs); 'all' both, and the E77 / B4 identities wherever they have been read."""
    ident = _read_json(identity_path() if path is None else path)
    if ident is None:
        raise RuntimeError('part 0 has not run (no identity.json): refusing to go on')
    if isinstance(wirings, str):
        wirings = [wirings]
    bad = []
    if which in ('point', 'all'):
        for w in wirings:
            pt = (ident.get('points') or {}).get(w)
            if pt is None or not pt.get('exact'):
                bad.append(f'identity point {w}: {"not run" if pt is None else "failed"}')
    if which in ('readers', 'all'):
        for k in ('base_population', 'e98_reading', 'e100_reading'):
            r = ident.get(k)
            if r is None or not r.get('passes'):
                bad.append(f'{k}: {"not run" if r is None else "failed"}')
    if which == 'all':
        for k in ('e77', 'b4'):
            for step, r in (ident.get(k) or {}).items():
                if not r.get('exact'):
                    bad.append(f'{k} {step}: failed')
    if bad:
        raise RuntimeError('refusing to go on past a failed identity: ' + '; '.join(bad))


# ---------------------------------------------------------------------------
# part B: one point (a step on a wiring)

def step_dir(step, root: Optional[Path] = None) -> Path:
    """The folder of one step's runs: `root`/B/<step_tag> (default outputs/E102_antennal_lobe_normalization/B/...)."""
    d = (Path(root) if root is not None else local_dir()) / 'B' / step_tag(step)
    d.mkdir(parents=True, exist_ok=True)
    return d


def point_breadth(wiring: str, kws: Mapping[str, float]) -> pd.Series:
    """The breadth behind the APL->KC weights at a step on one wiring (`E77_input_screen.handle_breadth` with the
    wiring's claw x bouton matrix)."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    return S.handle_breadth(dict(kws), wPNKC=E83.wiring_parts(wiring)['W'])


def hold_point(wiring: str, step) -> dict:
    """Both response rates re-held at a step on one wiring (`E77_input_screen.hold_rates` on E83's simulator with the
    step's PN kwargs and breadth; the wiring's fit's targets and parameters). Returns hold_rates' dict with 'pn_kws',
    'step', 'wiring' and 'breadth_series' (the breadth used)."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    P = E83.saved_fit(wiring)
    if P['size'] != S.c2_parts()['size']:
        raise RuntimeError(f"{wiring}: {P['size']} KC x odor entries, C2's {S.c2_parts()['size']} (hold_rates' rates)")
    kws = step_kws(step)
    breadth = point_breadth(wiring, kws)
    res = S.hold_rates(lambda: E83.build_simulator(wiring, breadth=breadth, params=P['p'], extra_kws=kws),
                       target_off=P['target_off'], target_on=P['target_on'], params=P['p'],
                       label=f'{wiring} {step_tag(step)}')
    res.update(pn_kws=kws, step=step_tag(step), wiring=wiring, breadth_series=breadth)
    return res


def run_extended(wiring: str, step, held: Mapping, breadth: pd.Series) -> Tuple[Dict[str, dict], pd.DataFrame, dict]:
    """The held point on E98's extended input (`E98.b1_input`), APL on (held params) and acutely off (hold_rates'
    params_off, as E98's B4): per condition tag ('on', 'acute_off') E98's outputs (responses, spike counts, APL's means
    with the settled level, PN window and pre-odor means, APL's trace) and E90's claw stage ('claw_win', 'claw_pre',
    `read_stages`, whose APL means must equal E98's). Returns (outputs, B1's stimulus table, the claws' metadata)."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E90_model_stages as E90
    from experiments import E98_apl_response_strength as E98
    from experiments import E101_apl_vs_kc_claws as E101
    orn, table = E98.b1_input()
    labels = E98.input_labels(orn)
    sim, info = E83.build_simulator(wiring, orn=orn, breadth=breadth, params=dict(held['params']), extra_kws=step_kws(step))
    mp, rv, order = E98._sim_cells(sim)
    W = info['W']
    gl = pd.Index(list(dict.fromkeys(W.columns.get_level_values('glomerulus'))), name='glomerulus')
    out = {}
    for c, ct in E98.CONDS:
        q = dict(held['params']) if c == 'on' else dict(held['params_off'])
        t1 = time.time()
        r, s = sim(q)
        o = E98.relabel(r, s, labels)
        o.update(E98.read_apl_pn(mp, rv, order, labels, gl))
        st = E90.read_stages(mp, rv, W, order, labels, E90.MODEL_WINDOW)
        d = (st['apl'][['win', 'pre']].astype(float) - o['apl'].loc[st['apl'].index, ['win', 'pre']].astype(float)).abs()
        if float(d.to_numpy().max()) > E98.APL_TOL:
            raise RuntimeError(f"{wiring} {step_tag(step)} {ct}: E90's and E98's readings of APL differ")
        o['claw_win'], o['claw_pre'] = st['claw_win'], st['claw_pre']
        o['params'], o['seconds'] = q, time.time() - t1
        out[ct] = o
        print(f'{wiring} {step_tag(step)} extended {ct}: {o["seconds"]:.0f} s', flush=True)
    meta = E101.unit_meta(W)
    del sim, info, mp, rv
    return out, table, meta


def _held_paths(d: Path, wiring: str) -> Dict[str, Tuple[Path, Path]]:
    return {c: (d / f'{wiring}_held_{ct}_responses.parquet', d / f'{wiring}_held_{ct}_spike_counts.parquet')
            for c, ct in (('on', 'on'), ('off', 'acute_off'))}


def load_held(step, wiring: str, root: Optional[Path] = None) -> Optional[Dict[str, tuple]]:
    """A point's 16-odor outputs ({'on', 'off'}: (responses, spike counts)), or None."""
    d = (Path(root) if root is not None else local_dir()) / 'B' / step_tag(step)
    ps = _held_paths(d, wiring)
    if not all(a.exists() and b.exists() for a, b in ps.values()):
        return None
    return {c: (pd.read_parquet(a), pd.read_parquet(b)) for c, (a, b) in ps.items()}


def load_point(step, wiring: str, root: Optional[Path] = None) -> Optional[dict]:
    """A point's record (<wiring>_point.json), or None."""
    return _read_json((Path(root) if root is not None else local_dir()) / 'B' / step_tag(step) / f'{wiring}_point.json')


def _differ(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    """Counts and responses of two (responses, spike counts) pairs: entries that differ, and whether the index and
    columns are the same (columns compared as abbreviated odor names)."""
    ra, sa = (H.abbrev_columns(x.astype(float)) for x in a)
    rb, sb = (H.abbrev_columns(x.astype(float)) for x in b)
    same = bool(sa.index.equals(sb.index) and list(map(str, sa.columns)) == list(map(str, sb.columns)))
    if not same:
        return dict(index_equal=False, counts_differ=-1, responses_differ=-1)
    return dict(index_equal=True, counts_differ=int((sa.to_numpy() != sb.to_numpy()).sum()),
                responses_differ=int((ra.to_numpy() != rb.to_numpy()).sum()))


def identity_point_checks(wiring: str, held: Mapping, breadth: pd.Series) -> dict:
    """The ladder's identity point on one wiring: factors 1, 0 count and response differences against the saved fit
    (16 odors, APL on and off), the breadth against the wiring's (largest absolute difference)."""
    from experiments import E83_cross_wiring_fit as E83
    P = E83.wiring_parts(wiring)
    res = dict(thr_factor=float(held['thr_factor']), apl_factor=float(held['apl_factor']))
    for c in ('on', 'off'):
        res[c] = _differ(held[c], P['saved'][c])
    b0 = P['breadth']
    res['breadth_index_equal'] = bool(breadth.index.equals(b0.index))
    res['breadth_max_abs_diff'] = float((breadth.to_numpy(float) - b0.to_numpy(float)).__abs__().max()) \
        if res['breadth_index_equal'] else float('nan')
    res['exact'] = bool(res['thr_factor'] == 1.0 and res['apl_factor'] == 1.0
                        and all(res[c]['index_equal'] and res[c]['counts_differ'] == 0 and res[c]['responses_differ'] == 0
                                for c in ('on', 'off')))
    return res


B1_COMPARE = ('win', 'pre', 'pre_settled', 'response_rate')


def b1_rows(step, wirings: Sequence[str] = ALL_WIRINGS, root: Optional[Path] = None) -> pd.DataFrame:
    """A step's extended runs as E98's `load_b1` gives them (one row per wiring, condition and stimulus; APL's window,
    pre-odor and settled means, the KCs' response rate, the parsed label); attrs['missing'] as E98's."""
    from experiments import E98_apl_response_strength as E98
    d = (Path(root) if root is not None else local_dir()) / 'B' / step_tag(step)
    return E98.load_b1(wirings, roots=[d])


def extended_identity(wiring: str, d: Path) -> dict:
    """The identity point's extended run against E98's saved B1 run (the real wiring) or E101's A3 rows (the other
    wirings): APL's window, pre-odor and settled means and the KCs' response rate per condition and stimulus (exact)."""
    from experiments import E98_apl_response_strength as E98
    mine = E98.load_b1([wiring], roots=[d])
    if E98.load_b1_run(wiring) is not None:
        ref, src = E98.load_b1([wiring]), "E98's B1"
    else:
        eb, _ = E98.load_b1_e101()
        ref, src = eb[eb['wiring'] == wiring], "E101's A3"
    a = mine.set_index(['cond', 'stimulus'])[list(B1_COMPARE)].sort_index()
    b = ref.set_index(['cond', 'stimulus'])[list(B1_COMPARE)].sort_index()
    same = bool(a.index.equals(b.index))
    dapl = float((a[list(B1_COMPARE[:3])] - b[list(B1_COMPARE[:3])]).abs().to_numpy().max()) if same else float('nan')
    drate = float((a['response_rate'] - b['response_rate']).abs().max()) if same else float('nan')
    return dict(source=src, index_equal=same, apl_max_abs_diff=dapl, response_rate_max_abs_diff=drate, n=int(len(a)),
                exact=bool(same and dapl <= E98.APL_TOL and drate == 0.0))


def e77_point_identity(step, held: Mapping) -> dict:
    """pn_inhadd 100 to 3141 on the real wiring against E77's screen (point S1_inhadd<a>): the absolute parameters with
    APL on and off (exact), the response counts, the spike counts against the archived screen outputs (0 differences),
    and the scores E77 recorded (E72's `effect_stats` under the calibrated readout and raw counts, within SCORE_TOL)."""
    from experiments import E72_wiring_population as E72
    from experiments import E77_input_screen as S
    pid = f'S1_inhadd{float(step):g}'
    res_dir = H.shared_dir() / 'results' / S.EXPERIMENT
    sp = json.loads((res_dir / 'screen_params.json').read_text())[pid]
    row = pd.read_csv(res_dir / 'screen.csv').set_index('point_id').loc[pid]
    out = dict(point=pid)
    for k, mine, theirs in (('params', held['params'], sp['params']), ('params_off', held['params_off'], sp['params_off'])):
        theirs = {a: float(v) for a, v in (json.loads(theirs) if isinstance(theirs, str) else theirs).items()}
        out[f'{k}_max_abs_diff'] = float(max(abs(float(mine[a]) - theirs[a]) for a in theirs))
        out[f'{k}_equal'] = bool(all(float(mine[a]) == theirs[a] for a in theirs) and set(mine) == set(theirs))
    out['n_on'], out['n_off'] = int(held['n_on']), int(held['n_off'])
    out['n_equal'] = bool(out['n_on'] == int(row['n_on']) and out['n_off'] == int(row['n_off']))
    arc = H.shared_dir() / 'data' / 'model_outputs' / S.EXPERIMENT / 'screen'
    local = H.outputs_dir(S.EXPERIMENT) / 'screen'
    saved = None
    for dd in (local, arc):
        fs = {c: (dd / f'C2_{pid}_{ct}_responses.parquet', dd / f'C2_{pid}_{ct}_spike_counts.parquet')
              for c, ct in (('on', 'on'), ('off', 'acute_off'))}
        if all(a.exists() and b.exists() for a, b in fs.values()):
            saved = {c: (pd.read_parquet(a), pd.read_parquet(b)) for c, (a, b) in fs.items()}
            out['saved_dir'] = str(dd)
            break
    for c in ('on', 'off'):
        out[c] = _differ(held[c], saved[c]) if saved is not None else dict(note='no saved screen outputs')
    o = {c: held[c] for c in ('on', 'off')}
    sc = {}
    for ro in ('raw', 'cal'):
        st = E72.effect_stats(E72.distances(o, ro), ro)
        for k in ('loss', 'mean', 'pattern', 'rho'):
            col = f'{ro}_{k}'
            if col in row.index:
                sc[col] = dict(mine=float(st[k]), screen=float(row[col]), diff=float(st[k] - float(row[col])))
    out['scores'] = sc
    out['scores_equal'] = bool(all(abs(v['diff']) <= SCORE_TOL for v in sc.values()))
    out['exact'] = bool(out['params_equal'] and out['params_off_equal'] and out['n_equal'] and saved is not None
                        and all(out[c].get('counts_differ') == 0 and out[c].get('responses_differ') == 0 for c in ('on', 'off'))
                        and out['scores_equal'])
    return out


SCORE_TOL = 1e-9


def b4_point_identity(step, d: Path, wiring: str = 'real') -> dict:
    """The extended run at 1000 or 3141 on the real wiring against E98's archived B4 run (E77_inhadd<a>): spike counts and
    responses on their common stimuli (0 differences), APL's window, pre-odor and settled means (within E98's APL_TOL)."""
    from experiments import E98_apl_response_strength as E98
    name = f'E77_inhadd{float(step):g}'
    roots = [E98.b_root() / 'B4', H.shared_dir() / 'data' / 'model_outputs' / E98.EXPERIMENT / 'B' / 'B4']
    ref = E98.load_b1_run(name, roots=roots)
    mine = E98.load_b1_run(wiring, roots=[d])
    out = dict(point=name)
    if ref is None or mine is None:
        out.update(exact=False, note='no saved B4 run' if ref is None else 'no extended run')
        return out
    out['ref_dir'] = ref['on']['dir']
    ok = True
    for _, ct in E98.CONDS:
        a, b = mine[ct], ref[ct]
        common = [x for x in b['spike_counts'].columns if x in set(a['spike_counts'].columns)]
        sa, sb = a['spike_counts'][common].astype(float), b['spike_counts'][common].astype(float)
        ra, rb = a['responses'][common].astype(float), b['responses'][common].astype(float)
        same = bool(np.array_equal(np.asarray(sa.index), np.asarray(sb.index)))
        nd = int((sa.to_numpy() != sb.to_numpy()).sum()) if same else -1
        nr = int((ra.to_numpy() != rb.to_numpy()).sum()) if same else -1
        cols = [k for k in ('win', 'pre', 'pre_settled') if k in b['apl'].columns]
        dapl = float((a['apl'].loc[common, cols].astype(float) - b['apl'].loc[common, cols].astype(float)).abs().to_numpy().max())
        out[ct] = dict(n_common=len(common), index_equal=same, counts_differ=nd, responses_differ=nr, apl_max_abs_diff=dapl)
        ok = ok and same and nd == 0 and nr == 0 and dapl <= E98.APL_TOL
    out['exact'] = bool(ok)
    return out


def run_point(wiring: str, step, root: Optional[Path] = None, force: bool = False, extended: bool = True) -> dict:
    """One point of part B (resumable): `hold_point`, the 16-odor outputs and the parameters saved, then `run_extended`
    saved in E98's B1 layout with E101's family table; checks that the extended run's 16 panel columns equal the held
    run's counts, and E98's amended blank identity on both conditions. At the identity step, the identity point's checks
    (`identity_point_checks`, `extended_identity`); on the real wiring at E77's S1 steps, `e77_point_identity`; at E98's
    B4 steps, `b4_point_identity`. A failed check is recorded and stops the run (RuntimeError). Any other step refuses to
    run until the wiring's identity point has passed (`require_identity`). Returns the point's record
    (<wiring>_point.json)."""
    import gc
    from experiments import E98_apl_response_strength as E98
    from experiments import E101_apl_vs_kc_claws as E101
    d = step_dir(step, root)
    rec_path = d / f'{wiring}_point.json'
    if not force and rec_path.exists():
        rec = _read_json(rec_path)
        if rec.get('complete') and (rec.get('extended_done') or not extended):
            return rec
    if not _is_identity(step):
        require_identity('point', [wiring])
    commits = run_commits()
    t0 = time.time()
    held = hold_point(wiring, step)
    breadth = held.pop('breadth_series')
    for c, (fr, fs) in _held_paths(d, wiring).items():
        held[c][0].to_parquet(fr)
        held[c][1].to_parquet(fs)
    kws = held['pn_kws']
    rec = dict(wiring=wiring, step=step_tag(step), pn_inhadd=kws['pn_inhadd'], pn_inhsc=kws['pn_inhsc'],
               params=held['params'], params_off=held['params_off'], thr_factor=float(held['thr_factor']),
               apl_factor=float(held['apl_factor']), n_on=int(held['n_on']), n_off=int(held['n_off']),
               target_on=int(held['target_on']), target_off=int(held['target_off']),
               at_bound_thr=bool(held['at_bound_thr']), at_bound_apl=bool(held['at_bound_apl']), n_sims=int(held['n_sims']),
               breadth_sha256=held['breadth_sha256'], sim_kws=held.get('kws'), seconds_hold=time.time() - t0,
               history=held.get('history', []), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
               threads=os.environ.get('OMP_NUM_THREADS'), thread_limit=os.environ.get('OMP_THREAD_LIMIT'), commits=commits,
               complete=False, extended_done=False)
    checks_ok = True
    if _is_identity(step):
        rec['identity_point'] = identity_point_checks(wiring, held, breadth)
        checks_ok &= rec['identity_point']['exact']
    if wiring == 'real' and not _is_identity(step) and float(step) in E77_S1:
        rec['e77_identity'] = e77_point_identity(step, held)
        checks_ok &= rec['e77_identity']['exact']
    _write_json(rec_path, rec)
    print(f"{wiring} {step_tag(step)}: held in {rec['seconds_hold']:.0f} s (thr x{rec['thr_factor']:.4f}, APL x"
          f"{rec['apl_factor']:.4f}; {rec['n_off']} / {rec['n_on']} responses off / on)", flush=True)
    if not checks_ok:
        raise RuntimeError(f'{wiring} {step_tag(step)}: an identity failed: '
                           + json.dumps({k: rec.get(k) for k in ('identity_point', 'e77_identity')}, default=str))
    if extended:
        t1 = time.time()
        out, table, meta = run_extended(wiring, step, held, breadth)
        tmp = d / f'.b1_stimuli_{wiring}.csv'
        table.to_csv(tmp, index=False)
        os.replace(tmp, d / 'b1_stimuli.csv')
        _write_json(d / f'{wiring}_units.json', meta)
        rows = []
        for c, ct in E98.CONDS:
            o = out[ct]
            E98.save_b1(d, wiring, ct, o)
            same = _differ(held[c], (o['responses'][list(H.APL16_ODORS)], o['spike_counts'][list(H.APL16_ODORS)]))
            rec[f'extended_vs_held_{c}'] = same
            rec[f'a0_identity_{c}'] = E98.a0_identity(o['apl'])
            rows.append(dict(wiring=wiring, step=step_tag(step), cond=ct, n_stimuli=int(o['spike_counts'].shape[1]),
                             runtime_s=o['seconds'], host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                             threads=os.environ.get('OMP_NUM_THREADS'), pn_inhadd=kws['pn_inhadd'], pn_inhsc=kws['pn_inhsc'],
                             **{f'param_{k}': v for k, v in o['params'].items()},
                             **{k: v for k, v in commits.items() if not isinstance(v, (dict, list))}))
        fam = E101.family_table(wiring, out, meta, list(H.APL16_ODORS))
        fam.to_csv(d / f'{wiring}_family.csv', index=False)
        pd.DataFrame(rows).to_csv(d / f'{wiring}_runs.csv', index=False)
        rec['seconds_extended'] = time.time() - t1
        ext_ok = all(rec[f'extended_vs_held_{c}'].get('counts_differ') == 0
                     and rec[f'extended_vs_held_{c}'].get('responses_differ') == 0 for c in ('on', 'off'))
        a0_ok = all(rec[f'a0_identity_{c}']['passes'] for c in ('on', 'off'))
        rec['extended_matches_held'], rec['a0_identity_passes'] = bool(ext_ok), bool(a0_ok)
        checks_ok &= ext_ok and a0_ok
        if _is_identity(step):
            rec['extended_identity'] = extended_identity(wiring, d)
            checks_ok &= rec['extended_identity']['exact']
        if wiring == 'real' and not _is_identity(step) and float(step) in E98_B4:
            rec['b4_identity'] = b4_point_identity(step, d, wiring)
            checks_ok &= rec['b4_identity']['exact']
        rec['extended_done'] = True
        del out
    rec['complete'] = bool(checks_ok)
    rec['seconds'] = time.time() - t0
    _write_json(rec_path, rec)
    gc.collect()
    print(f"{wiring} {step_tag(step)}: done in {rec['seconds']:.0f} s; checks {'passed' if checks_ok else 'FAILED'}",
          flush=True)
    if not checks_ok:
        raise RuntimeError(f'{wiring} {step_tag(step)}: a check failed (see {rec_path})')
    return rec


def part_b(wirings: Sequence[str], steps: Sequence, root: Optional[Path] = None, force: bool = False) -> None:
    """Part B's points on `wirings` x `steps` (resumable; every point's checks as `run_point`'s; the identity step
    belongs to part 0). Wiring by wiring, each wiring's steps in ladder order."""
    for w in wirings:
        for st in steps:
            run_point(w, st, root, force=force)


# ---------------------------------------------------------------------------
# part 0

def base_population() -> dict:
    """The unchanged fits' populations (E83's `score_runs` on the saved fits, no simulation) under the calibrated readout
    against all recordings, against the entry's BASE_PATTERN (S 0.0500, T 0.0486, all8 0.0482; BASE_PATTERN_TOL)."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E86_reduced_input_candidate as E86
    rd, _ = E83.score_runs({}, tg=E86.targets())
    base = rd[(rd.change == 'none') & (rd.readout == 'cal')].set_index('unit')
    out = {u: dict(pattern=float(base.loc[u, 'z8_pattern']), entry=v, diff=float(base.loc[u, 'z8_pattern'] - v))
           for u, v in BASE_PATTERN.items()}
    raw = rd[(rd.change == 'none') & (rd.readout == 'raw')].set_index('unit')
    for u in BASE_PATTERN:
        out[u].update(raw_loss=float(raw.loc[u, 'z8_loss']), raw_mean_effect=float(raw.loc[u, 'mean_effect']),
                      cal_mean_effect=float(base.loc[u, 'mean_effect']))
    out['passes'] = bool(all(abs(out[u]['diff']) <= BASE_PATTERN_TOL for u in BASE_PATTERN))
    return out


def e98_flies() -> Tuple[Dict[str, pd.DataFrame], dict, Optional[dict]]:
    """E98's flies' ratios under every reading set (READING_SETS), from its saved a3_ratios.csv (the rows of its
    `flies_ratios`, which no part here recomputes), with its saved pulse test and solvent rates (verdict.json)."""
    from experiments import E98_apl_response_strength as E98
    res = H.shared_dir() / 'results' / E98.EXPERIMENT
    a3 = pd.read_csv(res / 'a3_ratios.csv', float_precision='round_trip')
    flies = {}
    for k, (rd, sub) in E98.READING_SETS.items():
        f = a3[(a3['reading'] == E98.reading_label(rd)) & (a3['subset'] == sub)].reset_index(drop=True)
        if f.empty:
            raise RuntimeError(f"E98's a3_ratios.csv has no rows for {k} ({rd}, {sub})")
        flies[k] = f
    v = json.loads((res / 'verdict.json').read_text())
    return flies, v['pulse'], v.get('flies_solvent_rate')


def e98_verdict(b1: pd.DataFrame) -> Tuple[dict, pd.DataFrame, pd.DataFrame]:
    """E98's verdict (`verdict_core`) of extended runs (`b1`, E98's `load_b1` rows) at the model's three pre-odor levels
    (`with_a0`), against E98's saved flies' ratios, pulse test and solvent rates (`e98_flies`)."""
    from experiments import E98_apl_response_strength as E98
    flies, pulse, solvent = e98_flies()
    levels = {lv: E98.with_a0(b1, lv) for lv in E98.A0_DEFINITIONS}
    return E98.verdict_core(flies, levels, pulse, solvent)


def e98_reading_identity() -> dict:
    """The no-simulation identity of E98's reading: `e98_verdict` on E98's own B1 runs (`load_b1_all`) gives E98's saved
    V1, V2 and V2's held class, and the carrying rows of c_comparisons.csv."""
    from experiments import E98_apl_response_strength as E98
    v, comps, _ = e98_verdict(E98.load_b1_all(ALL_WIRINGS))
    saved = json.loads((H.shared_dir() / 'results' / E98.EXPERIMENT / 'verdict.json').read_text())
    sc = pd.read_csv(H.shared_dir() / 'results' / E98.EXPERIMENT / 'c_comparisons.csv')
    key = ['comparison', 'reading_set']
    m = comps.merge(sc, on=key, suffixes=('', '_saved'))
    cols = ['flies_lo90', 'flies_hi90', 'model_linear_min', 'model_linear_max', 'model_min', 'model_max']
    dmax = float(max((m[c] - m[f'{c}_saved']).abs().max() for c in cols)) if len(m) else float('nan')
    cat_same = bool(len(m) and (m['category'] == m['category_saved']).all() and (m['category_linear'] == m['category_linear_saved']).all())
    out = dict(V1=v['V1'], V1_saved=saved['V1'], V2=v['V2'], V2_saved=saved['V2'], V2_held=v['V2_held'],
               V2_held_saved=saved['V2_held'], n_rows=int(len(comps)), n_matched=int(len(m)), max_abs_diff=dmax,
               categories_same=cat_same)
    out['passes'] = bool(out['V1'] == out['V1_saved'] and out['V2'] == out['V2_saved'] and out['V2_held'] == out['V2_held_saved']
                         and out['n_matched'] == len(sc) and cat_same and dmax <= 1e-9)
    return out


def e100_w1(b1: pd.DataFrame) -> dict:
    """E100's W1 (the six weak-odor ratios over paa at 10^-3, the model's input paa at 10^-4) of extended runs (`b1`)
    under each pre-odor level, against E100's saved flies' intervals (b_ratios.csv, set W, Yang's table, the carrying
    reading): per ratio the held class (`held_class`; carrying level and partner), the classes at each level, the linear
    reading's, and W1 (`w1_verdict`)."""
    from experiments import E100_apl_weak_odors as E100
    ratios = pd.read_csv(H.shared_dir() / 'results' / E100.EXPERIMENT / 'b_ratios.csv')
    stims = list(E100.W1_RATIOS) + [E100.W1_DENOMINATOR]
    matched = E100.match_b1(b1, stims, E100.session_inputs())
    den = E100.estimate_labels(matched, E100.W1_DENOMINATOR)
    d = (den.get(next((k for k, v in den.items() if v in H.APL16_ODORS), None)) or next(iter(den.values()))) if den else None
    if d is None:
        raise RuntimeError('no B1 stimulus for the W1 denominator')
    res, rows = {}, []
    for a0 in E100.A0_DEFINITIONS:
        on = E100.family_on(b1, a0)
        for st in E100.W1_RATIOS:
            f = E100._flies_row(ratios, 'W', st, E100.W1_DENOMINATOR)
            nums = E100.estimate_labels(matched, st)
            r, wc = E100.read_ratio(on, float(f['lo90']), float(f['hi90']), nums, d)
            res[(a0, E100.stim_label(st))] = r
            lin = wc['ratio_linear'] if len(wc) else pd.Series(dtype=float)
            rows.append(dict(a0=a0, stimulus=E100.stim_label(st), reference=E100.stim_label(E100.W1_DENOMINATOR),
                             flies_ratio=float(f['ratio']), flies_lo90=float(f['lo90']), flies_hi90=float(f['hi90']),
                             cls=r['cls'], cls_linear=r['cls_linear'], per_estimate=json.dumps(r['per_estimate']),
                             linear_min=float(lin.min()) if len(lin) else np.nan,
                             linear_max=float(lin.max()) if len(lin) else np.nan,
                             family_min=float(wc['ratio_min'].min()) if len(wc) else np.nan,
                             family_max=float(wc['ratio_max'].max()) if len(wc) else np.nan))
    held = {E100.stim_label(st): E100.held_class({a0: res[(a0, E100.stim_label(st))]['cls'] for a0 in E100.A0_DEFINITIONS})
            for st in E100.W1_RATIOS}
    held_lin = {E100.stim_label(st): E100.held_class({a0: res[(a0, E100.stim_label(st))]['cls_linear']
                                                      for a0 in E100.A0_DEFINITIONS}) for st in E100.W1_RATIOS}
    return dict(W1=E100.w1_verdict(held), classes=held, W1_linear=E100.w1_verdict(held_lin), classes_linear=held_lin,
                rows=pd.DataFrame(rows))


def e100_reading_identity() -> dict:
    """The no-simulation identity of E100's reading: `e100_w1` on E98's B1 runs gives E100's saved W1 and its classes."""
    from experiments import E98_apl_response_strength as E98
    from experiments import E100_apl_weak_odors as E100
    r = e100_w1(E98.load_b1_all(ALL_WIRINGS))
    saved = json.loads((H.shared_dir() / 'results' / E100.EXPERIMENT / 'verdict.json').read_text())
    out = dict(W1=r['W1'], W1_saved=saved['W1'], classes=r['classes'], classes_saved=saved['W1_classes'],
               W1_linear=r['W1_linear'], W1_linear_saved=saved['W1_linear_beside']['W1'])
    out['passes'] = bool(r['W1'] == saved['W1'] and r['classes'] == saved['W1_classes']
                         and r['W1_linear'] == saved['W1_linear_beside']['W1'])
    return out


def e77_identity(root: Optional[Path] = None) -> dict:
    """pn_inhadd 100 to 3141 on the real wiring against E77's screen, as part B's points recorded it
    (`e77_point_identity`); steps not run yet are absent."""
    out = {}
    for a in E77_S1:
        rec = load_point(a, 'real', root)
        if rec is not None and 'e77_identity' in rec:
            out[step_tag(a)] = rec['e77_identity']
    return out


def b4_identity(root: Optional[Path] = None) -> dict:
    """The extended runs at 1000 and 3141 on the real wiring against E98's B4, as part B's points recorded it
    (`b4_point_identity`), with their amended blank identities."""
    out = {}
    for a in E98_B4:
        rec = load_point(a, 'real', root)
        if rec is not None and 'b4_identity' in rec:
            out[step_tag(a)] = dict(rec['b4_identity'], a0_identity_passes=rec.get('a0_identity_passes'))
            out[step_tag(a)]['exact'] = bool(rec['b4_identity'].get('exact') and rec.get('a0_identity_passes'))
    return out


def part_0(wirings: Sequence[str] = ALL_WIRINGS, runs: bool = True, root: Optional[Path] = None) -> dict:
    """identity.json: the no-simulation identities (the unchanged fits' populations against BASE_PATTERN, E98's and
    E100's readings of their own saved runs), the identity point on `wirings` (`run_point` at STEP_IDENTITY, with the
    extended run) when `runs`, and E77's and E98's B4 identities where part B's runs exist. Earlier entries are kept."""
    t0 = time.time()
    ident = update_identity(base_population=base_population(), e98_reading=e98_reading_identity(),
                            e100_reading=e100_reading_identity())
    print(f"readers: base population {ident['base_population']['passes']}, E98 {ident['e98_reading']['passes']}, "
          f"E100 {ident['e100_reading']['passes']} ({time.time() - t0:.0f} s)", flush=True)
    if runs:
        for w in wirings:
            try:
                rec = run_point(w, STEP_IDENTITY, root)
                ok = bool(rec.get('complete'))
            except RuntimeError as e:
                rec, ok = load_point(STEP_IDENTITY, w, root) or {}, False
                print(f'{w}: identity point FAILED: {e}', flush=True)
            ident = update_identity(points={w: dict(exact=ok, **{k: rec.get(k) for k in (
                'identity_point', 'extended_identity', 'extended_matches_held', 'a0_identity_passes', 'seconds')})})
            if not ok:
                raise RuntimeError(f'{w}: the identity point failed; the experiment stops (identity.json)')
    ident = update_identity(e77=e77_identity(root), b4=b4_identity(root))
    provenance('0', wirings=list(wirings), runs=runs)
    return dict(base_population=ident['base_population']['passes'], e98_reading=ident['e98_reading']['passes'],
                e100_reading=ident['e100_reading']['passes'],
                points={w: p['exact'] for w, p in (ident.get('points') or {}).items()},
                e77={k: v.get('exact') for k, v in ident['e77'].items()}, b4={k: v.get('exact') for k, v in ident['b4'].items()})


# ---------------------------------------------------------------------------
# part A (no simulations)

def _flat_input(stimuli: pd.DataFrame) -> pd.DataFrame:
    """B1's stimulus table's input vectors (stimulus x the input's glomeruli)."""
    cols = [c for c in stimuli.columns if c.startswith('orn_')]
    x = stimuli.set_index('label')[cols].astype(float)
    x.columns = [c[len('orn_'):] for c in cols]
    return x


def _pn_named(pn: pd.DataFrame, glomeruli: Sequence[str]) -> pd.DataFrame:
    """A PN table indexed by glomerulus: as saved, or (E98's B2 to B4 runs, a RangeIndex 'pn_row') given the wiring's
    glomeruli in order."""
    if pn.index.name == 'glomerulus' or not isinstance(pn.index, pd.RangeIndex):
        return pn
    if len(pn) != len(glomeruli):
        raise RuntimeError(f'{len(pn)} PN rows for {len(glomeruli)} glomeruli')
    out = pn.copy()
    out.index = pd.Index(list(glomeruli), name='glomerulus')
    return out


def chain_table(runs: Mapping[str, Mapping[str, dict]], stimuli: pd.DataFrame,
                glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """A1's chain of totals per setting and stimulus (APL on): the ORN input's total, positive total, median glomerulus,
    broad component (the number of the input's glomeruli times the median, the median at least 0) and the rest (the
    positive excess over the median); the PN total over all PN rows and over the input's glomeruli, as the window total
    minus the blank's (the carrying level) and as window minus pre-odor mean (E90's); APL's window mean, the blank's,
    the settled and E90's pre-odor means, and the response at the carrying level. `runs`: {setting: {'on': E98's
    `load_b1_run` dict}}; `glomeruli`: the PN rows' glomeruli when saved without names (the real wiring's)."""
    from experiments import E98_apl_response_strength as E98
    X = _flat_input(stimuli)
    gl_in = list(X.columns)
    med = X.median(axis=1).clip(lower=0.0)
    rows = []
    for setting, run in runs.items():
        o = run['on']
        apl = o['apl']
        if E98.BLANK_LABEL not in apl.index:
            raise RuntimeError(f'{setting}: no blank in the run')
        pw = _pn_named(o['pn_win'].astype(float), glomeruli or [])
        pp = _pn_named(o['pn_pre'].astype(float), glomeruli or [])
        pin = [g for g in pw.index if str(g) in set(gl_in)]
        bw = pw[E98.BLANK_LABEL]
        b_apl = float(apl.loc[E98.BLANK_LABEL, 'win'])
        for s in apl.index:
            if s not in X.index or s not in pw.columns:
                continue
            x = X.loc[s].to_numpy(float)
            rows.append(dict(
                setting=setting, stimulus=str(s), input_total=float(x.sum()), input_total_pos=float(np.clip(x, 0, None).sum()),
                input_median=float(med[s]), input_broad=float(len(x) * med[s]),
                input_rest=float(np.clip(x - med[s], 0, None).sum()), n_input_glomeruli=int(len(x)),
                pn_total_blank=float((pw[s] - bw).sum()), pn_total_blank_input=float((pw[s] - bw).loc[pin].sum()),
                pn_total_e90=float((pw[s] - pp[s]).sum()), pn_total_e90_input=float((pw[s] - pp[s]).loc[pin].sum()),
                pn_win_total=float(pw[s].sum()), apl_win=float(apl.loc[s, 'win']), apl_blank=b_apl,
                apl_pre_settled=float(apl.loc[s, 'pre_settled']) if 'pre_settled' in apl.columns else np.nan,
                apl_pre_e90=float(apl.loc[s, 'pre']), apl_response=float(apl.loc[s, 'win'] - b_apl),
                kc_rate_on=float(o['responses'][s].astype(float).mean()) if s in o['responses'].columns else np.nan))
    return pd.DataFrame(rows)


def fixed_ratios() -> List[dict]:
    """A1's fixed list of ratios: E98's V1 (paa over each estimate of ms's input at 10^-3), V2 and its two beside, V3 to V5
    (Hallem's input; where the run has it); E100's W1 six over paa (the model's input at 10^-4, standing for 10^-3; the
    kiwi-control session's estimate, E49's source, and the other sessions' beside); E101's r1 (pa over ms near 10^-3, the
    panel's). Each with where the flies' interval comes from."""
    from experiments import E98_apl_response_strength as E98
    out = []
    for est, lab in E98.V1_ESTIMATES.items():
        out.append(dict(ratio_id='V1', estimate=est, numerator=E98.V1_NUM, denominator=lab, flies='E98:V1'))
    out.append(dict(ratio_id='V2', estimate='validation2/diag', numerator=E98.V2_NUM, denominator=E98.V2_DEN, flies='E98:V2'))
    for sn, lab in E98.V2_BESIDE.items():
        out.append(dict(ratio_id='V2 beside', estimate=f'{sn}/diag', numerator=lab, denominator=E98.V2_DEN, flies='E98:V2'))
    for rid, (num, den) in E98.HALLEM_COMPARISONS.items():
        out.append(dict(ratio_id=rid, estimate='hallem', numerator=num, denominator=den, flies=f'E98:{rid}'))
    for o, c in (('2h', -6), ('ma', -7), ('4h-ol', -6), ('farn', -2), ('2-but', -6), ('paa', -5)):
        for sn in ('kiwi-control', 'megamat', 'validation2'):
            out.append(dict(ratio_id=f'W1 {o} @ {c}', estimate=f'{sn}/diag', numerator=E98.session_label(o, c, sn, E98.DIAG),
                            denominator='paa', flies=f'E100:{o} @ {c}'))
    out.append(dict(ratio_id='r1 pa over ms', estimate='the panel', numerator='pa', denominator='ms', flies='E101:r1'))
    return out


def flies_interval(source: str) -> Optional[dict]:
    """The flies' APL ratio and 90 % interval for a fixed ratio's source ('E98:<ratio id>' from E98's a3_ratios.csv,
    the carrying reading, all recordings, V1 the hull of its two readings of paa's concentration; 'E100:<stimulus>' set W
    over paa at 10^-3 from E100's b_ratios.csv, Yang's table, the carrying reading; 'E101:r1' X2's APL, carrying reading,
    set L, from E101's c_x2_flies.csv)."""
    from experiments import E98_apl_response_strength as E98
    kind, key = source.split(':', 1)
    if kind == 'E98':
        flies, _, _ = e98_flies()
        f = flies['carrying'].set_index('ratio_id')
        ids = ('V1_paa-4', 'V1_paa-2.954') if key == 'V1' else (key,)
        if any(i not in f.index for i in ids):
            return None
        lo, hi = E98.hull([(float(f.loc[i, 'lo90']), float(f.loc[i, 'hi90'])) for i in ids])
        ratio = float(np.exp(np.mean([np.log(float(f.loc[i, 'ratio'])) for i in ids])))
        return dict(flies_ratio=ratio, flies_lo90=lo, flies_hi90=hi, flies_source=f'E98 a3_ratios.csv carrying {"/".join(ids)}')
    if kind == 'E100':
        from experiments import E100_apl_weak_odors as E100
        ratios = pd.read_csv(H.shared_dir() / 'results' / E100.EXPERIMENT / 'b_ratios.csv')
        r = E100._flies_row(ratios, 'W', E100.parse_stim(key), E100.W1_DENOMINATOR)
        if r is None:
            return None
        return dict(flies_ratio=float(r['ratio']), flies_lo90=float(r['lo90']), flies_hi90=float(r['hi90']),
                    flies_source='E100 b_ratios.csv set W over paa @ -3')
    if kind == 'E101':
        from experiments import E101_apl_vs_kc_claws as E101
        t = read_table(H.shared_dir() / 'results' / E101.EXPERIMENT / 'c_x2_flies.csv')
        t = t[(t['reading'] == 'carrying') & (t['b_mult'] == 0.0) & (t['ratio_id'] == key) & (t['quantity'] == 'APL')
              & (t['flies'] == 'L')]
        if t.empty:
            return None
        lo, hi = float(t['lo'].iloc[0]), float(t['hi'].iloc[0])
        return dict(flies_ratio=float(np.sqrt(lo * hi)), flies_lo90=lo, flies_hi90=hi,
                    flies_source='E101 c_x2_flies.csv APL r1 carrying, set L (ratio: the geometric middle of the interval)')
    raise ValueError(source)


def _log(x: float) -> float:
    return float(np.log(x)) if (np.isfinite(x) and x > 0) else float('nan')


def chain_ratios(chain: pd.DataFrame, ratios: Optional[Sequence[dict]] = None) -> pd.DataFrame:
    """A1's ratios: per setting and fixed ratio (`fixed_ratios`), the input's ratio (total, positive total, the rest
    above the broad component), the PN totals' (carrying level), APL's (linear, carrying level; through E98's indicator
    family beside), the flies' APL ratio with its interval (`flies_interval`); the power-law exponents (log APL ratio
    over log input ratio, the model's and the flies'; the antennal lobe's, log PN over log input; the mushroom body's, log
    APL over log PN) and the log split of the model's flattening (log APL ratio minus log input ratio) into the antennal
    lobe's part (log PN minus log input) and the mushroom body's (log APL minus log PN), each as a share."""
    from experiments import E98_apl_response_strength as E98
    ratios = fixed_ratios() if ratios is None else ratios
    fl = {}
    rows = []
    for setting, g in chain.groupby('setting', sort=False):
        g = g.set_index('stimulus')
        for r in ratios:
            num, den = r['numerator'], r['denominator']
            if num not in g.index or den not in g.index:
                continue
            a, b = g.loc[num], g.loc[den]
            if r['flies'] not in fl:
                fl[r['flies']] = flies_interval(r['flies']) or {}
            f = fl[r['flies']]
            row = dict(setting=setting, **r, **f)
            for k in ('input_total', 'input_total_pos', 'input_rest', 'pn_total_blank', 'pn_total_blank_input', 'apl_response'):
                row[f'{k}_num'], row[f'{k}_den'] = float(a[k]), float(b[k])
                row[f'{k}_ratio'] = float(a[k] / b[k]) if b[k] != 0 else np.nan
            if a['apl_response'] > 0 and b['apl_response'] > 0:
                rs = E98.model_ratio_set(a['apl_win'], b['apl_win'], a['apl_blank'])['ratio']
                row.update(apl_indicator_min=float(rs.min()), apl_indicator_max=float(rs.max()))
                if f:
                    row['model_class'] = E98.classify(rs.to_numpy(float), f['flies_lo90'], f['flies_hi90']) \
                        if f['flies_lo90'] > 0 else ''
            li, lp, la = _log(row['input_total_ratio']), _log(row['pn_total_blank_ratio']), _log(row['apl_response_ratio'])
            lr = _log(row['input_rest_ratio'])
            lf = _log(f.get('flies_ratio', np.nan))
            row.update(exponent_model=la / li if li else np.nan, exponent_flies=lf / li if li else np.nan,
                       exponent_al=lp / li if li else np.nan, exponent_mb=la / lp if lp else np.nan,
                       exponent_model_rest=la / lr if lr else np.nan, exponent_flies_rest=lf / lr if lr else np.nan,
                       flattening_log=la - li, al_part_log=lp - li, mb_part_log=la - lp,
                       al_share=(lp - li) / (la - li) if (la - li) else np.nan,
                       mb_share=(la - lp) / (la - li) if (la - li) else np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def low_conc_rates(responses: Mapping[str, pd.DataFrame], labels: Sequence[str]) -> Dict[str, float]:
    """The KCs' response rates (two or more spikes) at B1 stimuli `labels`, per condition: {'<label>|<cond>': rate}."""
    out = {}
    for ct, r in responses.items():
        for lab in labels:
            if lab in r.columns:
                out[f'{lab}|{ct}'] = float(r[lab].astype(float).mean())
    return out


def flies_low(stim: Tuple[str, float]) -> dict:
    """The flies' KC response rates at a low-concentration stimulus (`E100.kc_data_rates`: E72 part F's table)."""
    from experiments import E100_apl_weak_odors as E100
    return E100.kc_data_rates(stim)


def low_conc_table(runs: Mapping[str, Mapping[str, dict]], wiring: str = 'real') -> pd.DataFrame:
    """A2 (and part B's per wiring): per setting and LOW_STIMULI (and LOW_BESIDE) stimulus the KCs' response rates with
    APL on and acutely off, their ratio, against the flies' (control flies with APL on, APL output blocked (TNT) with APL
    off; `flies_low`), and the blank's response rates."""
    from experiments import E98_apl_response_strength as E98
    rows = []
    stims = list(LOW_STIMULI) + [(lab, ('ma', -7.0), 'beside') for lab in LOW_BESIDE]
    fl = {}
    for setting, run in runs.items():
        for lab, st, role in stims + [(E98.BLANK_LABEL, ('pfo', 0.0), 'blank')]:
            if lab not in run['on']['responses'].columns:
                continue
            on = float(run['on']['responses'][lab].astype(float).mean())
            off = float(run['acute_off']['responses'][lab].astype(float).mean())
            if st not in fl:
                fl[st] = flies_low(st) if st[0] != 'pfo' else {}
            f = fl[st]
            rows.append(dict(setting=setting, wiring=wiring, stimulus=lab, flies_stimulus=f'{st[0]} @ {st[1]:g}', role=role,
                             rate_on=on, rate_off=off, off_over_on=off / on if on > 0 else np.nan,
                             flies_on=f.get('kc_ctrl', np.nan), flies_off=f.get('kc_tnt', np.nan),
                             flies_on_natmix=f.get('kc_ctrl_natmix', np.nan),
                             flies_off_over_on=(f.get('kc_tnt', np.nan) / f.get('kc_ctrl', np.nan)) if f.get('kc_ctrl') else np.nan,
                             flies_n_on=f.get('kc_ctrl_n'), flies_n_off=f.get('kc_tnt_n'),
                             flies_conc=f.get('kc_data_conc')))
    return pd.DataFrame(rows)


A3_POINTS = ('base', 'S1_identity') + tuple(f'S1_inhadd{a:g}' for a in E77_S1)


def s1_parts() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: E77's S1 points (screen_pairs.csv's per-pair distances) under the calibrated readout and on raw counts against
    all recordings (the KC data's per-pair APL effect, `H.yang_apl_effect`): the loss, the mean APL effect, the pattern
    part with its amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair APL effect; and the per-odor
    part of the pattern part's error (`H.apl_effect_pair_parts`, its per-odor component averaged over each odor's 15
    pairs; a lead, not a result)."""
    from experiments import E72_wiring_population as E72
    from experiments import E77_input_screen as S
    res = H.shared_dir() / 'results' / S.EXPERIMENT
    pairs = pd.read_csv(res / 'screen_pairs.csv', float_precision='round_trip')
    screen = pd.read_csv(res / 'screen.csv').set_index('point_id')
    rows, per = [], []
    for pid in A3_POINTS:
        g = pairs[pairs['point_id'] == pid]
        if g.empty:
            continue
        idx = pd.MultiIndex.from_arrays([g['odor_a'], g['odor_b']], names=['odor_a', 'odor_b'])
        inhadd = PN_INHADD_DEFAULT if pid in ('base', 'S1_identity') else float(pid.split('inhadd')[1])
        for ro in ('cal', 'raw'):
            on = pd.Series(g[f'{ro}_on'].to_numpy(float), index=idx)
            off = pd.Series(g[f'{ro}_off'].to_numpy(float), index=idx)
            data = E72.data_for(ro)
            st = H.apl_effect_stats(on, off, data)
            sp = H.apl_effect_spread_shape(on, off, data)
            row = dict(point_id=pid, pn_inhadd=inhadd, readout=ro, loss=st['apl_effect_rmse'], mean_effect=st['apl_effect_mean'],
                       pattern=st['apl_effect_rmse_centered'], amplitude=sp['apl_effect_amplitude'], shape=sp['apl_effect_shape'],
                       sd_model=sp['apl_effect_sd_model'], sd_data=sp['apl_effect_sd_data'], spread_ratio=sp['apl_effect_sd_ratio'],
                       rho=st['apl_effect_spearman'], pearson=st['apl_effect_pearson'])
            sc = f'{ro}_pattern'
            if pid in screen.index and sc in screen.columns:
                row['screen_pattern'] = float(screen.loc[pid, sc])
                row['screen_loss'] = float(screen.loc[pid, f'{ro}_loss'])
            for k in ('rate_on', 'rate_off', 'rate_rho_on', 'rate_rho_off', 'fixed_thr', 'realism_flag'):
                if pid in screen.index and k in screen.columns:
                    row[k] = screen.loc[pid, k]
            rows.append(row)
            pp = H.apl_effect_pair_parts(on, off, data)
            for o in H.APL16_ODORS:
                m = [(a == o or b == o) for a, b in pp.index]
                per.append(dict(point_id=pid, pn_inhadd=inhadd, readout=ro, odor=o,
                                per_odor_mean=float(pp.loc[m, 'per_odor'].mean()),
                                error_mean=float(pp.loc[m, 'error'].mean()),
                                error_rms=float(np.sqrt((pp.loc[m, 'error'] ** 2).mean()))))
    return pd.DataFrame(rows), pd.DataFrame(per)


def _gain_spread(inhadd: float, inh: Tuple[float, float] = (1300.0, 5600.0)) -> float:
    """How much the PN gain varies over the typical range of the lateral inhibition (E3's note: 0.25 inhA + 0.75 inhB is
    about 1300 to 5600 during odors): log of the gain at the lower end over the gain at the upper end, pn_inhsc not
    entering. log(4.23) at the default, 0 with the gain fixed."""
    return float(np.log((inhadd + inh[1]) / (inhadd + inh[0])))


def prediction(a3: pd.DataFrame, ratios: pd.DataFrame, base: Mapping[str, dict]) -> pd.DataFrame:
    """A4: the prediction for every ladder step on S, written before part B. Each quantity's change on the real wiring
    (from C2) is read from E77's S1 points (the pattern part under the calibrated readout, the mean APL effect on raw
    counts) and from E98's B1 and B4 runs (V2's linear ratio on the real wiring), as a function of the normalization's
    strength (`_gain_spread`: 1.44 at the default, 0 with the gain fixed), linearly interpolated, and extrapolated beyond
    3141 linearly from 1000 and 3141; S's population is predicted as its unchanged fits' value plus the real wiring's
    change (the same change on every wiring)."""
    pts = a3['point_id'].str.startswith('S1_') | (a3['point_id'] == 'base')
    cal = a3[(a3['readout'] == 'cal') & pts]
    raw = a3[(a3['readout'] == 'raw') & pts]

    def curve(t: pd.DataFrame, col: str) -> Tuple[np.ndarray, np.ndarray]:
        t = t[t['point_id'] != 'S1_identity']
        x = np.array([_gain_spread(v) for v in t['pn_inhadd']])
        y = t[col].to_numpy(float)
        o = np.argsort(x)
        return x[o], y[o]

    def at(x: np.ndarray, y: np.ndarray, x0: float) -> float:
        if x0 >= x.min():
            return float(np.interp(x0, x, y))
        # beyond the screen's weakest normalization (3141): the line through its two weakest points (1000 and 3141)
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return float(y[0] + slope * (x0 - x[0]))

    xc, yc = curve(cal, 'pattern')
    xr, yr = curve(raw, 'mean_effect')
    c2_cal = float(cal[cal['point_id'] == 'base']['pattern'].iloc[0])
    c2_raw = float(raw[raw['point_id'] == 'base']['mean_effect'].iloc[0])
    v2 = ratios[ratios['ratio_id'] == 'V2'].set_index('setting')
    lo, hi = float(v2['flies_lo90'].iloc[0]), float(v2['flies_hi90'].iloc[0])
    xv = np.array([_gain_spread(PN_INHADD_DEFAULT), _gain_spread(1000.0), _gain_spread(3141.0)])
    yv = np.array([float(v2.loc[k, 'apl_response_ratio']) for k in ('C2', 'inhadd1000', 'inhadd3141')])
    o = np.argsort(xv)
    xv, yv = xv[o], yv[o]
    rows = []
    for a in (PN_INHADD_DEFAULT,) + LADDER:
        g = _gain_spread(a)
        p_real, m_real, v_real = at(xc, yc, g), at(xr, yr, g), at(xv, yv, g)
        rows.append(dict(step=step_tag(STEP_IDENTITY if a == PN_INHADD_DEFAULT else a), pn_inhadd=a, gain_spread_log=g,
                         in_screen=bool(a <= 3141.0), pattern_real_cal=p_real, d_pattern=p_real - c2_cal,
                         pattern_S_cal=float(base['S']['pattern'] + p_real - c2_cal), mean_effect_real_raw=m_real,
                         d_mean_effect_raw=m_real - c2_raw,
                         mean_effect_S_raw=float(base['S']['raw_mean_effect'] + m_real - c2_raw),
                         v2_linear_real=v_real, v2_inside_linear=bool(lo <= v_real <= hi), flies_lo90=lo, flies_hi90=hi))
    return pd.DataFrame(rows)


A_SETTINGS = (('C2', 'B1', 'real', None), ('inhadd1000', 'B4', 'E77_inhadd1000', 1000.0),
              ('inhadd3141', 'B4', 'E77_inhadd3141', 3141.0))


def part_a_runs() -> Tuple[Dict[str, dict], pd.DataFrame]:
    """Part A's saved runs on the real wiring: C2 (E98's B1), pn_inhadd 1000 and 3141 (E98's B4), from E98's local
    outputs or its archive; and B1's stimulus table (the input vectors)."""
    from experiments import E98_apl_response_strength as E98
    arc = H.shared_dir() / 'data' / 'model_outputs' / E98.EXPERIMENT / 'B'
    runs = {}
    for setting, arm, name, _ in A_SETTINGS:
        r = E98.load_b1_run(name, roots=[E98.b_root() / arm, arc / arm])
        if r is None:
            raise RuntimeError(f"E98's {arm} run {name} not found")
        runs[setting] = r
    stim = None
    for d in (E98.b_root() / 'B1', arc / 'B1'):
        if (d / 'b1_stimuli.csv').exists():
            stim = pd.read_csv(d / 'b1_stimuli.csv', float_precision='round_trip')
            break
    if stim is None:
        raise RuntimeError("E98's b1_stimuli.csv not found")
    return runs, stim


def part_a() -> dict:
    """Part A (no simulations): A1 to A4 (module docstring); writes a1_chain.csv, a1_ratios.csv, a2_low_conc.csv,
    a3_s1_parts.csv, a3_s1_per_odor.csv, a4_prediction.csv, a4_prediction.json and a_summary.json."""
    from experiments import E83_cross_wiring_fit as E83
    out = out_dir()
    runs, stim = part_a_runs()
    gl = list(dict.fromkeys(E83.wiring_parts('real')['W'].columns.get_level_values('glomerulus')))
    chain = chain_table(runs, stim, gl)
    chain.to_csv(out / 'a1_chain.csv', index=False)
    ratios = chain_ratios(chain)
    write_table(ratios, out / 'a1_ratios.csv', caveat=True)
    low = low_conc_table(runs)
    write_table(low, out / 'a2_low_conc.csv', caveat=True)
    a3, per = s1_parts()
    write_table(a3, out / 'a3_s1_parts.csv', caveat=True)
    write_table(per, out / 'a3_s1_per_odor.csv', caveat=True)
    ident = _read_json(identity_path()) or {}
    base = ident.get('base_population') or base_population()
    pred = prediction(a3, ratios, base)
    pred.to_csv(out / 'a4_prediction.csv', index=False)
    first = [r for r in pred.to_dict('records') if r['v2_inside_linear'] and r['step'] != STEP_IDENTITY]
    pj = dict(written=time.strftime('%Y-%m-%d %H:%M'), host=_host(), before_part_b=True,
              method=prediction.__doc__.split('\n')[0] + ' (see a4_prediction.csv)',
              steps=pred.to_dict('records'),
              l_lin_predicted_on_real=first[0]['step'] if first else None)
    _write_json(out / 'a4_prediction.json', pj)
    v2 = ratios[ratios['ratio_id'] == 'V2'][['setting', 'input_total_ratio', 'input_rest_ratio', 'pn_total_blank_ratio',
                                              'apl_response_ratio', 'flies_ratio', 'exponent_model', 'exponent_flies',
                                              'al_share', 'mb_share']]
    summ = dict(caveat=YANG_QUALITY_CAVEAT, v2=v2.to_dict('records'),
                low=low[low['stimulus'] == LOW_CARRYING][['setting', 'rate_on', 'rate_off', 'flies_on', 'flies_off']].to_dict('records'),
                s1=a3[a3['readout'] == 'cal'][['point_id', 'pattern', 'amplitude', 'shape', 'spread_ratio', 'mean_effect']].to_dict('records'),
                prediction=pj)
    _write_json(out / 'a_summary.json', summ)
    provenance('A', settings=[list(x) for x in A_SETTINGS])
    return summ


# ---------------------------------------------------------------------------
# score, part C's levels, part D's classes

def equal_rule(d_pattern: Mapping[str, float], d_raw_loss: float, d_mean: float, n_flagged: int, n_at_bound: int) -> dict:
    """E94's equal rule on S: 'equal' if the pattern part is within EQUAL_PATTERN against every EQUAL_TARGETS, no run is
    flagged or at a bracket's bound, the raw-count loss within EQUAL_RAW and the mean APL effect within EQUAL_MEAN;
    'better' (beside, E94's) if the pattern part is lower against all three and the rest holds."""
    missing = [t for t in EQUAL_TARGETS if t not in d_pattern or not np.isfinite(d_pattern[t])]
    if missing:
        raise ValueError(f'the pattern part against {missing} is missing')
    rest = bool(int(n_flagged) == 0 and int(n_at_bound) == 0 and abs(float(d_raw_loss)) <= EQUAL_RAW
                and abs(float(d_mean)) <= EQUAL_MEAN)
    eq = bool(all(abs(float(d_pattern[t])) <= EQUAL_PATTERN for t in EQUAL_TARGETS) and rest)
    return dict(equal=eq, lower_on_all=bool(all(float(d_pattern[t]) < 0 for t in EQUAL_TARGETS) and rest), rest=rest)


def _first(steps: pd.DataFrame, mask: pd.Series, largest: bool = False) -> Optional[float]:
    x = steps.loc[mask & ~steps['flagged'].astype(bool), 'value']
    if x.empty:
        return None
    return float(x.max() if largest else x.min())


def levels(steps: pd.DataFrame) -> dict:
    """Part C's levels from one row per ladder step (value pn_inhadd; 'v2' E98's held V2 verdict; 'v2_linear' the
    carrying reading's linear category; 'equal_S' the step's equal on S; 'flagged' any of its runs outside the bands):
    L_apl, L_lin, L_free (None where none exists; a flagged step is never a level), and the flagged steps."""
    s = steps.copy()
    for c in ('value', 'v2', 'v2_linear', 'equal_S', 'flagged'):
        if c not in s.columns:
            raise ValueError(f'levels: column {c!r} missing')
    s = s[s['value'].astype(float) > PN_INHADD_DEFAULT]
    return dict(L_apl=_first(s, s['v2'] == 'inside'), L_lin=_first(s, s['v2_linear'] == 'inside'),
                L_free=_first(s, s['equal_S'].astype(bool), largest=True),
                flagged=[float(v) for v in s.loc[s['flagged'].astype(bool), 'value']])


def refit_levels(lv: Mapping, steps: pd.DataFrame) -> List[float]:
    """The levels the refits run at: L_apl, or the last (largest) unflagged step when it does not exist (the note then
    says the ladder does not reach the flies' APL), and L_lin if it exists and differs."""
    out = []
    if lv.get('L_apl') is not None:
        out.append(float(lv['L_apl']))
    else:
        s = steps[(steps['value'].astype(float) > PN_INHADD_DEFAULT) & ~steps['flagged'].astype(bool)]
        if len(s):
            out.append(float(s['value'].max()))
    if lv.get('L_lin') is not None and float(lv['L_lin']) not in out:
        out.append(float(lv['L_lin']))
    return out


def refit_class(d_pattern: Mapping[str, float], d_raw_loss: float, d_mean: float, n_flagged: int) -> str:
    """'equal', 'better' or 'price' for S's refitted population against the unchanged fits: `d_pattern` the pattern part's
    change (refit minus unchanged, calibrated readout) against 'z8' (all recordings) and halves 'halfA' to 'halfD'. Equal:
    within EQUAL_PATTERN against z8, A and B, the four fits inside the bands, the raw-count loss within EQUAL_RAW and the
    mean APL effect within EQUAL_MEAN. Better: lower by BETTER_MARGIN or more against z8 and lower against each half
    (better is read first: an improvement larger than the equal's margin). A price: anything else; as the rule is written
    (E94's, the change within EQUAL_PATTERN in absolute value), that includes a refit lower by more than BETTER_MARGIN
    against all recordings that is not lower against every half, which the note then names as such."""
    need = ('z8',) + BETTER_HALVES
    missing = [t for t in need if t not in d_pattern]
    if missing:
        raise ValueError(f'the pattern part against {missing} is missing')
    if float(d_pattern['z8']) <= -BETTER_MARGIN and all(float(d_pattern[h]) < 0 for h in BETTER_HALVES) and int(n_flagged) == 0:
        return 'better'
    if equal_rule(d_pattern, d_raw_loss, d_mean, n_flagged, 0)['equal']:
        return 'equal'
    return 'price'


def t_better_holds(d_pattern_t: Mapping[str, float], d_single: Mapping[str, float], d_raw_loss_t: float) -> bool:
    """Whether 'better' stands on T (E83's gate thresholds): T's population lower than its fits by BETTER_MARGIN or more
    against all recordings ('z8') and lower against each half (BETTER_HALVES), at least T_BETTER_MIN_WIRINGS of the eight
    single wirings lower (`d_single`: wiring -> its pattern part's change against all recordings), and the raw-count loss
    not higher."""
    return bool(float(d_pattern_t['z8']) <= -BETTER_MARGIN and all(float(d_pattern_t[h]) < 0 for h in BETTER_HALVES)
                and sum(float(v) < 0 for v in d_single.values()) >= T_BETTER_MIN_WIRINGS and float(d_raw_loss_t) <= 0)


def low_conc_class(model: float, flies: float, model_c2: float, factor: float = LOW_CONC_FACTOR) -> str:
    """'further', 'nearer' or 'about the same' (outcome (e)): the distance from the flies is |log(model / flies)|; a step
    is 'further' if its distance exceeds C2's by log(factor) or more, 'nearer' if C2's exceeds it by that much. A model
    rate of 0 is infinitely far."""
    def dist(m):
        return float('inf') if m <= 0 else abs(np.log(m / flies))
    if not (flies > 0):
        raise ValueError(f"the flies' rate must be positive: {flies}")
    d, d0 = dist(float(model)), dist(float(model_c2))
    lf = np.log(float(factor))
    if np.isinf(d) and np.isinf(d0):
        return 'about the same'
    if d - d0 >= lf - 1e-12:
        return 'further'
    if d0 - d >= lf - 1e-12:
        return 'nearer'
    return 'about the same'


def _effect_parts(pairs: pd.DataFrame, tg: Mapping[str, pd.Series]) -> pd.DataFrame:
    """Per (change, unit, readout) of E83's `score_runs` pairs and per target: the pattern part's amplitude and shape
    parts, the spread ratio and Pearson r on the per-pair APL effect (`H.apl_effect_spread_shape`, `H.apl_effect_stats`)."""
    rows = []
    for (ch, unit, ro), g in pairs.groupby(['change', 'unit', 'readout'], sort=False):
        idx = pd.MultiIndex.from_arrays([g['odor_a'], g['odor_b']], names=['odor_a', 'odor_b'])
        on, off = pd.Series(g['on'].to_numpy(float), index=idx), pd.Series(g['off'].to_numpy(float), index=idx)
        row = dict(change=ch, unit=unit, readout=ro)
        for t, de in tg.items():
            sp = H.apl_effect_spread_shape(on, off, de)
            st = H.apl_effect_stats(on, off, de)
            row.update({f'{t}_amplitude': sp['apl_effect_amplitude'], f'{t}_shape': sp['apl_effect_shape'],
                        f'{t}_sd_ratio': sp['apl_effect_sd_ratio'], f'{t}_pearson': st['apl_effect_pearson']})
        rows.append(row)
    return pd.DataFrame(rows)


def _zcell_readings(runs: Mapping[Tuple[str, str], Dict[str, tuple]], tg: Mapping[str, pd.Series]) -> pd.DataFrame:
    """E73's per-KC reading under the saturating ('zcell') and the linear map ('zcell_lin'; E95's
    `zcell_distances_map`), per run and per population (E83's POPULATIONS whose wirings all have a run), scored as E83's
    `effect_scores`."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E95_static_apl as E95
    rows = []
    by = {}
    for (ch, w), o in runs.items():
        for ro, amp in (('zcell', 'saturating'), ('zcell_lin', 'linear')):
            by[(ch, w, ro)] = E95.zcell_distances_map(o, amp)
    changes = list(dict.fromkeys(ch for ch, _ in runs))
    for ch in changes:
        for ro in ('zcell', 'zcell_lin'):
            units = [(w, 'wiring', (w,)) for w in ALL_WIRINGS if (ch, w, ro) in by]
            units += [(n, 'population', ws) for n, ws in E83.POPULATIONS.items() if all((ch, w, ro) in by for w in ws)]
            for unit, kind, ws in units:
                d = by[(ch, unit, ro)] if kind == 'wiring' else E83.population([by[(ch, w, ro)] for w in ws])
                rows.append(dict(change=ch, unit=unit, kind=kind, n_wirings=len(ws), readout=ro, **E83.effect_scores(d, tg)))
    return pd.DataFrame(rows)


def available_steps(root: Optional[Path] = None) -> List[str]:
    """The step tags with at least one complete point under part B's folder, in ladder order."""
    base = (Path(root) if root is not None else local_dir()) / 'B'
    tags = [step_tag(STEP_IDENTITY)] + [step_tag(a) for a in LADDER]
    return [t for t in tags if any((_read_json(f) or {}).get('complete') for f in (base / t).glob('*_point.json'))]


def step_readings(tag: str, root: Optional[Path] = None, wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """One step's readings of its extended runs (no simulation): E98's verdict (`e98_verdict`) with its comparison rows,
    E100's W1 (`e100_w1`), E101's X2 on APL and Q for r1 and r3 (`E101.model_x2`, `x2_model_classes`, against E101's
    saved flies' intervals, set L, b 0), the KCs at low concentration per wiring (`low_conc_table`) and the PN stage
    against the PN dendrites (E90's `pn_statistics`, on the real wiring's run: the PN stage does not depend on the
    wiring)."""
    from experiments import E90_model_stages as E90
    from experiments import E98_apl_response_strength as E98
    from experiments import E101_apl_vs_kc_claws as E101
    step = step_of_tag(tag)
    d = (Path(root) if root is not None else local_dir()) / 'B' / tag
    have = [w for w in wirings if (d / f'{w}_on_apl.csv').exists()]
    out = dict(step=tag, wirings=have)
    if not have:
        return out
    b1 = E98.load_b1(have, roots=[d])
    v, comps, wrows = e98_verdict(b1)
    out['e98'] = v
    out['e98_comparisons'] = comps.assign(step=tag)
    out['e98_wirings'] = wrows.assign(step=tag)
    w1 = e100_w1(b1)
    out['e100'] = {k: w1[k] for k in ('W1', 'classes', 'W1_linear', 'classes_linear')}
    out['e100_rows'] = w1['rows'].assign(step=tag)
    fam = pd.concat([pd.read_csv(d / f'{w}_family.csv') for w in have if (d / f'{w}_family.csv').exists()], ignore_index=True)
    sets = E101.model_x2(fam, baselines=('blank', 'settled'))
    t = read_table(H.shared_dir() / 'results' / E101.EXPERIMENT / 'c_x2_flies.csv')
    t = t[(t['reading'] == 'carrying') & (t['b_mult'] == 0.0) & (t['flies'] == 'L')]
    x2 = []
    for rid in ('r1', 'r3'):
        for q in ('APL', 'Q'):
            f = t[(t['ratio_id'] == rid) & (t['quantity'] == q)]
            if f.empty:
                continue
            lo, hi = float(f['lo'].iloc[0]), float(f['hi'].iloc[0])
            by = {a0: E101.x2_model_classes(sets, rid, q, lo, hi, a0) for a0 in ('blank', 'settled')}
            held = by['blank']['cls'] if by['blank']['cls'] == by['settled']['cls'] else 'reading-dependent'
            x2.append(dict(step=tag, ratio_id=rid, quantity=q, flies_lo90=lo, flies_hi90=hi, cls_blank=by['blank']['cls'],
                           cls_settled=by['settled']['cls'], held=held,
                           per_estimate=json.dumps(by['blank']['per_estimate'])))
    s = E101.x2_model_summary(sets)
    out['e101'] = pd.DataFrame(x2)
    out['e101_summary'] = s.assign(step=tag)
    lows = []
    for w in have:
        run = E98.load_b1_run(w, roots=[d])
        lows.append(low_conc_table({tag: run}, w))
    out['low'] = pd.concat(lows, ignore_index=True)
    ref = 'real' if 'real' in have else have[0]
    run = E98.load_b1_run(ref, roots=[d])
    win, pre = run['on']['pn_win'].astype(float), run['on']['pn_pre'].astype(float)
    pn = E90.pn_statistics(win[list(H.APL16_ODORS)], pre[list(H.APL16_ODORS)])
    out['pn'] = dict(step=tag, pn_inhadd=step_value(step), wiring=ref, **pn)
    return out


def _step_table(rd: pd.DataFrame, realism: pd.DataFrame, apl: Mapping[str, dict], base_tag: str = 'none') -> pd.DataFrame:
    """One row per step for part C: pn_inhadd, V2's held verdict and its linear category, the equal on S at fixed scales
    (`equal_rule` against the unchanged fits, E83's change 'none'), the number of flagged runs and runs at a bound."""
    rows = []
    cal = rd[rd['readout'] == 'cal'].set_index(['change', 'unit'])
    raw = rd[rd['readout'] == 'raw'].set_index(['change', 'unit'])
    for tag, a in apl.items():
        if 'e98' not in a:
            continue
        r = realism[realism['step'] == tag]
        rs = r[r['wiring'].isin(S_WIRINGS)]
        comps = a['e98_comparisons']
        c2 = comps[comps['comparison'].str.startswith('V2, ') & (comps['reading_set'] == 'carrying')]
        row = dict(step=tag, value=step_value(step_of_tag(tag)), v2=a['e98']['V2'], v2_held=a['e98']['V2_held'],
                   v2_linear=str(c2['category_linear'].iloc[0]) if len(c2) else 'not read',
                   v2_linear_min=float(c2['model_linear_min'].iloc[0]) if len(c2) else np.nan,
                   v2_linear_max=float(c2['model_linear_max'].iloc[0]) if len(c2) else np.nan,
                   v1=a['e98']['V1'], w1=a['e100']['W1'], n_wirings=len(a['wirings']),
                   n_flagged=int((r['realism_flag'].fillna('').astype(str) != '').sum()),
                   n_flagged_S=int((rs['realism_flag'].fillna('').astype(str) != '').sum()),
                   n_at_bound=int((r['at_bound_thr'] | r['at_bound_apl']).sum()),
                   n_at_bound_S=int((rs['at_bound_thr'] | rs['at_bound_apl']).sum()))
        row['flagged'] = bool(row['n_flagged'] > 0)
        if (tag, 'S') in cal.index and (base_tag, 'S') in cal.index:
            dp = {t: float(cal.loc[(tag, 'S'), f'{t}_pattern'] - cal.loc[(base_tag, 'S'), f'{t}_pattern'])
                  for t in EQUAL_TARGETS + BETTER_HALVES}
            draw = float(raw.loc[(tag, 'S'), 'z8_loss'] - raw.loc[(base_tag, 'S'), 'z8_loss'])
            dme = float(cal.loc[(tag, 'S'), 'mean_effect'] - cal.loc[(base_tag, 'S'), 'mean_effect'])
            eq = equal_rule(dp, draw, dme, row['n_flagged_S'], row['n_at_bound_S'])
            row.update(equal_S=eq['equal'], **{f'd_pattern_S_{t}': v for t, v in dp.items()}, d_raw_loss_S=draw,
                       d_mean_effect_cal_S=dme, pattern_S=float(cal.loc[(tag, 'S'), 'z8_pattern']))
        else:
            row['equal_S'] = False
        rows.append(row)
    return pd.DataFrame(rows)


def part_score(root: Optional[Path] = None, out: Optional[Path] = None) -> dict:
    """The score (no simulation; from part B's saved runs, every step with a complete point): b_readings.csv (E83's
    `score_runs` per step, unit and readout against every target of E86's, with the amplitude and shape parts, the spread
    ratio and Pearson r, and E73's per-KC readings), b_realism.csv (each run's bands, E83's `realism`, and whether its
    bisections stopped at a bound), b_e98.csv / b_e98_wirings.csv (E98's comparisons), b_e100.csv (W1), b_e101.csv (X2's
    APL and Q), b_low_conc.csv, b_pn.csv, b_steps.csv (one row per step) and levels.json (part C). Refuses to run past a
    failed identity (the E77 and B4 identities are first read again from part B's points). `out`: the folder of the
    tables (default results/E102_antennal_lobe_normalization/)."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E86_reduced_input_candidate as E86
    out = out_dir() if out is None else Path(out)
    out.mkdir(parents=True, exist_ok=True)
    update_identity(e77=e77_identity(root), b4=b4_identity(root))
    tags = available_steps(root)
    runs, real_rows = {}, []
    for tag in tags:
        for w in ALL_WIRINGS:
            rec = load_point(step_of_tag(tag), w, root)
            o = load_held(step_of_tag(tag), w, root)
            if rec is None or o is None or not rec.get('complete'):
                continue
            runs[(tag, w)] = o
            real_rows.append(dict(step=tag, wiring=w, pn_inhadd=rec['pn_inhadd'], pn_inhsc=rec['pn_inhsc'],
                                  thr_factor=rec['thr_factor'], apl_factor=rec['apl_factor'],
                                  at_bound_thr=rec['at_bound_thr'], at_bound_apl=rec['at_bound_apl'],
                                  **{f'param_{k}': v for k, v in rec['params'].items()}, **E83.realism(o)))
    require_identity('all', sorted({w for _, w in runs}))
    tg = E86.targets()
    meta = {t: dict(step=t, pn_inhadd=step_value(step_of_tag(t))) for t in tags}
    rd, pairs = E83.score_runs(runs, meta, tg=tg)
    parts = _effect_parts(pairs, tg)
    rd = rd.merge(parts, on=['change', 'unit', 'readout'], how='left')
    zc = _zcell_readings(runs, tg)
    rd = pd.concat([rd, zc.assign(**{'step': zc['change']})], ignore_index=True) if len(zc) else rd
    write_table(rd, out / 'b_readings.csv', caveat=True)
    realism = pd.DataFrame(real_rows)
    write_table(realism, out / 'b_realism.csv', caveat=True)
    apl = {t: step_readings(t, root) for t in tags}
    for key, name in (('e98_comparisons', 'b_e98.csv'), ('e98_wirings', 'b_e98_wirings.csv'), ('e100_rows', 'b_e100.csv'),
                      ('e101', 'b_e101.csv'), ('e101_summary', 'b_e101_summary.csv'), ('low', 'b_low_conc.csv')):
        frames = [a[key] for a in apl.values() if key in a and len(a[key])]
        if frames:
            write_table(pd.concat(frames, ignore_index=True), out / name, caveat=True)
    pn = pd.DataFrame([a['pn'] for a in apl.values() if 'pn' in a])
    write_table(pn, out / 'b_pn.csv', caveat=True)
    steps = _step_table(rd, realism, apl)
    write_table(steps, out / 'b_steps.csv', caveat=True)
    lv = levels(steps[steps['n_wirings'] == len(ALL_WIRINGS)]) if len(steps) else {}
    lv.update(steps_scored=tags, complete_steps=[t for t in tags if any(
        (s['step'] == t and s['n_wirings'] == len(ALL_WIRINGS)) for s in steps.to_dict('records'))],
              refit_levels=refit_levels(lv, steps) if lv else [], rules=dict(
                  L_apl="smallest pn_inhadd with E98's V2 inside (6 of 8, carrying reading and level, held under partners)",
                  L_lin="smallest pn_inhadd with V2's linear reading inside on 6 of 8 (carrying reading and level)",
                  L_free="largest pn_inhadd that is an equal at fixed scales on S by E94's rule",
                  flagged='a step with a run outside the realism bands is not a level'),
              time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(out / 'levels.json', lv)
    if out == out_dir():
        provenance('score', steps=tags, n_runs=len(runs))
    return lv


def refit_commands(level: float, wirings: Sequence[str] = S_WIRINGS, root: Optional[Path] = None) -> List[str]:
    """Part D's commands: E10_joint_fit.py with C2's options (REFIT_ARGS, REFIT_MAXFEV evaluations) on each wiring at the
    level, from its point of part B (its held parameters as --init-params), with the step's PN parameters. E10_joint_fit.py
    does not take the PN parameters yet (the entry's "New: --pn-inhadd, Python only"); the commands carry
    `--pn-inhadd` / `--pn-inhsc` for when it does."""
    from experiments import E83_cross_wiring_fit as E83
    kws = step_kws(level)
    pfile = local_dir('D') / f'init_{step_tag(level)}.json'
    params = {}
    cmds = []
    for w in wirings:
        rec = load_point(level, w, root)
        if rec is None:
            raise RuntimeError(f'{w}: no part B point at {step_tag(level)}')
        params[w] = rec['params']
        wk = E83.wiring_kws(w)
        extra = []
        if 'wiring' in wk:
            extra = ['--wiring', wk['wiring'], '--wiring-seed', str(wk['wiring_seed'])]
        elif 'wiring_matrix' in wk:
            seed = wk['wiring_tag'].split('seed')[1]
            f = H.shared_dir() / 'data' / 'model_outputs' / 'E59_coinput_shuffle' / f'wPNKC_curveball_within_type-seed{seed}_matched.parquet'
            extra = ['--wiring-file', str(f), '--wiring-tag', wk['wiring_tag']]
        cmds.append(' '.join(['./venv/bin/python', '-u', 'experiments/E10_joint_fit.py', *REFIT_ARGS, '--maxfev',
                              str(REFIT_MAXFEV), '--experiment', EXPERIMENT, '--tag', f'refit_{step_tag(level)}_{w}',
                              '--pn-inhadd', repr(kws['pn_inhadd']), '--pn-inhsc', repr(kws['pn_inhsc']),
                              '--init-params', f'{pfile}:{w}', *extra]))
    _write_json(pfile, params)
    return cmds


# ---------------------------------------------------------------------------
# the CLI

def _wirings(arg: str) -> Tuple[str, ...]:
    if arg in WIRING_SETS:
        return WIRING_SETS[arg]
    w = tuple(x.strip() for x in arg.split(',') if x.strip())
    bad = [x for x in w if x not in ALL_WIRINGS]
    if bad:
        raise SystemExit(f'unknown wirings {bad}')
    return w


def _steps(arg: str) -> Tuple:
    if arg == 'ladder':
        return LADDER
    out = []
    for x in arg.split(','):
        x = x.strip()
        out.append(STEP_IDENTITY if x == STEP_IDENTITY else float(x))
    return tuple(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'score', 'D'])
    ap.add_argument('--wirings', default='all8', help="'S', 'T', 'all8' or a comma-separated list")
    ap.add_argument('--steps', default='ladder', help="'ladder' or comma-separated pn_inhadd values ('identity' allowed)")
    ap.add_argument('--root', default=None, help="part B: put the runs under another folder (default the local outputs)")
    ap.add_argument('--force', action='store_true', help='part B: run points that are saved again')
    ap.add_argument('--level', type=float, default=None, help='part D: the level (pn_inhadd) to print the refits for')
    ap.add_argument('--no-runs', action='store_true', help='part 0: only the identities that need no simulation')
    ap.add_argument('--out', default=None, help='score: write the tables to another folder (e.g. a smoke test)')
    a = ap.parse_args()
    root = Path(a.root) if a.root else None
    if a.part == '0':
        print(json.dumps(part_0(_wirings(a.wirings), runs=not a.no_runs, root=root), indent=1, default=str))
    elif a.part == 'A':
        print(json.dumps(part_a(), indent=1, default=str))
    elif a.part == 'B':
        steps = _steps(a.steps)
        if any(_is_identity(x) for x in steps):
            raise SystemExit('the identity step belongs to part 0')
        part_b(_wirings(a.wirings), steps, root, force=a.force)
    elif a.part == 'score':
        print(json.dumps(part_score(root, Path(a.out) if a.out else None), indent=1, default=str))
    elif a.part == 'D':
        if a.level is None:
            raise SystemExit('--part D needs --level')
        print('\n'.join(refit_commands(a.level, _wirings(a.wirings), root)))


if __name__ == '__main__':
    main()
