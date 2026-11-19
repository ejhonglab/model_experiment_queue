#!/usr/bin/env python3
"""E111: a residual APL output in the blocked flies in place of KC-KC inhibition: whether a weakened APL can do in the
best family what KC-KC does, and what would tell the two apart.

The entry is item E111 of experiment_priorities.md (proposed by viki 2026-11-18 15:21; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-18 14:49, continued 15:18; moved into the list 15:24 on Tom's word of 15:17,
relayed by hal). "The family" is C2's (the breadth-2 APL-effect fit of E67) on E72 part B's eight wirings, each at its own
fit (E83's `saved_fit`); "the family without KC-KC" is E110's part-B refit (`E110_no_kckc.counted_b`: KC-KC held at 0,
the other three parameters fitted, the lower of its starts); "the residue" is the fraction f of APL's output that the
blocked flies keep: the blocked condition with APL->KC and APL->bouton at f times their APL-on scales (E10_joint_fit.py's
`--apl-residual`, `apl_residual_params`; the inputs to APL, the threshold and KC-KC as APL on); "the residue family" is
part C's fits; "the two ρ's" are ρ on odor-odor distances with APL on against the control flies and with APL off against
APL output blocked (TNT), on raw counts (`E106_rho_gap_parts.two_rhos`). The family's fits and E110's part-B fits are
references and are not refitted. Below, "APL off" and "the blocked condition" of a run with a residue are its blocked
condition at its f.

Every fit is E72's / E67's part-B command, `E10_joint_fit.py -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl
--apl-readout both --tag _cap3_<part>` (150 evaluations), with the lab's input and the wiring's options (E110's
`wiring_args`), run as its own process into its own folder (results and outputs of `--experiment E111_residual_apl/fits/
<name>`), as E110's are.

Parts (CLI `--part`; `--wirings` takes a subset; `--lane` names the lane's tables, default the host's short name; a point
is claimed by `mkdir` of results/E111_residual_apl/claims/<name>/ and skipped when done or claimed, `--reclaim` taking over
a claim left by a process that died; `--wait` retries every WAIT_SECONDS points whose prerequisites another lane or host
is still making; `--part commands` prints every command):
- 0, the option's two identities on the model, each one E10_joint_fit.py process of one evaluation (`--maxfev 1`, so that
  its fits.csv row and its saved runs are written) from C2's exact start on the hemibrain wiring (`--init-from-exact`,
  KC-KC free as in the family): f0_real (`--apl-residual 0`): its first history row equals C2's history at its best
  evaluation on every column they share (E110's `exact_check`), its fits.csv row equals C2's (every `fit_<term>` and
  `loss` exactly; every other shared column but the run's own metadata, META_SKIP, compared and listed), and its saved
  APL-on and blocked responses and spike counts equal C2's saved ones exactly; f1_real (`--apl-residual 1`): its saved
  blocked run equals its saved APL-on run exactly, and that equals C2's. check.json; `--checks` runs a subset (f0_real
  first; f1_real may follow later). The fits and simulations below need f0_real passed and no check failed.
- B (simulations): the dose at fixed scales. E110's part-B optimum of each wiring (its parameters and its saved APL-on
  run), the blocked condition at f = 0.01, 0.02, 0.05, 0.1 and 0.2, nothing re-held (`E83.build_simulator(w, params)`;
  on the hemibrain wiring first f = 0, which must reproduce E110's saved blocked run exactly). Waits for E110's counted
  part-B fit of the wiring and its saved runs (this host's outputs/E110_no_kckc/fits/ or the archive). Rows in
  runs_B_<lane>.csv; runs under outputs/E111_residual_apl/B/ (tag B<f>_<w>, the APL-on run E110's copied beside).
- C (fits): the residue fitted. KC-KC held at 0 (`--kckc-fixed 0`); the somatic threshold, APL->KC (KC->APL and
  bouton->APL tied), APL->bouton (cap 3 times the base's) and f (`--apl-residual 0.05 --fit-apl-residual`, bounds 0 to
  0.5) fitted; started at E110's part-B optimum (starts.json:C_<w>, all five scales); a second start C2 at f = 0.2 on real,
  P1, P2 and M1. Of two starts the lower `fit_apl_cal_rmse` counts (E110's rule). C on P3, P4, M2 and M3 runs only when
  E110's verdict is final and "needed" or "not determined" (the entry: otherwise part C runs on real, P1, P2, M1 only and
  gives no class); `--force` runs them anyway.
- D (fits): both free, on real, P1, P2 and M1: KC-KC and f fitted together (`--apl-residual 0.02 --fit-apl-residual`,
  KC-KC free) from the family's fit (`--init-from-exact`). Does not wait for E110.
- E (simulations, after C): what each family predicts, no fit. Per wiring, with f_C the residue family's fitted f: the
  family at its fit with the blocked condition at f_C (Efam_<w>), the family without KC-KC likewise (Eb_<w>), and the
  residue family with its residue removed (f = 0, Ec_<w>; on the hemibrain wiring first its blocked run at f_C, which
  must reproduce part C's saved one exactly). Rows in runs_E_<lane>.csv; runs under outputs/E111_residual_apl/E/.
- score (no simulations): runs.csv, readings.csv (E110's readings per kind, unit and readout, with the block's
  pair-specific remainder), changes_C_family.csv and changes_C_B.csv (the verdict's measures), bands_C.csv, e_table.csv
  (part E per family and wiring and as the population), told_apart.csv, d_reading.csv, b_dose.csv, data_intervals.json,
  verdict.json, comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive and kcq/apl_effect.py.
- commands: every fit's E10_joint_fit.py command and the lanes' lines; starts: (re)writes starts.json.

The verdict (fixed in the entry before any run; part C, the population and the mean paired difference over the eight
wirings, E110's margins; `verdict`, `told_apart`):
- "a residue replaces KC-KC": against the family, every condition of E110's "not needed" holds (the pattern part under
  the calibrated readout not higher by more than 0.003, the raw-count loss not higher by more than 0.005, the mean APL
  effect not lower by more than 0.01 under either reading, each ρ not lower by more than 0.005, response rate ρ with APL
  on and with APL off not lower by more than 0.05: each on the population, the mean paired difference and 6 of the 8
  wirings; every run inside the four bands), and f not at a bound;
- "a residue helps but does not replace": some of these fail, and against the family without KC-KC the pattern part is
  lower by 0.003 or more, or a ρ higher by 0.010 or more, on the population and on 6 of the 8 wirings (E110's gains);
- "a residue does nothing": neither, or the fitted f is below 0.01 on at least 6 of the 8 wirings;
- "told apart by the recordings on hand": on at least 6 of the 8 wirings, on one of part E's measures of the blocked
  condition (TOLD_MEASURES), the flies' 90 % interval over recordings contains one family's value and excludes the
  other's (the residue family against the family; against the family without KC-KC reported beside). Otherwise "not told
  apart", and part E's removal (the change in the blocked condition's response rate and mean pair distance when the
  residue is removed) is the prediction for a fly experiment;
- if E110's verdict is "not needed" or "better without": part C runs on real, P1, P2 and M1 only and gives no class;
  parts B, D and E are reported;
- reported apart (pending Tom's answers of 2026-11-17 14:44): the block, the natural-source alignment on, off and on
  minus off, the structure against wild-type somata (E110's `pending_table`). E106's yardsticks hold as in E110. No
  correlation chooses anything, no pair is weighted, and no threshold moves after a run is read.

Deviations and choices (fixed here before any run):
- The option is a condition hook in E10_joint_fit.py (`fit_one`): optimize_mb_model_params' overrides are static, and a
  residue follows the fitted APL-on scales. With the option off the call is unchanged (test/test_e111_residual_apl.py);
  part 0 checks f = 0 and f = 1 on the model. The breadth run and its pre-passes never take the option (APL's output a
  true 0 there), so the breadth weights are the family's.
- f is fitted linearly within [0, 0.5] (optimize_mb_model_params' rule for a lower bound of 0: a log scale needs a
  positive bound). "f at a bound": within 1 % of the range (F_BOUND_TOL, 0.005) of 0 or of 0.5; "f not at a bound" holds
  when at most 2 of the 8 wirings are at a bound (the 6-of-8 share). "f below 0.01": strictly below F_NOTHING.
- E110's class (`e110_class`): its verdict once final, or its provisional verdict when that is "needed" or "not
  determined" (E110's provisional class changes only to "not determined", when its part C's return fits are flat, so
  either way part C runs on all eight wirings and gives a class); a provisional "not needed" or "better without" waits
  for the final verdict (part C on P3, P4, M2 and M3 and the class wait; `--force` runs those fits anyway).
- Part C's second start (C2) has the 150 evaluations of the entry's command (E110's second start had 200 by its entry);
  part D likewise.
- "The residue family" is part C's counted fit per wiring; part E reads it per wiring at its own f_C. "The step from f to
  0" is, for each family, its blocked condition at f_C minus its blocked condition at 0, pair by pair (E76's step is the
  previous point's distance minus the next's, so a step toward less APL is positive): for the residue family part C's
  blocked run minus Ec's; for the family Efam's minus its fit's; for the family without KC-KC Eb's minus E110's. Its ρ
  over the 120 pairs is against the recordings' step from the six typical recordings to the two strongest (E76's
  `_steps`, step 2: the six's mean distance minus the two's) in Yang's z and in baseline-noise units, on raw counts and
  under the calibrated readout with K held at the run's APL-on calibration; the two recordings' own agreement (E76's
  rho_own) beside. "The change if the residue were removed" is the blocked condition at 0 minus at f_C (response rate,
  mean pair distance raw and under the readout), for each family with the same f_C.
- The measures of "told apart" (TOLD_MEASURES) are those of part E with a value of the flies' own: the block (the
  registry's 90 % interval of `block`), the block's pair-specific remainder (the mean over the nine block pairs of the
  isotonic residual beyond additive odor effects fitted on the other 111 pairs, in cosine similarity, as
  analyses/yang_pair_deviations.py computes it, with its 90 % interval from results/yang_pair_deviations/summary.csv),
  the APL-off response rate (two spikes in the model; the flies' rates in baseline-noise units per recording, E76's
  recordings.csv, the 90 % interval of their mean over 4000 resamples of the eight recordings with APL output blocked
  (TNT), seed 0) and odors per KC with APL off (the total variation distance from the flies' mean distribution, against
  the band, the most atypical recording's distance: the flies' recordings lie inside it by construction, so the interval
  is [0, band]). Response rate ρ with APL off is reported but not a "told apart" measure (it is an agreement with the
  flies, with no value of the flies' own). The APL-off / APL-on response rate ratio with the registry's interval
  (noise_rate_ratio_fraction) is reported beside.
- Part D's reading (outcome (e)): per wiring, the fit "keeps KC-KC" when its KC-KC ends above 10 (E110's return
  threshold) and "uses a residue" when its f ends at 0.01 or above; reported with counts, choosing nothing.
- check_passed: the fits and simulations need part 0's f0_real passed and no check failed (f1_real may be run later on
  a loaded machine; check.json then says it is pending).
- The fits' and the simulations' child processes get this checkout's root first on PYTHONPATH (`_child_env`), so a run
  from a worktree loads the worktree's code; the venv's editable install points at the main checkout otherwise. The
  fits are launched with this process's interpreter (sys.executable).

Provenance: as E110 (`E108._commits_clean`, `E110._fit_clean`; each row records `H.code_commits()` with
'queue_dirty_paths'; provenance_<part>_<lane>.json has `H.model_provenance` per point; each fit's folder has
E10_joint_fit.py's own provenance and the launch record launch.json). The holds (CLAUDE.md, "Connectome data (Tom,
2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of it; no other connectome table; no raw-only
recording. Nothing is published to claude.ai.

Usage (from the repository root; threads = odors, 16, or fewer on a shared machine; under nice). Part 0 runs once (its
check.json is shared); then lanes on any hosts (the claims split the points):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E111_residual_apl.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E111_residual_apl.py --part D --lane hal1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E111_residual_apl.py --part B --lane hal1 --wait
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E111_residual_apl.py --part C --lane hal1 --wait
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E111_residual_apl.py --part E --lane hal1 --wait
    ./venv/bin/python -u experiments/E111_residual_apl.py --part commands
then, on every host that ran a part:
    ./venv/bin/python archive_tools/archive_model_outputs.py E111_residual_apl --exclude 'check,distances_cache.pkl,fits/*/*base_acute_off_spike_counts.parquet'
    ./venv/bin/python archive_tools/migrate_model_output_params.py E111_residual_apl --apply
and, on any host:
    ./venv/bin/python -u experiments/E111_residual_apl.py --part score
    ./venv/bin/python -u experiments/E111_residual_apl.py --part variants
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import subprocess
import sys
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E10_joint_fit as J
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108
from experiments import E110_no_kckc as E110

EXPERIMENT = 'E111_residual_apl'
WIRINGS = E110.WIRINGS
S4 = E110.S_WIRINGS                                # real, P1, P2, M1
T4 = E83.T_WIRINGS                                 # P3, P4, M2, M3
C2_WIRINGS = S4
D_WIRINGS = S4
F_DOSE = (0.01, 0.02, 0.05, 0.1, 0.2)
C_F0, C2_F0, D_F0 = 0.05, 0.2, 0.02
F_BOUNDS = J.APL_RESIDUAL_BOUNDS                   # (0, 0.5)
F_NOTHING = 0.01
F_BOUND_TOL = 0.01 * (F_BOUNDS[1] - F_BOUNDS[0])
D_KCKC_KEEP = E110.C_RETURN_KCKC                   # 10
K = E110.K
MAXFEV = 150
CHOOSE_COL = E110.CHOOSE_COL
PARAM_KEYS = E110.PARAM_KEYS
READOUTS = E110.READOUTS
UNITS = dict(E110.UNITS)
N_WIRINGS = E110.N_WIRINGS                         # 6 of the 8
E110_NO_CLASS = ('not needed', 'better without')
E110_RUN_ALL = ('needed', 'not determined')
CLASSES = ('a residue replaces KC-KC', 'a residue helps but does not replace', 'a residue does nothing')
TOLD = ('told apart by the recordings on hand', 'not told apart')
FIT_EXPERIMENT_PREFIX = f'{EXPERIMENT}/fits'
CHECK_EXPERIMENT_PREFIX = f'{EXPERIMENT}/check'
WAIT_SECONDS = 300
# part 0: the run's own metadata in its fits.csv row (not compared with C2's row)
META_SKIP = ('tag', 'variant', 'maxfev', 'n_evals', 'seconds', 'initial_loss', 'init_from', 'init_from_experiment',
             'init_from_exact', 'init_params', 'init_thr', 'kckc0', 'thr0_factor', 'method', 'apl_residual',
             'apl_residual_fitted', 'apl_residual_init', 'apl_residual_lo', 'apl_residual_hi', 'off_wAPLKC_scale',
             'off_wAPLPN_scale')
# the bootstrap of the flies' APL-off response rate
RATE_BOOT, RATE_BOOT_SEED = 4000, 0

# part E's measures of the blocked condition that have a value of the flies' own: key, label, e_table column
ToldMeasure = namedtuple('ToldMeasure', 'key label col')
TOLD_MEASURES = (
    ToldMeasure('block', 'the block (raw counts)', 'block'),
    ToldMeasure('block_ps', "the block's pair-specific remainder (raw counts, cosine similarity)", 'block_ps'),
    ToldMeasure('rate_off', 'the APL-off response rate', 'rate_off'),
    ToldMeasure('opk_off_tv', 'odors per KC with APL off (distance from the flies)', 'opk_off_tv'),
)
FAMILY_KINDS = {'family': 'the family', 'B110': 'the family without KC-KC', 'C': 'the residue family'}

_CACHE: Dict[object, object] = {}

FitSpec = namedtuple('FitSpec', 'name part wiring kckc_fixed f0 maxfev start suffix')
# start: ('params', <starts.json key>) or ('exact',) (the wiring's family fit, --init-from-exact)


# ---------------------------------------------------------------------------
# folders and small helpers

def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    return H.outputs_dir(f'{EXPERIMENT}/{part}' if part else EXPERIMENT)


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


_write_json, _read_json, _host, _now = E110._write_json, E110._read_json, E110._host, E110._now


def _python() -> str:
    return sys.executable


def _child_env() -> dict:
    """E110's child environment (OMP_THREAD_LIMIT set to OMP_NUM_THREADS when only that is set, unbuffered), with this
    checkout's root first on PYTHONPATH, so that the child loads the code of the checkout this process runs from (a
    worktree's; the venv's editable install points at the main checkout)."""
    env = E110._child_env()
    pp = [p for p in env.get('PYTHONPATH', '').split(os.pathsep) if p]
    env['PYTHONPATH'] = os.pathsep.join([str(H.HERE)] + [p for p in pp if Path(p).resolve() != H.HERE.resolve()])
    return env


def residual_spec(f0: float, fitted: bool) -> dict:
    return J.apl_residual_spec(float(f0), fitted)


# ---------------------------------------------------------------------------
# the fits' commands

def fit_tag(spec: FitSpec, rate: Optional[float] = None) -> str:
    """The tag E10_joint_fit.py gives a fit of this experiment (its prefix rules: the wiring, rate<R>_, aplboth_,
    kckc<S>_, aplresfit<f0>_, then k2 and the --tag suffix)."""
    r = float(H.data_target('ctrl_frac')) if rate is None else float(rate)
    kc = '' if spec.kckc_fixed is None else f'kckc{float(spec.kckc_fixed):g}_'
    return (f'{E110.wiring_prefix(spec.wiring)}rate{r:.3f}_aplboth_{kc}{J.apl_residual_tag(residual_spec(spec.f0, True))}'
            f'k{float(K):g}_cap3_{spec.suffix}')


def fit_specs(part: str, wirings: Optional[Sequence[str]] = None, with_t4: bool = True) -> List[FitSpec]:
    """The fits of a part ('C' or 'D'), in order; `wirings` keeps a subset; `with_t4` False drops part C's fits on P3, P4,
    M2 and M3 (run only on E110's verdict "needed" or "not determined")."""
    out: List[FitSpec] = []
    if part == 'C':
        for w in WIRINGS:
            if w in T4 and not with_t4:
                continue
            out.append(FitSpec(f'C1_{w}', 'C1', w, 0.0, C_F0, MAXFEV, ('params', f'C_{w}'), 'C1'))
            if w in C2_WIRINGS:
                out.append(FitSpec(f'C2_{w}', 'C2', w, 0.0, C2_F0, MAXFEV, ('params', f'C_{w}'), 'C2'))
    elif part == 'D':
        out = [FitSpec(f'D_{w}', 'D', w, None, D_F0, MAXFEV, ('exact',), 'D') for w in D_WIRINGS]
    else:
        raise ValueError(f'no fits in part {part!r}')
    return [s for s in out if wirings is None or s.wiring in wirings]


def starts_path() -> Path:
    return out_dir() / 'starts.json'


def starts_spec(key: str) -> str:
    return f'results/{EXPERIMENT}/starts.json:{key}'


def fit_command(spec: FitSpec, python: Optional[str] = None) -> List[str]:
    """The E10_joint_fit.py command of one fit (E72's part-B command, the residue fitted from spec.f0, this fit's start),
    into its own folders."""
    cmd = [python or _python(), '-u', 'experiments/E10_joint_fit.py', '-k', f'{float(K):g}', '--aplpn-max', '3',
           '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_cap3_{spec.suffix}',
           '--maxfev', str(int(spec.maxfev)), '--experiment', f'{FIT_EXPERIMENT_PREFIX}/{spec.name}'] + E110.wiring_args(spec.wiring)
    if spec.kckc_fixed is not None:
        cmd += ['--kckc-fixed', repr(float(spec.kckc_fixed))]
    cmd += ['--apl-residual', repr(float(spec.f0)), '--fit-apl-residual']
    if spec.start[0] == 'params':
        cmd += ['--init-params', starts_spec(spec.start[1])]
    elif spec.start[0] == 'exact':
        folder, tag = E83.members()[spec.wiring]
        cmd += ['--init-from', tag, '--init-from-experiment', folder, '--init-from-exact']
    else:
        raise ValueError(f'unknown start {spec.start}')
    return cmd


def fit_results_dir(name: str, check: bool = False) -> Path:
    return H.shared_dir() / 'results' / (CHECK_EXPERIMENT_PREFIX if check else FIT_EXPERIMENT_PREFIX) / name


def fit_row(name: str, tag: str) -> Optional[pd.Series]:
    p = fit_results_dir(name) / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    r = f[f['tag'] == tag]
    return r.iloc[-1] if len(r) else None


def load_fit_outputs(name: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A fit's saved runs: this host's outputs/E111_residual_apl/fits/<name>/, else the archive's."""
    for d in (H.HERE / 'outputs' / FIT_EXPERIMENT_PREFIX / name, archive_dir() / 'fits' / name):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def params_of_row(r: Mapping) -> Dict[str, float]:
    return E110.params_of_row(r)


def f_of_row(r: Mapping) -> float:
    return float(r['apl_residual'])


# ---------------------------------------------------------------------------
# starts, claims, launches

def c_start(row: Mapping, name: str) -> dict:
    """Part C's start: E110's counted part-B optimum (a fits.csv row; KC-KC 0), all five scales as fitted."""
    p = params_of_row(row)
    if float(p['wKCKC_scale']) != 0.0:
        raise ValueError(f"E110's part-B fit {name} has KC-KC {p['wKCKC_scale']}")
    return E110.start_entry(p, source=f"E110's counted part-B optimum {name}:{row['tag']} (results/E110_no_kckc/fits/{name})")


def write_starts() -> dict:
    """starts.json from what is saved: C_<w> for every wiring whose E110 part-B counted fit is done."""
    entries = {}
    for w in WIRINGS:
        b = E110.counted_b(w)
        if b is not None:
            entries[f'C_{w}'] = c_start(b['row'], b['name'])
    return E110.merge_starts(starts_path(), entries)


def claim(name: str, lane: str, reclaim: bool = False) -> bool:
    d = out_dir() / 'claims' / name
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(name=name, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(name: str) -> None:
    shutil.rmtree(out_dir() / 'claims' / name, ignore_errors=True)


def launch(cmd: List[str], res: Path, log: Path, cc: dict, name: str) -> dict:
    """Runs one E10_joint_fit.py process (log `log`), with its launch record res/launch.json (as E110's `launch`)."""
    env = _child_env()
    rec = dict(name=name, command=cmd, cwd=str(H.HERE), host=_host(), started=_now(), commits=cc,
               queue_dirty_paths=str(cc.get('queue_dirty_paths') or '').split(), log=str(log.relative_to(H.HERE)),
               threads=dict(OMP_NUM_THREADS=env.get('OMP_NUM_THREADS'), OMP_THREAD_LIMIT=env.get('OMP_THREAD_LIMIT')),
               pythonpath=env.get('PYTHONPATH'), launcher='experiments/E111_residual_apl.py')
    res.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    _write_json(res / 'launch.json', rec)
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT, env=env)
    rec.update(finished=_now(), returncode=r.returncode, seconds=time.time() - t0)
    _write_json(res / 'launch.json', rec)
    return rec


def runs_table_of(kind: str) -> pd.DataFrame:
    """Every lane's rows of a kind ('B', 'E' or 'fits'): runs_<kind>_<lane>.csv, last row per tag."""
    ps = sorted(out_dir().glob(f'runs_{kind}_*.csv'))
    if not ps:
        return pd.DataFrame(columns=['tag', 'wiring', 'name'])
    return pd.concat([E108._read_rows(p) for p in ps], ignore_index=True, sort=False).drop_duplicates('tag', keep='last')


# ---------------------------------------------------------------------------
# part 0: the option's two identities

CHECKS = ('f0_real', 'f1_real')


def check_specs() -> Dict[str, Tuple[str, List[str]]]:
    """Part 0's two E10_joint_fit.py runs (one evaluation each) on the hemibrain wiring from C2's exact start: {name: (tag,
    command)}."""
    folder, tag = E83.members()['real']
    out = {}
    for name, f, suffix in (('f0_real', 0.0, 'check0'), ('f1_real', 1.0, 'check1')):
        cmd = [_python(), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3', '--rate-target', 'ctrl_frac',
               '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_cap3_{suffix}', '--maxfev', '1', '--experiment',
               f'{CHECK_EXPERIMENT_PREFIX}/{name}', '--init-from', tag, '--init-from-experiment', folder, '--init-from-exact',
               '--apl-residual', repr(f)]
        r = float(H.data_target('ctrl_frac'))
        out[name] = (f'rate{r:.3f}_aplboth_{J.apl_residual_tag(residual_spec(f, False))}k2_cap3_{suffix}', cmd)
    return out


def _check_local(name: str) -> Path:
    return H.HERE / 'outputs' / CHECK_EXPERIMENT_PREFIX / name


def same_outputs(a: Dict[str, tuple], b: Dict[str, tuple], ca: str, cb: str) -> dict:
    """Condition `ca` of `a` against `cb` of `b`: index and columns equal, and the number of responses and spike counts
    that differ (0 and 0: exact)."""
    (ra, sa), (rb, sb) = a[ca], b[cb]
    shape = bool(sa.index.equals(sb.index) and sa.columns.equals(sb.columns) and ra.index.equals(rb.index)
                 and ra.columns.equals(rb.columns))
    ns = int((sa.to_numpy(dtype=float) != sb.to_numpy(dtype=float)).sum()) if sa.shape == sb.shape else -1
    nr = int((ra.to_numpy(dtype=float) != rb.to_numpy(dtype=float)).sum()) if ra.shape == rb.shape else -1
    return dict(same_labels=shape, n_counts_differ=ns, n_responses_differ=nr, exact=bool(shape and ns == 0 and nr == 0))


def row_check(row: Mapping, ref: Mapping) -> dict:
    """The check run's fits.csv row against C2's: every shared column but META_SKIP, exactly (E110's `_num_same`); the
    fit_* columns and 'loss' apart (the row E110's check compares)."""
    cols = [c for c in row.keys() if c in ref.keys() and c not in META_SKIP]
    diff = {c: [row[c], ref[c]] for c in cols if not E110._num_same(row[c], ref[c])}
    fit_cols = [c for c in cols if c.startswith('fit_') or c == 'loss']
    return dict(n_columns=len(cols), n_fit_columns=len(fit_cols), differs=diff,
                fit_differs={c: v for c, v in diff.items() if c in fit_cols}, exact=bool(cols and not diff),
                fit_exact=bool(fit_cols and not any(c in diff for c in fit_cols)))


def evaluate_check(name: str, tag: str) -> dict:
    """Part 0's reading of one check run (the module docstring); 'passes'."""
    res = {}
    hp = fit_results_dir(name, check=True) / f'history_{tag}.csv'
    o = E83._load_pair(_check_local(name), tag)
    if not hp.exists() or o is None:
        return dict(passes=False, why='history or saved runs missing')
    h = pd.read_csv(hp)
    F = E83.saved_fit('real')
    res['n_evaluations'] = int(len(h))
    res['seconds_per_evaluation'] = float(h['time_s'].mean()) if 'time_s' in h else None
    res['apl_residual_history'] = sorted({float(x) for x in h[J.APL_RESIDUAL_KEY]}) if J.APL_RESIDUAL_KEY in h else None
    res['on_vs_c2_on'] = same_outputs(o, F['saved'], 'on', 'on')
    if name == 'f0_real':
        folder, ftag = E83.members()['real']
        hist = pd.read_csv(H.results_dir(folder) / f'history_{ftag}.csv')
        best = hist.loc[hist['loss'].idxmin()]
        ref_row = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[ftag]
        ex = E110.exact_check(h.iloc[0], best, ref_row)
        res['history_vs_c2_best'] = dict(reference=f'{folder}:{ftag} (history row {int(best["eval"])}, the best)', **ex)
        row = pd.read_csv(fit_results_dir(name, check=True) / 'fits.csv').set_index('tag').loc[tag]
        res['row_vs_c2_row'] = row_check(row, ref_row)
        res['off_vs_c2_off'] = same_outputs(o, F['saved'], 'off', 'off')
        res['passes'] = bool(ex['exact'] and res['row_vs_c2_row']['fit_exact'] and res['on_vs_c2_on']['exact']
                             and res['off_vs_c2_off']['exact'])
    else:
        res['off_vs_own_on'] = same_outputs(o, o, 'off', 'on')
        h0 = h.iloc[0]
        res['history_rates'] = dict(rate_on=float(h0['rate_on']), rate_off=float(h0['rate_off']),
                                    apl_raw_mean=float(h0['apl_raw_mean']), apl_cal_mean=float(h0['apl_cal_mean']))
        res['passes'] = bool(res['off_vs_own_on']['exact'] and res['on_vs_c2_on']['exact']
                             and float(h0['rate_on']) == float(h0['rate_off']))
    return res


def part_0(checks: Sequence[str] = CHECKS, rerun: bool = False) -> dict:
    """Part 0: runs each check of `checks` unless its history and launch record exist with the same command (`rerun`
    runs it again), then reads it; check.json ('f0_real', 'f1_real', 'all_pass', 'f0_pass', 'pending'). Raises if a check
    run fails."""
    specs = check_specs()
    res = dict(_read_json(out_dir() / 'check.json'))
    cc = None
    for name in checks:
        tag, cmd = specs[name]
        d = fit_results_dir(name, check=True)
        hp = d / f'history_{tag}.csv'
        rec = _read_json(d / 'launch.json')
        if rerun or not hp.exists() or rec.get('returncode') is None or not E110.same_command(rec.get('command') or [], cmd):
            cc = cc or E110._fit_clean()
            print(f'E111 0 {name}: ' + ' '.join(cmd), flush=True)
            rec = launch(cmd, d, local_dir('logs') / f'check_{name}.log', cc, name)
        r = dict(tag=tag, command=' '.join(rec['command']), launch={k: rec.get(k) for k in
                 ('host', 'started', 'finished', 'returncode', 'seconds', 'log', 'launcher', 'threads', 'pythonpath')},
                 commits=rec.get('commits'), queue_dirty_paths=rec.get('queue_dirty_paths'))
        if rec.get('returncode') != 0:
            r.update(passes=False, why=f"E10_joint_fit.py exited {rec.get('returncode')}")
        else:
            r.update(evaluate_check(name, tag))
        res[name] = r
    res['f0_pass'] = bool(res.get('f0_real', {}).get('passes'))
    res['pending'] = [n for n in CHECKS if n not in res]
    res['failed'] = [n for n in CHECKS if n in res and not res[n].get('passes')]
    res['all_pass'] = bool(not res['pending'] and not res['failed'])
    res['checked'] = _now()
    res['evaluated_by'] = dict(script='experiments/E111_residual_apl.py', host=_host(), **H.code_commits())
    _write_json(out_dir() / 'check.json', res)
    print(json.dumps({n: {k: res[n].get(k) for k in ('passes', 'n_evaluations', 'seconds_per_evaluation')} for n in CHECKS
                      if n in res}, indent=1, default=str), flush=True)
    if res['failed']:
        raise RuntimeError(f"E111 part 0 fails: {res['failed']} (check.json); nothing else runs until it passes")
    return res


def check_passed() -> bool:
    c = _read_json(out_dir() / 'check.json')
    return bool(c.get('f0_pass') and not c.get('failed'))


# ---------------------------------------------------------------------------
# E110's references and the single-condition runs (parts B and E)

def e110_b(w: str) -> Optional[dict]:
    """E110's counted part-B fit of a wiring with its saved runs: {'name', 'tag', 'row', 'p', 'o'}; None until both exist."""
    b = E110.counted_b(w)
    if b is None:
        return None
    o = E110.load_fit_outputs(b['name'], b['tag'])
    if o is None:
        return None
    return dict(b, p=params_of_row(b['row']), o=o)


def e110_verdict() -> dict:
    """E110's verdict.json's class: {'cls', 'final' (part C of E110 complete), 'source'}."""
    v = _read_json(H.results_dir('E110_no_kckc') / 'verdict.json')
    vv = v.get('verdict') or {}
    return dict(cls=vv.get('cls'), final=bool(vv) and not vv.get('c_pending', True), generated=v.get('generated'),
                source='results/E110_no_kckc/verdict.json')


def _save_pair(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, archive_dir() / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def _run_row(part: str, tag: str, w: str, f: float, p: Dict[str, float], o: Dict[str, tuple], cc: dict, lane: str,
             **extra) -> dict:
    off_p = J.apl_residual_params(p, f)
    row = dict(part=part, tag=tag, wiring=w, apl_residual=float(f), **{k: float(p[k]) for k in PARAM_KEYS},
               off_wAPLKC_scale=float(off_p['wAPLKC_scale']), off_wAPLPN_scale=float(off_p['wAPLPN_scale']),
               params_on=json.dumps({k: float(v) for k, v in p.items()}),
               params_off=json.dumps({k: float(v) for k, v in off_p.items()}), K_run=E83.readout_K(o), host=_host(),
               lane=lane, **extra)
    row.update(E83.realism(o))
    row.update(E110.run_extras(o))
    row.update(cc)
    return row


def _simulate(sim, p: Dict[str, float], f: float) -> tuple:
    return sim(J.apl_residual_params(p, f))


def part_b(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False, wait: bool = False
           ) -> pd.DataFrame:
    """Part B (the module docstring), resumable: a wiring is done when its dose rows are in any lane's runs_B table."""
    if not check_passed():
        raise RuntimeError('run --part 0 first (check.json: f0_real must pass)')
    lane = lane or _host()
    cc = E108._commits_clean()
    path, prov_path = out_dir() / f'runs_B_{lane}.csv', out_dir() / f'provenance_B_{lane}.json'
    while True:
        waiting = []
        for w in wirings:
            tags = [f'B{f:g}_{w}' for f in F_DOSE]
            if set(tags) <= set(runs_table_of('B')['tag']):
                print(f'E111 B {w}: done', flush=True)
                continue
            ref = e110_b(w)
            if ref is None:
                print(f"E111 B {w}: not yet, E110's counted part-B fit and its saved runs", flush=True)
                waiting.append(w)
                continue
            if not claim(f'B_{w}', lane, reclaim):
                print(f'E111 B {w}: claimed by another process, skipped', flush=True)
                continue
            try:
                E108._still_clean(cc)
                t0 = time.time()
                sim, info = E83.build_simulator(w, params=ref['p'])
                ident = None
                if w == 'real':
                    o0 = {'on': ref['o']['on'], 'off': _simulate(sim, ref['p'], 0.0)}
                    ident = same_outputs(o0, ref['o'], 'off', 'off')
                    if not ident['exact']:
                        raise RuntimeError(f"E111 B real: f = 0 does not reproduce E110's saved blocked run ({ident})")
                    print('E111 B real: f = 0 reproduces E110\'s blocked run exactly', flush=True)
                prov = _read_json(prov_path) or dict(part='B', lane=lane, script='experiments/E111_residual_apl.py', points={})
                for f, tag in zip(F_DOSE, tags):
                    o = {'on': ref['o']['on'], 'off': _simulate(sim, ref['p'], f)}
                    _save_pair('B', tag, o)
                    row = _run_row('B', tag, w, f, ref['p'], o, cc, lane, ref_name=ref['name'], ref_tag=ref['tag'],
                                   identity_f0=(None if ident is None else ident['exact']))
                    E83._append_row(path, row)
                    prov['points'][tag] = dict(H.model_provenance(orn_deltas=info['orn']), change=(
                        f"E110's counted part-B fit {ref['name']}:{ref['tag']} at fixed scales; the blocked condition with "
                        f"APL->KC and APL->bouton at {f:g} x the APL-on scales (E10_joint_fit.apl_residual_params); the "
                        "APL-on run is E110's saved one"), kws=info['kws'])
                    print(f"E111 B {w} f {f:g}: rate off {row['rate_off']:.4f} (x{row['ratio']:.2f}); flag "
                          f"'{row.get('realism_flag') or ''}'", flush=True)
                _write_json(prov_path, prov)
                print(f'E111 B {w}: {time.time() - t0:.0f}s', flush=True)
                del sim, info
                gc.collect()
            except BaseException:
                release(f'B_{w}')
                raise
        if not (wait and waiting):
            return runs_table_of('B')
        print(f'E111 B: waiting for {waiting} ({WAIT_SECONDS} s)', flush=True)
        time.sleep(WAIT_SECONDS)


def counted_c(w: str) -> Optional[dict]:
    """Part C's counted fit of a wiring: of C1_<w> and C2_<w> (where planned) the lower CHOOSE_COL ({'name', 'tag',
    'row', 'n_starts'}); None while any planned start is not done."""
    cands = []
    for s in fit_specs('C', [w]):
        r = fit_row(s.name, fit_tag(s))
        if r is None:
            return None
        cands.append(dict(name=s.name, tag=fit_tag(s), row=r))
    best = E110.choose_counted([c['row'] for c in cands])
    return dict(cands[best], n_starts=len(cands))


def e110_class() -> Optional[str]:
    """E110's class as far as this experiment needs it: its verdict once final, or its provisional one when that is
    "needed" or "not determined" (E110's provisional class changes only to "not determined", when its part C is flat,
    so either way part C runs on all eight and gives a class); None otherwise."""
    v = e110_verdict()
    return v['cls'] if (v['final'] or v['cls'] in E110_RUN_ALL) else None


def c_wirings() -> Tuple[str, ...]:
    """The wirings part C gives the verdict on: all eight unless E110's final verdict is "not needed" or "better
    without" (then real, P1, P2, M1)."""
    v = e110_verdict()
    return S4 if (v['final'] and v['cls'] in E110_NO_CLASS) else WIRINGS


def part_e(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False, wait: bool = False
           ) -> pd.DataFrame:
    """Part E (the module docstring), resumable: a wiring is done when its three rows are in any lane's runs_E table."""
    if not check_passed():
        raise RuntimeError('run --part 0 first (check.json: f0_real must pass)')
    lane = lane or _host()
    cc = E108._commits_clean()
    path, prov_path = out_dir() / f'runs_E_{lane}.csv', out_dir() / f'provenance_E_{lane}.json'
    while True:
        waiting = []
        for w in wirings:
            tags = [f'Efam_{w}', f'Eb_{w}', f'Ec_{w}']
            if set(tags) <= set(runs_table_of('E')['tag']):
                print(f'E111 E {w}: done', flush=True)
                continue
            if w not in c_wirings():
                print(f"E111 E {w}: no part C fit planned (E110's final verdict), skipped", flush=True)
                continue
            ref_b, c = e110_b(w), counted_c(w)
            oc = None if c is None else load_fit_outputs(c['name'], c['tag'])
            if ref_b is None or c is None or oc is None:
                print(f"E111 E {w}: not yet, E110's part-B fit and part C's counted fit with their saved runs", flush=True)
                waiting.append(w)
                continue
            if not claim(f'E_{w}', lane, reclaim):
                print(f'E111 E {w}: claimed by another process, skipped', flush=True)
                continue
            try:
                E108._still_clean(cc)
                t0 = time.time()
                F = E83.saved_fit(w)
                fc, pc = f_of_row(c['row']), params_of_row(c['row'])
                sim0, info0 = E83.build_simulator(w, params=ref_b['p'])
                if w == 'real':
                    ident = same_outputs({'off': _simulate(sim0, pc, fc)}, oc, 'off', 'off')
                    if not ident['exact']:
                        raise RuntimeError(f"E111 E real: part C's saved blocked run not reproduced at f_C ({ident})")
                    print("E111 E real: part C's blocked run reproduced exactly at f_C", flush=True)
                runs = {f'Eb_{w}': (ref_b['p'], fc, ref_b['o']['on'], f"E110's {ref_b['name']} blocked at f_C"),
                        f'Ec_{w}': (pc, 0.0, oc['on'], f"part C's {c['name']} with its residue removed (f = 0)")}
                outs = {t: {'on': on, 'off': _simulate(sim0, p, f)} for t, (p, f, on, _) in runs.items()}
                del sim0
                gc.collect()
                simF, infoF = E83.build_simulator(w)
                runs[f'Efam_{w}'] = (F['p'], fc, F['saved']['on'], f"the family's fit {F['folder']}:{F['tag']} blocked at f_C")
                outs[f'Efam_{w}'] = {'on': F['saved']['on'], 'off': _simulate(simF, F['p'], fc)}
                del simF
                gc.collect()
                prov = _read_json(prov_path) or dict(part='E', lane=lane, script='experiments/E111_residual_apl.py', points={})
                for tag in tags:
                    p, f, _, what = runs[tag]
                    _save_pair('E', tag, outs[tag])
                    row = _run_row('E', tag, w, f, p, outs[tag], cc, lane, f_c=fc, c_name=c['name'], c_tag=c['tag'],
                                   what=what)
                    E83._append_row(path, row)
                    prov['points'][tag] = dict(H.model_provenance(orn_deltas=infoF['orn']), change=what +
                                               ' (E10_joint_fit.apl_residual_params); the APL-on run is the saved one')
                    print(f"E111 E {tag}: f {f:g}, rate off {row['rate_off']:.4f}", flush=True)
                _write_json(prov_path, prov)
                print(f'E111 E {w}: {time.time() - t0:.0f}s', flush=True)
            except BaseException:
                release(f'E_{w}')
                raise
        if not (wait and waiting):
            return runs_table_of('E')
        print(f'E111 E: waiting for {waiting} ({WAIT_SECONDS} s)', flush=True)
        time.sleep(WAIT_SECONDS)


# ---------------------------------------------------------------------------
# the fits (parts C and D)

def fit_done(spec: FitSpec) -> bool:
    return fit_row(spec.name, fit_tag(spec)) is not None


def _prerequisite(spec: FitSpec) -> Optional[str]:
    """Why a fit cannot start yet (None when it can)."""
    if spec.start[0] == 'params' and spec.start[1] not in _read_json(starts_path()):
        return f"start {spec.start[1]} not in starts.json (E110's counted part-B fit first)"
    return None


def run_fits(specs: Sequence[FitSpec], lane: Optional[str] = None, reclaim: bool = False, dry_run: bool = False,
             wait: bool = False, gate=None) -> List[str]:
    """Runs the fits one process after another (`fit_command`), skipping those done or claimed (as E110's `run_fits`); a
    row per fit in runs_fits_<lane>.csv. `gate` (a callable returning (specs to add, why not yet or None)) adds fits that
    become due while the lane waits (part C on P3, P4, M2 and M3 on E110's final verdict)."""
    if not dry_run and not check_passed():
        raise RuntimeError('run --part 0 first (check.json: f0_real must pass)')
    cc = None if dry_run else E110._fit_clean()
    specs = list(specs)
    while True:
        gated = None
        if gate is not None:
            add, gated = gate()
            specs += [s for s in add if s.name not in {x.name for x in specs}]
        done, waiting = _run_fits_pass(specs, lane, reclaim, dry_run, cc)
        if gated:
            waiting.append(gated)
        if not (wait and waiting) or dry_run:
            if gated:
                print(f'E111: not run: {gated}', flush=True)
            return done
        print(f'E111: waiting for {waiting} ({WAIT_SECONDS} s)', flush=True)
        time.sleep(WAIT_SECONDS)


def _run_fits_pass(specs: Sequence[FitSpec], lane: Optional[str], reclaim: bool, dry_run: bool, cc: Optional[dict]
                   ) -> Tuple[List[str], List[str]]:
    lane = lane or _host()
    path = out_dir() / f'runs_fits_{lane}.csv'
    done, waiting = [], []
    for spec in specs:
        tag = fit_tag(spec)
        if fit_done(spec):
            print(f'E111 {spec.name}: done ({tag})', flush=True)
            done.append(spec.name)
            continue
        if not dry_run:
            write_starts()
        why = _prerequisite(spec)
        cmd = fit_command(spec)
        if dry_run:
            print(f'E111 {spec.name}{" (" + why + ")" if why else ""}: ' + ' '.join(cmd), flush=True)
            continue
        if why:
            print(f'E111 {spec.name}: not yet, {why}', flush=True)
            waiting.append(spec.name)
            continue
        if not claim(spec.name, lane, reclaim):
            print(f'E111 {spec.name}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            print(f'E111 {spec.name}: ' + ' '.join(cmd), flush=True)
            rec = launch(cmd, fit_results_dir(spec.name), local_dir('logs') / f'{spec.name}.log', cc, spec.name)
            if rec['returncode']:
                raise RuntimeError(f"E111 {spec.name}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
            r = fit_row(spec.name, tag)
            if r is None:
                raise RuntimeError(f"E111 {spec.name}: no row {tag} in its fits.csv (log {rec['log']})")
        except BaseException:
            release(spec.name)
            raise
        row = dict(name=spec.name, part=spec.part, wiring=spec.wiring, kckc_fixed=spec.kckc_fixed, f0=spec.f0,
                   maxfev=spec.maxfev, start=':'.join(spec.start), tag=tag, folder=f'{FIT_EXPERIMENT_PREFIX}/{spec.name}',
                   host=_host(), lane=lane, seconds=rec['seconds'], returncode=rec['returncode'], command=' '.join(cmd),
                   **params_of_row(r), apl_residual=f_of_row(r), loss=float(r['loss']), **{CHOOSE_COL: float(r[CHOOSE_COL])},
                   n_evals=int(r['n_evals']), **cc)
        E83._append_row(path, row)
        print(f"E111 {spec.name}: loss {float(r['loss']):.4f}, {CHOOSE_COL} {float(r[CHOOSE_COL]):.4f}, f "
              f"{f_of_row(r):.4g}, KC-KC {float(r['kc_kc_inh']):.3g}, {rec['seconds'] / 60:.0f} min", flush=True)
        done.append(spec.name)
    return done, waiting


def c_gate(wirings: Optional[Sequence[str]] = None, force: bool = False):
    """Part C's gate for P3, P4, M2 and M3: (their specs, None) once E110's verdict is final and "needed" or "not
    determined" (or `force`), ([], why) otherwise."""
    def gate():
        t4 = [s for s in fit_specs('C', wirings) if s.wiring in T4]
        if force or e110_class() in E110_RUN_ALL:
            return t4, None
        if not t4 or e110_verdict()['final']:
            return [], None
        return [], f"part C on {T4}: E110's verdict not in or not final ({e110_verdict()['cls']})"
    return gate


# ---------------------------------------------------------------------------
# readings

def block_pair_specific(d_on: pd.Series, d_off: pd.Series) -> dict:
    """The block's pair-specific remainder of one unit's distances, as analyses/yang_pair_deviations.py computes the
    flies' (in cosine similarity, 1 - distance): the isotonic residual of APL on on APL off, additive odor effects fitted
    on the 111 pairs outside the block, the mean over the nine block pairs of the residual beyond them ('block_ps'), and
    the block's own mean residual ('block_resid', minus the block on distances)."""
    from analyses import yang_pair_deviations as YPD
    on, off = E83._cs(d_on), E83._cs(d_off).reindex(E83._cs(d_on).index)
    idx = list(on.index)
    C, T = 1.0 - on.to_numpy(dtype=float), 1.0 - off.to_numpy(dtype=float)
    block = np.array([H.block_label(a, b) == H.FRUITY_SMALL_BLOCK for a, b in idx])
    res = YPD.iso_resid(C, T)
    X = YPD.odor_design(idx, H.APL16_ODORS)
    _, rem = YPD.odor_fit(res, X, ~block)
    return dict(block_ps=float(rem[block].mean()), block_resid=float(res[block].mean()))


def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (`E83.pair_distances`), cached in this process and in outputs/E111_residual_apl/
    distances_cache.pkl by the sha256 of the run's spike counts and the readout (as E110's)."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        try:
            _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist_store'] = {}
    store = _CACHE['dist_store']
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in store:
        store[key] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    return store[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            tmp = local_dir() / f'.distances_cache.{os.getpid()}.tmp'
            tmp.write_bytes(pickle.dumps(_CACHE['dist_store']))
            os.replace(tmp, local_dir() / 'distances_cache.pkl')
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


RUN_COLS = E110.RUN_COLS + ('apl_residual',)


def collect_runs() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run there is: ({(kind, wiring): outputs}, runs table). Kinds: 'family' (the saved fits), 'B110' (E110's
    counted part-B fits), 'B<f>' (part B's dose), 'C1', 'C2', 'C' (part C's counted), 'D', 'Efam', 'Eb', 'Ec'. The table
    has one row per run: kind, wiring, tag, folder, the absolute parameters, f, the fit's loss and CHOOSE_COL, realism,
    run_extras, the bands and 'flagged'."""
    outs: Dict[Tuple[str, str], Dict[str, tuple]] = {}
    rows = []

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        rows.append(dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E110.run_extras(o)))

    for w in WIRINGS:
        F = E83.saved_fit(w)
        fr = E110.family_fit_row(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], **F['p'], apl_residual=0.0, loss=float(fr['loss']),
            **{CHOOSE_COL: float(fr[CHOOSE_COL])}, n_evals=int(fr['n_evals']))
    for w in WIRINGS:
        b = e110_b(w)
        if b is not None:
            add('B110', w, b['o'], tag=b['tag'], folder=f"E110_no_kckc/fits/{b['name']}", name=b['name'], **b['p'],
                apl_residual=0.0, loss=float(b['row']['loss']), **{CHOOSE_COL: float(b['row'][CHOOSE_COL])},
                n_evals=int(b['row']['n_evals']))
    for part in ('B', 'E'):
        rt = runs_table_of(part)
        for _, r in rt.iterrows():
            o = load_run(part, r['tag'])
            if o is None:
                print(f"E111 score: {r['tag']} has a row but no outputs here or in the archive", flush=True)
                continue
            kind = (f"B{float(r['apl_residual']):g}" if part == 'B' else str(r['tag']).split('_', 1)[0])
            add(kind, r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/{part}', **{k: float(r[k]) for k in PARAM_KEYS},
                apl_residual=float(r['apl_residual']), **({'f_c': float(r['f_c'])} if 'f_c' in r and pd.notna(r['f_c']) else {}))
    for s in fit_specs('C') + fit_specs('D'):
        tag = fit_tag(s)
        r = fit_row(s.name, tag)
        if r is None:
            continue
        o = load_fit_outputs(s.name, tag)
        if o is None:
            print(f'E111 score: {s.name} has a fits.csv row but no outputs here or in the archive', flush=True)
            continue
        add(s.part, s.wiring, o, tag=tag, folder=f'{FIT_EXPERIMENT_PREFIX}/{s.name}', name=s.name, **params_of_row(r),
            apl_residual=f_of_row(r), loss=float(r['loss']), **{CHOOSE_COL: float(r[CHOOSE_COL])}, n_evals=int(r['n_evals']))
    for w in WIRINGS:
        c = counted_c(w)
        if c is None:
            continue
        src = c['name'].split('_', 1)[0]
        if (src, w) not in outs:
            continue
        meta = next(dict(r) for r in rows if r['kind'] == src and r['wiring'] == w)
        meta = {k: v for k, v in meta.items() if k not in ('kind', 'wiring')}
        meta.update(counted_from=src, n_starts=c['n_starts'])
        outs[('C', w)] = outs[(src, w)]
        rows.append(dict(kind='C', wiring=w, **meta))
    runs = pd.DataFrame(rows)
    if len(runs):
        bands = E110.bands_table(runs)
        runs = pd.concat([runs.reset_index(drop=True), bands.reset_index(drop=True)], axis=1)
        runs['flagged'] = [E108._flagged(x) for x in runs['realism_flag']]
    return outs, runs


def readings_table(outs: Mapping[Tuple[str, str], Dict[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run, and each population of UNITS whose wirings all have one) and readout:
    E110's `scores_of`, the block's pair-specific remainder (`block_pair_specific`), the runs' means of RUN_COLS
    ('run_<col>'), 'n_flagged', and the change from the family's same unit and readout ('d_<col>')."""
    rows = []
    kinds = list(dict.fromkeys(k for k, _ in outs))
    for kind in kinds:
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in UNITS.items() if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                single = len(ws) == 1 and unit in dd
                d = dd[unit] if single else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, unit_kind='wiring' if single else 'population', n_wirings=len(ws),
                                 readout=ro, n_flagged=int(s['flagged'].sum()),
                                 **{f'run_{c}': float(s[c].astype(float).mean()) for c in RUN_COLS if c in s.columns},
                                 **E110.scores_of(d), **block_pair_specific(d['on'], d['off'])))
    _save_dist_cache()
    rd = pd.DataFrame(rows)
    base = rd[rd['kind'] == 'family'].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c not in ('kind', 'unit', 'unit_kind', 'n_wirings', 'readout')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, r), np.nan) for u, r, v in zip(rd['unit'], rd['readout'], rd[c])]
    return rd


# ---------------------------------------------------------------------------
# part E's table, the flies' intervals and "told apart"

def data_step() -> Dict[str, pd.Series]:
    """The recordings' step from the six typical recordings with APL output blocked (TNT) to the two strongest (E76's
    `_steps`' step 2: the six's mean distance minus the two's), in Yang's z and in baseline-noise units, canonical pair
    order; and the two recordings' own agreement in z ('rho_own')."""
    if 'data_step' not in _CACHE:
        from experiments import E76_data_side_target as E76
        out = {}
        two = E76._actual_two('TNT')
        for u in ('z', 'noise'):
            R = E76.reading_rdms(u)
            s1, s2, own = E76._steps(R, two, 'TNT')
            out[u] = E83._cs(H.canonical_pairs(s2))
            if u == 'z':
                out['rho_own'] = float(H.spearman(own[0], own[1]))
        _CACHE['data_step'] = out
    return _CACHE['data_step']


def data_intervals() -> dict:
    """The flies' values and 90 % intervals of TOLD_MEASURES (the module docstring), and the ratio's beside:
    {key: {'value', 'lo', 'hi', 'source'}}."""
    from kcq import apl_effect as A
    out = {}
    lo, hi = H.data_target_ci('block')
    out['block'] = dict(value=H.data_target('block'), lo=lo, hi=hi, source="data-target registry 'block' (90 % interval)")
    s = pd.read_csv(H.shared_dir() / 'results' / 'yang_pair_deviations' / 'summary.csv').set_index('statistic')
    r = s.loc['block_pair_specific']
    out['block_ps'] = dict(value=float(r['value']), lo=float(r['lo']), hi=float(r['hi']),
                           source='results/yang_pair_deviations/summary.csv block_pair_specific (analyses/yang_pair_deviations.py)')
    rec = pd.read_csv(H.results_dir('E76_data_side_target') / 'recordings.csv')
    x = rec.loc[rec['genotype'] == 'TNT', 'rate_noise'].to_numpy(dtype=float)
    rng = np.random.default_rng(RATE_BOOT_SEED)
    boots = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(RATE_BOOT)])
    blo, bhi = np.percentile(boots, [5, 95])
    out['rate_off'] = dict(value=float(x.mean()), lo=float(blo), hi=float(bhi), n_recordings=int(len(x)),
                           registry_noise_tnt_frac=H.data_target('noise_tnt_frac'),
                           source=f'E76 recordings.csv rate_noise of the {len(x)} TNT recordings, mean, 90 % of {RATE_BOOT} '
                                  f'resamples (seed {RATE_BOOT_SEED})')
    band = float(A.opk_band('off', 'noise'))
    out['opk_off_tv'] = dict(value=None, lo=0.0, hi=band, source="the odors-per-KC band with APL off (kcq.apl_effect.opk_band"
                             "('off', 'noise'): the most atypical recording's distance from the flies' mean)")
    rlo, rhi = H.data_target_ci('noise_rate_ratio_fraction')
    out['ratio'] = dict(value=H.data_target('noise_rate_ratio_fraction'), lo=rlo, hi=rhi, reported_only=True,
                        source="data-target registry 'noise_rate_ratio_fraction' (90 % interval)")
    return out


def _step_rhos(d_f: Dict[str, pd.Series], d_0: Dict[str, pd.Series]) -> dict:
    """ρ between a model's step (blocked at f minus blocked at 0, pair by pair) and the recordings' step, in z and in
    baseline-noise units; the step's mean."""
    stp = E83._cs(d_f['off']) - E83._cs(d_0['off']).reindex(E83._cs(d_f['off']).index)
    ds = data_step()
    return dict(step_mean=float(stp.mean()), step_rho_z=float(H.spearman(stp, ds['z'])),
                step_rho_noise=float(H.spearman(stp, ds['noise'])))


# (family kind, its kind blocked at f_C, its kind blocked at 0): the step and the removal of each family
E_PAIRS = {'family': ('Efam', 'family'), 'B110': ('Eb', 'B110'), 'C': ('C', 'Ec')}


def e_table(outs: Mapping[Tuple[str, str], Dict[str, tuple]], runs: pd.DataFrame, rd: pd.DataFrame) -> pd.DataFrame:
    """Part E per family (FAMILY_KINDS), unit (each wiring with every run, and the population of the eight) and readout:
    the blocked condition's block and its pair-specific remainder (raw counts), response rate, ratio, response rate ρ and
    odors per KC with APL off (as fitted: f = 0 for the family and the family without KC-KC, f_C for the residue family);
    the step from f_C to 0 and its ρ against the recordings' step; the removal (blocked at 0 minus at f_C: response rate,
    mean pair distance)."""
    rows = []
    for fam, (k_f, k_0) in E_PAIRS.items():
        have = [w for w in WIRINGS if all((k, w) in outs for k in (fam, k_f, k_0))]
        units = [(w, (w,)) for w in have] + ([('all8', WIRINGS)] if len(have) == len(WIRINGS) else [])
        R = runs.set_index(['kind', 'wiring'])
        for ro in READOUTS:
            for unit, ws in units:
                if len(ws) == 1:
                    d_f, d_0 = _dist(outs[(k_f, unit)], ro), _dist(outs[(k_0, unit)], ro)
                else:
                    d_f = E83.population([_dist(outs[(k_f, w)], ro) for w in ws])
                    d_0 = E83.population([_dist(outs[(k_0, w)], ro) for w in ws])
                g = lambda kind, col: float(np.mean([float(R.loc[(kind, w), col]) for w in ws]))  # noqa: E731
                r_fit = rd[(rd['kind'] == fam) & (rd['unit'] == unit) & (rd['readout'] == 'raw')]
                row = dict(family=fam, label=FAMILY_KINDS[fam], unit=unit, readout=ro, n_wirings=len(ws),
                           f_c=float(np.mean([float(R.loc[('C', w), 'apl_residual']) for w in ws])) if ('C', ws[0]) in R.index else np.nan,
                           f_fit=g(fam, 'apl_residual'),
                           block=float(r_fit['block'].iloc[0]) if len(r_fit) else np.nan,
                           block_ps=float(r_fit['block_ps'].iloc[0]) if len(r_fit) else np.nan,
                           rate_off=g(fam, 'rate_off'), ratio=g(fam, 'ratio'), rate_rho_off=g(fam, 'rate_rho_off'),
                           opk_off_tv=g(fam, 'opk_off_tv'),
                           removal_rate_off=g(k_0, 'rate_off') - g(k_f, 'rate_off'),
                           removal_mean_dist_off=float(E83._cs(d_0['off']).mean() - E83._cs(d_f['off']).mean()),
                           **_step_rhos(d_f, d_0))
                rows.append(row)
    _save_dist_cache()
    t = pd.DataFrame(rows)
    if len(t):
        t['rho_own_data'] = data_step()['rho_own']
    return t


def told_apart(values: pd.DataFrame, intervals: Mapping[str, Mapping], a: str = 'C', b: str = 'family',
               measures: Sequence[ToldMeasure] = TOLD_MEASURES, n_need: Optional[int] = None) -> dict:
    """"Told apart by the recordings on hand" of families `a` and `b`: `values` has one row per family and wiring
    (columns 'family', 'unit' and the measures' columns); per measure and wiring, whether the flies' interval
    (`intervals[key]`: 'lo', 'hi') contains exactly one of the two families' values. Told apart when one measure does so on
    at least `n_need` wirings (default E110's `wirings_needed` of the wirings both have). Returns dict('told', 'cls', 'per'
    (DataFrame: measure, wiring, value_a, value_b, inside_a, inside_b, apart), 'counts' ({key: n}), 'measures_apart',
    'favoured' ({key: the family inside on most of its wirings apart}), 'n_needed', 'n_wirings')."""
    va = values[values['family'] == a].set_index('unit')
    vb = values[values['family'] == b].set_index('unit')
    wirings = [w for w in WIRINGS if w in va.index and w in vb.index]
    need = E110.wirings_needed(len(wirings)) if n_need is None else int(n_need)
    per, counts, favoured = [], {}, {}
    for m in measures:
        lo, hi = float(intervals[m.key]['lo']), float(intervals[m.key]['hi'])
        n, n_a = 0, 0
        for w in wirings:
            x, y = float(va.loc[w, m.col]), float(vb.loc[w, m.col])
            ia, ib = bool(np.isfinite(x) and lo <= x <= hi), bool(np.isfinite(y) and lo <= y <= hi)
            apart = ia != ib
            n += int(apart)
            n_a += int(apart and ia)
            per.append(dict(measure=m.key, wiring=w, value_a=x, value_b=y, lo=lo, hi=hi, inside_a=ia, inside_b=ib, apart=apart))
        counts[m.key] = n
        if n:
            favoured[m.key] = a if n_a * 2 >= n else b
    apart_keys = [k for k, n in counts.items() if n >= need and wirings]
    told = bool(apart_keys)
    return dict(told=told, cls=TOLD[0] if told else TOLD[1], per=pd.DataFrame(per), counts=counts,
                measures_apart=apart_keys, favoured={k: favoured[k] for k in apart_keys}, n_needed=need,
                n_wirings=len(wirings), a=a, b=b)


# ---------------------------------------------------------------------------
# the verdict

def f_at_bound(f: float, bounds: Tuple[float, float] = F_BOUNDS, tol: float = F_BOUND_TOL) -> bool:
    return bool(f <= bounds[0] + tol or f >= bounds[1] - tol)


def verdict(ch_family: pd.DataFrame, ch_b: pd.DataFrame, bands: pd.DataFrame, f: Mapping[str, float],
            e110_cls: Optional[str]) -> dict:
    """The class of part C (the module docstring) from E110's `changes_table` of the residue family against the family
    (`ch_family`) and against the family without KC-KC (`ch_b`), the residue family's bands (`bands`: index the wirings,
    columns E110's BANDS' keys) and its fitted f per wiring (`f`); `e110_cls` E110's verdict.

    Returns dict: 'cls' (one of CLASSES, or None when E110's verdict gives no class or is not in), 'reason', 'margins'
    (E110's `verdict` against the family: its not_needed_fails, groups, per_measure), 'gains_vs_b' (E110's gains against
    the family without KC-KC), 'f' (per wiring), 'n_f_below', 'n_f_at_bound', 'f_inside', 'n_needed', 'e110_cls'."""
    wirings = list(bands.index)
    n_need = E110.wirings_needed(len(wirings))
    v_fam = E110.verdict(ch_family, bands)
    v_b = E110.verdict(ch_b, bands)
    fv = {w: float(f[w]) for w in wirings}
    n_below = int(sum(x < F_NOTHING for x in fv.values()))
    n_bound = int(sum(f_at_bound(x) for x in fv.values()))
    f_inside = len(wirings) - n_bound >= n_need
    margins_hold = not v_fam['not_needed_fails']
    out = dict(margins=dict(not_needed_fails=v_fam['not_needed_fails'], groups=v_fam['groups'],
                            per_measure=v_fam['per_measure']),
               gains_vs_b=v_b['gains'], f=fv, n_f_below=n_below, n_f_at_bound=n_bound, f_inside=bool(f_inside),
               margins_hold=bool(margins_hold), n_needed=n_need, e110_cls=e110_cls, n_wirings=len(wirings))
    if e110_cls is None:
        out.update(cls=None, reason="E110's verdict is not in yet")
    elif e110_cls in E110_NO_CLASS:
        out.update(cls=None, reason=f'E110\'s verdict is "{e110_cls}": part C gives no class (the entry)')
    elif n_below >= n_need:
        out.update(cls=CLASSES[2], reason=f'the fitted f is below {F_NOTHING:g} on {n_below} of {len(wirings)} wirings')
    elif margins_hold and f_inside:
        out.update(cls=CLASSES[0], reason=("every condition of E110's \"not needed\" holds against the family and f is not "
                                           f"at a bound ({n_bound} of {len(wirings)} wirings at one)"))
    elif v_b['gains']:
        out.update(cls=CLASSES[1], reason=('a condition fails against the family ('
                                           + ('f at a bound' if margins_hold else ', '.join(v_fam['not_needed_fails']))
                                           + '), and against the family without KC-KC '
                                           + ', '.join(f"{g['label']} {g['pop']:+.4f} on {g['n_wirings']} wirings"
                                                       for g in v_b['gains'].values())))
    else:
        out.update(cls=CLASSES[2], reason=('a condition fails against the family, and no gain against the family without '
                                           'KC-KC on the population and on enough wirings'))
    return out


def d_reading(runs: pd.DataFrame) -> pd.DataFrame:
    """Part D per wiring: the fitted f and KC-KC beside the family's KC-KC; 'keeps_kckc' (KC-KC above D_KCKC_KEEP),
    'uses_residue' (f at F_NOTHING or above), choosing nothing."""
    rows = []
    for w in D_WIRINGS:
        d = runs[(runs['kind'] == 'D') & (runs['wiring'] == w)]
        fam = runs[(runs['kind'] == 'family') & (runs['wiring'] == w)]
        if not len(d) or not len(fam):
            continue
        r, rf = d.iloc[0], fam.iloc[0]
        rows.append(dict(wiring=w, f_d=float(r['apl_residual']), kckc_d=float(r['wKCKC_scale']),
                         kckc_family=float(rf['wKCKC_scale']), keeps_kckc=bool(float(r['wKCKC_scale']) > D_KCKC_KEEP),
                         uses_residue=bool(float(r['apl_residual']) >= F_NOTHING),
                         d_choose=float(r[CHOOSE_COL]) - float(rf[CHOOSE_COL])))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# part score

def _fmt(x, spec='.3f') -> str:
    return E110._fmt(x, spec)


COMPARISON_ROWS = E110.COMPARISON_ROWS + (("the block's pair-specific remainder (raw)", 'block_ps', 'raw'),
                                          ('residue f', 'run_apl_residual', 'raw'))


def comparison_table(rd: pd.DataFrame, v: Optional[dict], told: Optional[dict], et: pd.DataFrame, dr: pd.DataFrame,
                     intervals: dict) -> str:
    """comparison_table.md: the population of the eight per kind, the verdict, "told apart", part E and part D."""
    L = ['# E111: a residual APL output in the blocked flies in place of KC-KC inhibition', '',
         "Each column the mean matrix over the eight wirings (the family: E72 part B's fits; B110: E110's part-B refit, "
         "KC-KC 0; B<f>: B110 with the blocked condition at f, nothing re-held; C: the residue fitted, KC-KC 0, the counted "
         "start; Efam / Eb: the family / B110 blocked at the residue family's f; Ec: C with its residue removed). Run "
         'statistics are means over the eight runs.', '']
    order = ['family', 'B110'] + [f'B{f:g}' for f in F_DOSE] + ['C', 'C1', 'C2', 'Efam', 'Eb', 'Ec']
    kinds = [k for k in order if len(rd[(rd['kind'] == k) & (rd['unit'] == 'all8')])]
    if kinds:
        L += ['| statistic | ' + ' | '.join(kinds) + ' |', '|---|' + '---|' * len(kinds)]
        for lab, col, ro in COMPARISON_ROWS:
            sp = '.0f' if col in ('n_flagged', 'run_n_kc_ge20_off', 'run_max_count_off') else ('.4g' if 'scale' in col else '.3f')
            L.append(f'| {lab} | ' + ' | '.join(_fmt(E110._v(rd, k, 'all8', col, ro), sp) for k in kinds) + ' |')
    if v is not None:
        L += ['', f"## The verdict (part C against the family): {v['cls']}", '', v['reason'] + '.', '',
              f"f per wiring: {', '.join(f'{w} {x:.4g}' for w, x in v['f'].items())}; below {F_NOTHING:g} on "
              f"{v['n_f_below']}, at a bound on {v['n_f_at_bound']}.", '']
        for g, items in v['margins']['groups'].items():
            L += [f'- against the family, {g}: ' + '; '.join(items)]
        if v['gains_vs_b']:
            L += ['- gains against the family without KC-KC: ' + '; '.join(
                f"{x['label']} {x['pop']:+.4f} ({x['n_wirings']} wirings)" for x in v['gains_vs_b'].values())]
    if told is not None:
        L += ['', f"## Told apart (the residue family against the family): {told['cls']}", '',
              '| measure | flies 90 % | wirings apart | favoured |', '|---|---|---|---|']
        for m in TOLD_MEASURES:
            iv = intervals[m.key]
            L.append(f"| {m.label} | {_fmt(iv['lo'], '.4f')} to {_fmt(iv['hi'], '.4f')} | {told['counts'].get(m.key, 0)} of "
                     f"{told['n_wirings']} | {told['favoured'].get(m.key, '-')} |")
    if len(et):
        L += ['', '## Part E (the population of the eight, raw counts)', '',
              '| family | f_C | block | pair-specific | rate off | ratio | rate ρ off | odors per KC off | step ρ (z) | '
              'step ρ (noise) | removal: rate off | removal: mean distance off |', '|---|' + '---|' * 11]
        for r in et[(et['unit'] == 'all8') & (et['readout'] == 'raw')].itertuples():
            L.append(f'| {r.label} | {_fmt(r.f_c, ".4f")} | {_fmt(r.block)} | {_fmt(r.block_ps, ".4f")} | {_fmt(r.rate_off, ".4f")} '
                     f'| {_fmt(r.ratio, ".2f")} | {_fmt(r.rate_rho_off)} | {_fmt(r.opk_off_tv)} | {_fmt(r.step_rho_z)} | '
                     f'{_fmt(r.step_rho_noise)} | {_fmt(r.removal_rate_off, "+.4f")} | {_fmt(r.removal_mean_dist_off, "+.4f")} |')
    if len(dr):
        L += ['', f"## Part D (both free): keeps KC-KC on {int(dr['keeps_kckc'].sum())}, uses a residue on "
              f"{int(dr['uses_residue'].sum())} of {len(dr)}", '', '| wiring | f | KC-KC | KC-KC (family) |', '|---|---|---|---|']
        for r in dr.itertuples():
            L.append(f'| {r.wiring} | {r.f_d:.4g} | {r.kckc_d:.3g} | {r.kckc_family:.3g} |')
    return '\n'.join(L) + '\n'


def b_dose(rd: pd.DataFrame) -> pd.DataFrame:
    """Part B's dose on the population of the eight and per wiring: per f the blocked condition's readings beside B110's."""
    keep = ('kind', 'unit', 'readout', 'run_apl_residual', 'run_rate_off', 'run_ratio', 'run_rate_rho_off', 'run_opk_off_tv',
            'mean_effect', 'all_loss', 'all_pattern', 'all_amplitude', 'all_shape', 'all_sd_ratio', 'rho_on', 'rho_off',
            'block', 'block_ps', 'mean_off', 'n_flagged')
    sub = rd[rd['kind'].isin(['B110'] + [f'B{f:g}' for f in F_DOSE])]
    return sub[[c for c in keep if c in sub.columns]]


def score() -> dict:
    """Part score (no simulations; the module docstring)."""
    out = out_dir()
    outs, runs = collect_runs()
    runs.to_csv(out / 'runs.csv', index=False)
    rd = readings_table(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    b_dose(rd).to_csv(out / 'b_dose.csv', index=False)
    iv = data_intervals()
    _write_json(out / 'data_intervals.json', iv)
    E108.e106_yardsticks().to_csv(out / 'yardsticks_e106.csv', index=False)
    et = e_table(outs, runs, rd)
    et.to_csv(out / 'e_table.csv', index=False)
    dr = d_reading(runs)
    dr.to_csv(out / 'd_reading.csv', index=False)
    res = dict(e110=e110_verdict())
    if len(dr):
        res['d'] = dict(n_keeps_kckc=int(dr['keeps_kckc'].sum()), n_uses_residue=int(dr['uses_residue'].sum()),
                        n_wirings=int(len(dr)))
        if set(dr['wiring']) == set(D_WIRINGS):
            E110.changes_table(rd, 'D', 'family', D_WIRINGS, 'S4').to_csv(out / 'changes_D_family.csv')
    ws = c_wirings()
    have_c = all(('C', w) in outs for w in ws) and all(('B110', w) in outs for w in ws)
    v = told = None
    if have_c:
        pop = 'all8' if len(ws) == len(WIRINGS) else 'S4'
        chf = E110.changes_table(rd, 'C', 'family', ws, pop)
        chb = E110.changes_table(rd, 'C', 'B110', ws, pop)
        chf.to_csv(out / 'changes_C_family.csv')
        chb.to_csv(out / 'changes_C_B.csv')
        bands = runs[runs['kind'] == 'C'].set_index('wiring').loc[list(ws), [b.key for b in E110.BANDS]]
        bands.to_csv(out / 'bands_C.csv')
        f = runs[runs['kind'] == 'C'].set_index('wiring')['apl_residual'].to_dict()
        v = verdict(chf, chb, bands, f, e110_class())
        pt = E110.pending_table(rd, 'C', 'family', ws, pop)
        pt.to_csv(out / 'pending_C.csv')
        v['pending'] = {k: dict(label=r['label'], pop=float(r['pop']), mean_paired=float(r['mean_paired']),
                                fails=bool(r['fails'])) for k, r in pt.iterrows()}
        res['verdict'] = v
    if len(et):
        ew = et[(et['readout'] == 'raw') & (et['unit'] != 'all8')]
        if set(ew.loc[ew['family'] == 'C', 'unit']) and set(ew.loc[ew['family'] == 'family', 'unit']):
            told = told_apart(ew, iv, 'C', 'family')
            told_b = told_apart(ew, iv, 'C', 'B110')
            pd.concat({'C_vs_family': told['per'], 'C_vs_B110': told_b['per']}).to_csv(out / 'told_apart.csv')
            res['told_apart'] = {k: v_ for k, v_ in told.items() if k != 'per'}
            res['told_apart_vs_b'] = {k: v_ for k, v_ in told_b.items() if k != 'per'}
    (out / 'comparison_table.md').write_text(comparison_table(rd, v, told, et, dr, iv))
    _write_json(out / 'verdict.json', dict(res, generated=_now(), **H.code_commits()))
    _write_json(out / 'provenance_score.json', dict(part='score', script='experiments/E111_residual_apl.py', generated=_now(),
                                                    host=_host(), **H.code_commits()))
    print(json.dumps(dict(verdict=None if v is None else {k: v[k] for k in ('cls', 'reason', 'n_f_below', 'n_f_at_bound')},
                          told=None if told is None else told['cls'], d=res.get('d'), e110=res['e110']), indent=1,
                     default=str), flush=True)
    return res


# ---------------------------------------------------------------------------
# the archive's tables

KIND_TEXT = {'B': "E110's part-B refit (KC-KC 0) at fixed scales, the blocked condition at f x the APL-on output scales",
             'C1': 'KC-KC held at 0; threshold, APL->KC, APL->bouton and the residue f fitted from E110\'s part-B optimum, f from 0.05',
             'C2': 'KC-KC held at 0; threshold, APL->KC, APL->bouton and the residue f fitted from E110\'s part-B optimum, f from 0.2',
             'D': "KC-KC and the residue f both fitted from the family's fit, f from 0.02",
             'Efam': "the family's fit with the blocked condition at the residue family's f",
             'Eb': "E110's part-B refit with the blocked condition at the residue family's f",
             'Ec': 'the residue family with its residue removed (the blocked condition at f = 0)'}


def write_variants() -> pd.DataFrame:
    """variants.csv (one row per run of this experiment, absolute parameters and f) and model_pair_distances_long.csv
    (raw spike counts, the harness's primary model RDM) for the model-output archive and kcq/apl_effect.py."""
    out = out_dir()
    outs, runs = collect_runs()
    rows, parts = [], []
    for r in runs[~runs['kind'].isin(['family', 'B110', 'C'])].itertuples():
        o = outs[(r.kind, r.wiring)]
        base_kind = 'B' if str(r.kind).startswith('B') else r.kind
        subdir = ('B' if base_kind == 'B' else 'E' if base_kind in ('Efam', 'Eb', 'Ec') else r.folder.split('/', 1)[1])
        f = float(r.apl_residual)
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.kind, base='bouton', k=K, wiring=r.wiring,
                         kind=f"C2's family's {r.wiring}: {KIND_TEXT[base_kind]}", ref_experiment=E83.members()[r.wiring][0],
                         ref_fit=E83.members()[r.wiring][1], fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale, apl_residual=f,
                         off_wAPLKC_scale=f * float(r.wAPLKC_scale), off_wAPLPN_scale=f * float(r.wAPLPN_scale),
                         rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on, max_count_off=r.max_count_off,
                         realism_flag=r.realism_flag, subdir=subdir))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# commands

def lane_line(part: str, lane: str, threads: int, wirings: Optional[Sequence[str]] = None, wait: bool = True) -> str:
    w = '' if wirings is None else ' --wirings ' + ' '.join(wirings)
    return (f'OMP_NUM_THREADS={threads} OMP_THREAD_LIMIT={threads} nice -n 10 ./venv/bin/python -u '
            f'experiments/E111_residual_apl.py --part {part}{w} --lane {lane}' + (' --wait' if wait else ''))


def print_commands() -> None:
    """Every fit's E10_joint_fit.py command and example lane lines."""
    for part in ('C', 'D'):
        for s in fit_specs(part):
            print(f'# {s.name} -> {fit_tag(s)}')
            print(' '.join(['./venv/bin/python' if i == 0 else x for i, x in enumerate(fit_command(s))]))
    print('# lanes (each lane runs the parts in this order; the claims split the points):')
    for lane, t in (('atlas', 8), ('hal1', 16), ('hal2', 16)):
        for part in ('D', 'B', 'C', 'E'):
            print(lane_line(part, lane, t, wait=part != 'D'))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'B', 'C', 'D', 'E', 'score', 'variants', 'commands', 'starts'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS), help="a subset of the part's wirings")
    ap.add_argument('--lane', default=None, help="the lane's name (its tables runs_<B|E|fits>_<lane>.csv; default the host)")
    ap.add_argument('--reclaim', action='store_true', help='take over claims (a process that died left them)')
    ap.add_argument('--dry-run', action='store_true', help="print the fits' commands, run nothing")
    ap.add_argument('--force', action='store_true', help="part C on P3, P4, M2 and M3 whatever E110's verdict")
    ap.add_argument('--rerun', action='store_true', help='part 0: run the checks again')
    ap.add_argument('--checks', nargs='+', default=list(CHECKS), choices=list(CHECKS), help='part 0: the checks to run')
    ap.add_argument('--wait', action='store_true',
                    help=f'retry points whose prerequisites (E110\'s or part C\'s fits) are not saved yet every {WAIT_SECONDS} s')
    a = ap.parse_args()
    if a.part == '0':
        part_0(a.checks, a.rerun)
    elif a.part == 'B':
        part_b(a.wirings or WIRINGS, a.lane, a.reclaim, a.wait)
    elif a.part == 'C':
        run_fits(fit_specs('C', a.wirings, with_t4=False), a.lane, a.reclaim, a.dry_run, a.wait, gate=c_gate(a.wirings, a.force))
    elif a.part == 'D':
        run_fits(fit_specs('D', a.wirings), a.lane, a.reclaim, a.dry_run, a.wait)
    elif a.part == 'E':
        part_e(a.wirings or WIRINGS, a.lane, a.reclaim, a.wait)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'commands':
        print_commands()
    elif a.part == 'starts':
        print(json.dumps(write_starts(), indent=1))


if __name__ == '__main__':
    main()
