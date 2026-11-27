#!/usr/bin/env python3
"""E147's runner: the reads of steps 1 and 2 of the megamat stage error read on model flies with each ORN fly's own
input (experiment_priorities.md, E147; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-27
09:01; hal's part, beside atlas's experiments/E147_megamat_input_flies.py, which gives the input flies, the points and
tags, the folders, the lock and the score, and is used here as it is).

A read is one model fly: one wiring with one input fly's input (atlas's `input_fly`), at a setting's parameters, both
response rates re-held at the panel's operating point as E136's `hold_at` holds a read (`hold_with_input`: its steps,
with the input and a start in place of the panel's input and the held megamat point; E77's `hold_rates` on E83's
`build_simulator` with the input: the somatic threshold with APL off, then the APL->KC scale with its tied scales; E134's
`held_targets` on megamat, the wiring's own fit's counts on APL16), the breadth the lab input's, then built again at the held parameters and run with every stage saved
(E134's `build_family` and `run_stages`), which must give the hold's outputs exactly. The wiring keeps its seed for
pooling and noise, since the run is saved under its wiring and read as a wiring is read. The hold starts where the
setting's reading at the lab's input started, so that an input fly with no deviation reproduces that reading exactly:
- 'family' (the unrefitted family): from the wiring's fit (E83's `saved_fit`), as E134 held its megamat runs and as the
  fits on APL16 were made;
- 'A1' (E136's A1 at its classed point, the lock's A1_classed_z, which must be E136's `fitted_z('A1')`): from the
  wiring's held megamat point with A1's factors (E136's `start_params`), as E136's A1 reads were held.
The row's thr_factor and apl_factor are the hold's, against that start (its 'start' column); thr_vs_own and apl_vs_own
are against the setting's reading at the lab's input on the panel (E134's held megamat run, the wiring's fit on APL16,
E136's A1 read).

The free reading ('familyfree', 'A1free', beside, classing nothing): each model fly at the parameters of its setting's
reading at the lab's input on the panel, nothing re-held, one run of each condition.

"Held" in the rows is the lock's (revised 09:45; atlas's `is_held`, which the score applies itself): both counts within
HOLD_TOL (E136's, 0.5 %) of their targets and neither factor at a bound. Beside it, 'held_exact': both counts at their
targets, which E77's bisection misses by a count or two when its 20 steps end between two counts (E134's held runs on
real, M1 and P4; 12 of E85's 72 held model flies). Unheld model flies are saved and counted, never dropped; the free
reading is not held.

Parts (CLI `--part`):
- identity: before any read, the input checks (`input_checks`: the zero-deviation input fly is the lab's input exactly,
  every input fly has the lab input's glomeruli and odors and no NaN, the mean over the nine is the lab's input, and on
  APL16 paa, euc, -aPine and B-myr are the lab's in every fly), then the zero-deviation input fly read as a read on
  `--wirings` (default IDENTITY_WIRINGS) for each of `--settings` and `--panels` (default the family on megamat): its
  counts, responses and every stage matrix reproduce the setting's reading at the lab's input exactly (E134's held
  megamat run; the wiring's fit and E136's A1 reads where asked). Records in results/.../identity/identity_<lane>.json;
  RuntimeError on any difference. The reads refuse to start until the family's megamat identity is exact on one wiring
  and no identity has failed (`require_identity`).
- reads: one lane of a step's reads (`--step 1|2`, atlas's `planned_points`: step 1 the family on megamat, the nine
  input flies on the eight wirings, the T wirings first, 72; step 2 on T, A1 on megamat, then A1 and the family on
  APL16, 36 each), `--free` for the same points under the free reading.

Lanes and claims (E136's reads): a lane is '<host>' or '<host>-<n>' (`lane_name`), so that two lanes never write one
table; each read is claimed by mkdir of results/.../claims/<tag>/ before it starts and skipped when claimed elsewhere or
done (its row in any lane's reads_*.csv, atlas's `reads_table`); `--reverse` takes the list from its other end, so that
two lanes meet in the middle; `--reclaim` retakes this lane's own stale claims only (a lane stopped by a crash), never
another lane's. A failed read releases its claim. One process runs a lane's reads one after another; each read's
simulators are dropped and collected before the next is built (pyolfsysm's fixed-parameter simulator is a reference
cycle; E136's `_rebuild`, queue c51ee3e).

Outputs:
- outputs/E147_megamat_input_flies/step<k>/<setting>/<wiring>/<fly>/ (atlas's `run_dir`; local, large): each read in
  E134's step-3 layout (every stage, the units' metadata, the input), also in the model-output archive's twin
  (atlas's `run_dirs`; `--no-archive` keeps it local); the identity's runs under outputs/.../identity/.
- <shared>/results/E147_megamat_input_flies/ (compact): reads_<lane>.csv (atlas's READ_COLUMNS first, then this
  runner's columns), provenance_reads_<lane>.json, claims/; identity/ (identity_<lane>.json, its rows, claims and
  provenance).

Holds (CLAUDE.md, the entry): no olfsysm change (the wirings keep 'first' claws, pinned, as the family's fits have
them); nothing searched, fitted or weighted (the inputs are atlas's, the recorded flies' own deviations); validation2 is
not read; Pratyush's hemibrain only.

Usage (from the repository root, with -m, so that a worktree's own experiments/ package is the one imported and not the
venv's editable install; a megamat read takes about 280 s at 16 threads on viki or disco, several times that on hal; a
read is mostly single-threaded, so a host may run two lanes):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E147_reads --part identity --lane <host>
    ... --part reads --step 1 --lane <host>[-<n>] [--reverse]
    ... --part reads --step 2 --lane <host>[-<n>] [--reverse]
    ... --part reads --step 1 --free --lane <host>[-<n>]        # the free reading, beside
"""
from __future__ import annotations

import argparse
import contextlib
import gc
import json
import os
import re
import socket
import sys
import time
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E147_megamat_input_flies'
SCRIPT = 'experiments/E147_reads.py'
# the measure under which E136's A1 is fitted and read (its measure_strength/ folder)
E136_MEASURE = 'strength'
# where each re-held setting's hold starts (the module docstring)
STARTS: Dict[str, str] = {'family': 'fit', 'A1': 'held_megamat_point'}
# the lock's HOLD_TOL (E136's): a read within this relative miss of both targets, no factor at a bound, is held
HOLD_TOL = 0.005
IDENTITY_WIRINGS: Tuple[str, ...] = ('P3',)
IDENTITY_FLY = 'zero'
IDENTITY_DIR = 'identity'
# the mean over the input flies against the lab's input (atlas: "exactly", up to float rounding of the sums)
INPUT_MEAN_TOL = 1e-9
# a claim folder without its claim.json this old (s) is an orphan that `--reclaim` may retake
ORPHAN_CLAIM_S = 600.0


# ---------------------------------------------------------------------------
# atlas's module, E136 on its measure, the lane

def atlas():
    """atlas's experiments/E147_megamat_input_flies.py: the input flies (`input_fly`, `input_flies`, `lab_input`,
    `zero_deviation_input`), the points and tags (`planned_points`, `point_tag`), the folders (`results_dir`,
    `local_dir`, `run_dirs`), the rows (`reads_table`, READ_COLUMNS, PARAM_KEYS) and the lock (`require_lock`)."""
    from experiments import E147_megamat_input_flies as M
    return M


def input_of(fly: str, panel: str) -> pd.DataFrame:
    """Input fly `fly`'s model input on `panel` (atlas's `input_fly`)."""
    return atlas().input_fly(fly, panel)


def lab_of(panel: str) -> pd.DataFrame:
    """The zero-deviation input fly on `panel` (atlas's `zero_deviation_input`, the lab's input)."""
    return atlas().zero_deviation_input(panel)


@contextlib.contextmanager
def e136() -> Iterator[object]:
    """E136's module with its measure set to E136_MEASURE for the block (restored after), so that its A1 fit record
    and reads are those of measure_strength/."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        yield F
    finally:
        F.set_measure(old)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


_COMMITS: Dict[str, object] = {}


def commits() -> Dict[str, object]:
    """The checkouts' commits and dirty flags (`H.code_commits()`, scalars only), once per process: the code a lane
    runs is the code it loaded."""
    if not _COMMITS:
        _COMMITS.update({k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    return dict(_COMMITS)


def lane_name(lane: Optional[str] = None) -> str:
    """The lane: '<host>' (the default) or '<host>-<n>' of this host. ValueError otherwise, since a lane's table and
    provenance are written by its own process only and a lane name of another host's would share them."""
    h = _host()
    lane = h if lane is None else str(lane)
    if not re.fullmatch(rf'{re.escape(h)}(-\d+)?', lane):
        raise ValueError(f"a lane is '<host>' or '<host>-<n>' on this host ({h!r}), not {lane!r}")
    return lane


def claim(tag: str, lane: str, root: Path, reclaim: bool = False) -> bool:
    """A read claimed by mkdir of <root>/claims/<tag>/ (E134's `claim`, atomic on NFS); False when another process
    holds it. `reclaim` retakes a claim of this lane only (its claim.json names this lane and host), or one without a
    claim.json older than ORPHAN_CLAIM_S (its process stopped between the mkdir and the write); never another lane's."""
    from experiments import E134_population_baseline as B
    if B.claim(tag, lane, False, root=root):
        return True
    if not reclaim:
        return False
    d = Path(root) / 'claims' / tag
    held = _read_json(d / 'claim.json')
    if not held:
        try:
            if time.time() - d.stat().st_mtime < ORPHAN_CLAIM_S:
                return False
        except FileNotFoundError:
            return B.claim(tag, lane, False, root=root)
    elif held.get('lane') != lane or held.get('host') != _host():
        return False
    return B.claim(tag, lane, True, root=root)


# ---------------------------------------------------------------------------
# the settings

def split_setting(setting: str) -> Tuple[str, bool]:
    """(the re-held setting, whether the reading is free) of a setting of atlas's SETTINGS."""
    if setting not in atlas().SETTINGS:
        raise ValueError(f'{setting=} not in {atlas().SETTINGS}')
    return (setting[:-len('free')], True) if setting.endswith('free') else (setting, False)


def a1_z(lock: Mapping) -> Dict[str, float]:
    """A1's classed point: the lock's settings.A1_classed_z, which must be E136's `fitted_z('A1')` (measure_strength);
    RuntimeError otherwise, or when E135 hands A1 a table of PN rates (E136's `a1_stage` caches that table per wiring
    and panel, not per input, so another input would read the lab input's table)."""
    z = {k: float(v) for k, v in lock['settings']['A1_classed_z'].items()}
    with e136() as F:
        zf = F.fitted_z('A1')
        mode = F.require_e135()['mode']
    if set(z) != set(zf) or any(z[k] != zf[k] for k in z):
        raise RuntimeError(f"the lock's A1_classed_z {z} is not E136's fitted_z('A1') {zf}")
    if mode != 'search':
        raise RuntimeError(f"A1's PN stage is E135's table (mode {mode!r}), which E136 does not build per input")
    return z


def setting_model(base: str, wiring: str, lock: Mapping) -> dict:
    """A re-held setting on a wiring: dict(arm: E136's arm name, z: its coordinates, p0: the parameters its hold starts
    from (None: E136's `start_params`), start: where that is)."""
    from experiments import E83_cross_wiring_fit as E83
    if base == 'family':
        with e136() as F:
            arm = F.FAMILY
        return dict(arm=arm, z={}, p0={k: float(v) for k, v in E83.saved_fit(wiring)['p'].items()},
                    start=STARTS['family'])
    if base == 'A1':
        return dict(arm='A1', z=a1_z(lock), p0=None, start=STARTS['A1'])
    raise ValueError(f'{base=}')


def lab_reading(base: str, panel: str, wiring: str) -> dict:
    """The setting's reading at the lab's input on a panel and wiring, which an input fly with no deviation reproduces:
    the family on megamat E134's held run (megamat_held_<wiring>), on APL16 the wiring's fit (E83's `saved_fit`, its
    APL-off condition the APL->KC and APL->bouton scales at 0); A1 E136's read of it (measure_strength,
    A1_<panel>_<wiring>). Returns dict(params, params_off, outputs: {'on', 'off'}: (responses, spike counts), run:
    the SavedRun with its stage matrices or None, source). FileNotFoundError where it is not saved."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E134_population_baseline as B
    if base == 'family' and panel == 'megamat':
        tag = B.tag_of(B.Point('megamat', 'held', wiring))
        run = B.load_run(tag)
        if run is None:
            raise FileNotFoundError(f"E134's {tag} is not saved")
        source = f'E134 {tag}'
    elif base == 'family' and panel == 'apl16':
        S = E83.saved_fit(wiring)
        p = {k: float(v) for k, v in S['p'].items()}
        return dict(params=p, params_off=dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0), outputs=S['saved'], run=None,
                    source=f"the wiring's fit (E83 saved_fit, {S['folder']}:{S['tag']})")
    elif base == 'A1':
        with e136() as F:
            pt = F.ReadPoint('A1', panel, wiring, 'held' if panel == 'megamat' else None)
            run = F.load_read(pt)
            tag = F.read_tag(pt)
        if run is None:
            raise FileNotFoundError(f"E136's {tag} (measure_strength) is not saved")
        source = f'E136 {tag}'
    else:
        raise ValueError(f'{base=} {panel=}')
    return dict(params={k: float(v) for k, v in run.info['params'].items()},
                params_off={k: float(v) for k, v in run.info['params_off'].items()},
                outputs={c: (run.conditions[c]['responses'], run.conditions[c]['spike_counts']) for c in ('on', 'off')},
                run=run, source=source)


_OWN: Dict[Tuple[str, str, str], Optional[dict]] = {}


def own_params(base: str, panel: str, wiring: str) -> Optional[dict]:
    """`lab_reading`'s params, params_off and source (cached); None where it is not saved."""
    key = (base, panel, wiring)
    if key not in _OWN:
        try:
            r = lab_reading(base, panel, wiring)
            _OWN[key] = dict(params=r['params'], params_off=r['params_off'], source=r['source'])
        except FileNotFoundError:
            _OWN[key] = None
    return _OWN[key]


# ---------------------------------------------------------------------------
# one read

def n_responses(responses: pd.DataFrame) -> int:
    """The number of responses (two or more spikes) of a responses table."""
    return int((responses.to_numpy(dtype=float) > 0.5).sum())


def held_flags(n_on: int, n_off: int, target_on: int, target_off: int, at_bound: bool,
               tol: float = HOLD_TOL) -> Tuple[bool, bool]:
    """(held_exact, held): both counts at their targets and no factor at a bound, and both within `tol` of their
    targets (relative) and no factor at a bound (the lock's, atlas's `is_held`)."""
    exact = n_on == target_on and n_off == target_off and not at_bound
    near = abs(n_on - target_on) <= tol * target_on and abs(n_off - target_off) <= tol * target_off and not at_bound
    return bool(exact), bool(near)


def hold_with_input(arm: str, z: Mapping[str, float], wiring: str, panel: str, orn: pd.DataFrame,
                    p0: Optional[Mapping[str, float]], lock134: Mapping, label: str = 'E147') -> dict:
    """E136's `hold_at` with an input and a start of the caller's: E136's arm `arm` on `wiring` at coordinates `z`
    applied as log factors (E136's FACTORS) to `p0` (None: the wiring's held megamat point, E136's `start_params`),
    both response rates held at the panel's operating point (E136's `read_targets`) by E77's `hold_rates` on E83's
    `build_simulator` with `orn` and the lab input's breadth (`wiring_parts(wiring)['breadth']`), the arm's kwargs
    (E136's `arm_kws`, A1's `a1_kws`), APL bracket (E136's HOLD_RANGE, or HOLD_RANGE_DIVISIVE_APL) and APL-off
    parameters (`off_extra_of`). Line for line `hold_at`'s steps, which read the panel's input and the held megamat
    point only; E136 is left as it is (its lanes run). Returns `hold_at`'s dict."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import hold_rates
    with e136() as F:
        base = F.base_arm(arm)
        a = F.ARMS[base]
        p0 = F.start_params(wiring)[0] if p0 is None else {k: float(v) for k, v in p0.items()}
        if 'kckc_fixed' in a.kws:
            p0['wKCKC_scale'] = float(a.kws['kckc_fixed'])
        for k, v in z.items():
            if k in F.FACTORS:
                for kk in F.FACTORS[k]:
                    p0[kk] = p0[kk] * float(np.exp(v))
        ct = F.CLAW_THR_UNIT * float(z['claw_thr']) if 'claw_thr' in z else None
        kws = dict(F.arm_kws(base), **({'claw_thr': ct} if ct else {}))
        if base == 'A1':
            kws.update(F.a1_kws(wiring, panel, orn, z))
        Pp = E83.wiring_parts(wiring)
        targets = F.read_targets(wiring, panel, orn.shape[1], lock134)
        rng = F.HOLD_RANGE_DIVISIVE_APL if kws.get('apl_inh_divisive') else F.HOLD_RANGE
        off_extra = F.off_extra_of(base, p0)
    hold = hold_rates(lambda: E83.build_simulator(wiring, orn=orn, breadth=Pp['breadth'], params=dict(p0),
                                                  extra_kws=kws or None),
                      target_off=targets[0], target_on=targets[1], params=dict(p0), label=label, apl_bracket=rng,
                      off_extra=off_extra)
    return dict(hold=hold, kws=kws, orn=orn, breadth=Pp['breadth'], claw_thr=ct, targets=tuple(targets), p0=p0,
                off_extra=off_extra)


def simulate(setting: str, panel: str, wiring: str, orn: pd.DataFrame, lock: Mapping, lock134: Mapping,
             label: str = 'E147') -> dict:
    """One model fly's runs: the setting's model on the wiring with input `orn`, both rates re-held (`setting_model`'s
    start; `hold_with_input`) or, for a free setting, at its lab-input reading's parameters
    (`own_params`); built at those parameters and both conditions run with every stage saved. Returns dict(out:
    E134's `run_stages` output, built, params, params_off, hold: hold_rates' dict or None, kws: the model kwargs
    beside the family's, targets: (target_off, target_on, fraction_off, fraction_on), arm, z, start, p_start: the
    parameters the hold started from (the free reading's own), free). RuntimeError where the run differs from the
    hold's outputs."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    base, free = split_setting(setting)
    m = setting_model(base, wiring, lock)
    with e136() as F:
        if free:
            own = own_params(base, panel, wiring)
            if own is None:
                raise FileNotFoundError(f"{setting} {panel} {wiring}: the setting's reading at the lab's input is not "
                                        'saved')
            p, poff = dict(own['params']), dict(own['params_off'])
            kws = dict(F.arm_kws(m['arm']), **(F.a1_kws(wiring, panel, orn, m['z']) if base == 'A1' else {}))
            targets = tuple(F.read_targets(wiring, panel, orn.shape[1], lock134))
            hold, start, p_start = None, own['source'], dict(p)
        else:
            h = hold_with_input(m['arm'], m['z'], wiring, panel, orn, m['p0'], lock134, label=label)
            hold, kws, targets = h['hold'], dict(h['kws']), tuple(h['targets'])
            p, poff = dict(hold['params']), dict(hold['params_off'])
            start, p_start = m['start'], dict(h['p0'])
        built = B.build_family(wiring, panel, p, orn=orn, extra_kws=kws or None)
        out = B.run_stages(B.Point(F.PANEL, 'held', wiring), built, p, poff, claw_reader=F.claw_reader_of(m['arm']))
    if hold is not None:
        for c in R.CONDITIONS:
            if B._differ(out[c]['spike_counts'], hold[c][1]) or B._differ(out[c]['responses'], hold[c][0]):
                raise RuntimeError(f"{label} {c}: the saved run differs from the hold's outputs")
    return dict(out=out, built=built, params=p, params_off=poff, hold=hold, kws=kws, targets=targets, arm=m['arm'],
                z=dict(m['z']), start=start, p_start=p_start, free=free)


def _ratio(a, b) -> float:
    try:
        return float(a) / float(b) if float(b) != 0.0 else float('nan')
    except (TypeError, ValueError):
        return float('nan')


def read_row(tag: str, step: int, setting: str, panel: str, fly: str, wiring: str, sim: Mapping, lane: str,
             seconds: float) -> dict:
    """A read's row: atlas's READ_COLUMNS first, in its order ('held' atlas's `is_held`, False for the free reading),
    then this runner's: the reading ('held' / 'free'), held_exact and the tolerance, E136's arm and A1's coordinates
    (z_<name>) and normalization, where the hold started (start, start_fixed_thr, start_wAPLKC_scale), the setting's
    reading at the lab's input (own_source, own_fixed_thr, own_wAPLKC_scale, thr_vs_own, apl_vs_own), the APL-off /
    APL-on rate ratio, the other parameters and the commits."""
    from experiments.E77_input_screen import frame_sha256
    M = atlas()
    base, free = split_setting(setting)
    out, p, poff, hold = sim['out'], sim['params'], sim['params_off'], sim['hold']
    n_on, n_off = n_responses(out['on']['responses']), n_responses(out['off']['responses'])
    t_off, t_on, f_off, f_on = sim['targets']
    if hold is not None:
        thr_f, apl_f = float(hold['thr_factor']), float(hold['apl_factor'])
        b_thr, b_apl, n_sims = bool(hold['at_bound_thr']), bool(hold['at_bound_apl']), int(hold['n_sims'])
        held_exact, held = held_flags(n_on, n_off, int(t_on), int(t_off), b_thr or b_apl)
    else:
        thr_f, apl_f, b_thr, b_apl, n_sims, held_exact, held = 1.0, 1.0, False, False, 2, False, False
    rate = {c: float((out[c]['responses'].to_numpy(dtype=float) > 0.5).mean()) for c in ('on', 'off')}
    vals = dict(tag=tag, step=int(step), setting=setting, panel=panel, fly=fly, wiring=wiring,
                **{k: float(p[k]) for k in M.PARAM_KEYS}, **{f'off_{k}': float(poff[k]) for k in M.PARAM_KEYS},
                thr_factor=thr_f, apl_factor=apl_f, n_on=n_on, n_off=n_off, target_on=int(t_on), target_off=int(t_off),
                frac_on_target=float(f_on), frac_off_target=float(f_off), at_bound_thr=b_thr, at_bound_apl=b_apl,
                held=held, rate_on=rate['on'], rate_off=rate['off'],
                max_count_on=float(out['on']['spike_counts'].to_numpy(dtype=float).max()),
                max_count_off=float(out['off']['spike_counts'].to_numpy(dtype=float).max()),
                input_sha256=frame_sha256(sim['built'].orn), hold_n_sims=n_sims, seconds=round(float(seconds), 1),
                host=_host(), lane=lane, time=_now(), threads=os.environ.get('OMP_NUM_THREADS'), **commits())
    missing = [c for c in M.READ_COLUMNS if c not in vals]
    if missing:
        raise RuntimeError(f"the row lacks atlas's READ_COLUMNS {missing}")
    row = {c: vals[c] for c in M.READ_COLUMNS}
    if not free and M.is_held(row) != held:
        raise RuntimeError(f"{tag}: the runner's held ({held}) is not atlas's is_held")
    own = own_params(base, panel, wiring)
    extra = dict(reading='free' if free else 'held', held_exact=held_exact, hold_tol=HOLD_TOL, arm=sim['arm'])
    extra.update({f'z_{k}': float(v) for k, v in sim['z'].items()})
    extra.update({k: float(v) for k, v in sim['kws'].items() if k in ('pn_inhadd', 'pn_inhsc')})
    extra.update(start=sim['start'], start_fixed_thr=float(sim['p_start']['fixed_thr']),
                 start_wAPLKC_scale=float(sim['p_start']['wAPLKC_scale']),
                 own_source=None if own is None else own['source'],
                 own_fixed_thr=np.nan if own is None else float(own['params']['fixed_thr']),
                 own_wAPLKC_scale=np.nan if own is None else float(own['params']['wAPLKC_scale']),
                 thr_vs_own=np.nan if own is None else _ratio(p['fixed_thr'], own['params']['fixed_thr']),
                 apl_vs_own=np.nan if own is None else _ratio(p['wAPLKC_scale'], own['params']['wAPLKC_scale']),
                 rate_ratio=_ratio(rate['off'], rate['on']), thread_limit=os.environ.get('OMP_THREAD_LIMIT'))
    extra.update({k: float(v) for k, v in p.items() if k not in M.PARAM_KEYS and np.isscalar(v)})
    extra.update({f'off_{k}': float(v) for k, v in poff.items() if k not in M.PARAM_KEYS and np.isscalar(v)})
    extra.update({k: v for k, v in vals.items() if k not in row})
    row.update({k: v for k, v in extra.items() if k not in row})
    return row


def read_one(tag: str, step: int, setting: str, panel: str, fly: str, wiring: str, orn_of: Callable[[], pd.DataFrame],
             lane: str, root: Path, dirs: Sequence[Path], reclaim: bool = False, done: Optional[set] = None,
             lock: Optional[Mapping] = None, label: str = 'E147') -> Optional[dict]:
    """One read, claimed under `root` (`claim`; None if done, i.e. in `done`, or claimed elsewhere): `simulate` on the
    input `orn_of()` (called once the claim is held), the run saved in each of `dirs` (E134's step-3 layout, with the
    read's step, setting, fly, start and hold in its info), its row (`read_row`) appended to <root>/reads_<lane>.csv
    and its provenance to <root>/provenance_reads_<lane>.json. The claim is released on failure; the read's simulators
    are dropped and collected before it returns."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    from experiments.E77_input_screen import frame_sha256
    from kcq import stage_readings as R
    root = Path(root)
    if done is not None and tag in done:
        return None
    if not claim(tag, lane, root, reclaim):
        print(f'{label} {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    sim = None
    try:
        lock = atlas().require_lock() if lock is None else lock
        lock134 = R.require_lock()
        orn = orn_of()
        sim = simulate(setting, panel, wiring, orn, lock, lock134, label=f'{label} {tag}')
        hold = sim['hold']
        hold_rec = None if hold is None else {k: hold.get(k) for k in ('thr_factor', 'apl_factor', 'n_on', 'n_off',
                                                                        'target_on', 'target_off', 'at_bound_thr',
                                                                        'at_bound_apl', 'n_sims', 'history')}
        info = dict(model=B.FAMILY, wiring=wiring, panel=panel, reading='free' if sim['free'] else 'read',
                    experiment=EXPERIMENT, step=int(step), setting=setting, fly=fly, arm=sim['arm'], z=sim['z'],
                    params=sim['params'], params_off=sim['params_off'], start=sim['start'], p_start=sim['p_start'],
                    targets=list(sim['targets']), hold=hold_rec, tag=tag, host=_host(), lane=lane, time=_now(),
                    input_sha256=frame_sha256(orn))
        R.save_step3_run(tag, sim['out'], sim['built'].meta, sim['built'].orn, info, dirs)
        row = read_row(tag, step, setting, panel, fly, wiring, sim, lane, time.time() - t0)
        F._append_row(root / f'reads_{lane}.csv', row)
        prov_path = root / f'provenance_reads_{lane}.json'
        prov = _read_json(prov_path) or dict(part='reads', lane=lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(info={k: v for k, v in sim['built'].info.items() if not str(k).startswith('_')},
                                   params=sim['params'], params_off=sim['params_off'], hold=hold_rec,
                                   model=H.model_provenance(orn_deltas=sim['built'].orn))
        _write_json(prov_path, prov)
        print(f"{label} {tag}: {row['seconds']:.0f} s, n on {row['n_on']} / {row['target_on']}, off {row['n_off']} / "
              f"{row['target_off']}, held {row['held']} (exact {row['held_exact']})", flush=True)
        return row
    except BaseException:
        B.release(tag, root=root)
        raise
    finally:
        sim = None
        gc.collect()


# ---------------------------------------------------------------------------
# the input checks and the identity

def input_checks(panel: str) -> dict:
    """The input flies on a panel before any read (atlas's functions): the zero-deviation input fly is the lab's input
    exactly (values, index, columns); every input fly has the lab input's glomeruli and odors and no NaN; the mean over
    the nine is the lab's input within INPUT_MEAN_TOL; on APL16 the odors kept at the lab's input (paa, euc, -aPine,
    B-myr) are the lab's exactly in every fly. RuntimeError otherwise. Returns the checks' record (each fly's input
    sha256 and number of entries that differ from the lab's)."""
    from experiments.E77_input_screen import frame_sha256
    M = atlas()
    lab = M.lab_input(panel)
    zero = lab_of(panel)
    if not (zero.index.equals(lab.index) and zero.columns.equals(lab.columns)
            and np.array_equal(zero.to_numpy(dtype=float), lab.to_numpy(dtype=float))):
        raise RuntimeError(f'{panel}: the zero-deviation input fly is not the lab input')
    ins = {f: input_of(f, panel) for f in M.input_flies(panel)}
    for f, x in ins.items():
        if not (x.index.equals(lab.index) and x.columns.equals(lab.columns)):
            raise RuntimeError(f"{panel} {f}: the input's glomeruli or odors are not the lab input's")
        if not np.isfinite(x.to_numpy(dtype=float)).all():
            raise RuntimeError(f'{panel} {f}: NaN or inf in the input')
    L = lab.to_numpy(dtype=float)
    mean_diff = float(np.abs(np.mean([x.to_numpy(dtype=float) for x in ins.values()], axis=0) - L).max())
    if mean_diff > INPUT_MEAN_TOL:
        raise RuntimeError(f"{panel}: the mean over the input flies differs from the lab's input by {mean_diff:.3g}")
    kept = []
    if panel == 'apl16':
        kept = [j for j, c in enumerate(lab.columns) if H.odor_abbrev(c) in M.APL16_KEPT_AT_LAB]
        if len(kept) != len(M.APL16_KEPT_AT_LAB):
            raise RuntimeError(f'apl16: {len(kept)} of the kept odors {M.APL16_KEPT_AT_LAB} in the lab input')
        for f, x in ins.items():
            if not np.array_equal(x.to_numpy(dtype=float)[:, kept], L[:, kept]):
                raise RuntimeError(f"apl16 {f}: a kept odor differs from the lab's input")
    return dict(panel=panel, lab_sha256=frame_sha256(lab), zero_sha256=frame_sha256(zero), max_mean_diff=mean_diff,
                mean_tol=INPUT_MEAN_TOL, kept_odors=[H.odor_abbrev(lab.columns[j]) for j in kept],
                flies={f: dict(sha256=frame_sha256(x), n_entries_changed=int((x.to_numpy(dtype=float) != L).sum()))
                       for f, x in ins.items()})


def identity_tag(setting: str, panel: str, wiring: str) -> str:
    """'identity_<setting>_<panel>_<wiring>'."""
    return f'identity_{setting}_{panel}_{wiring}'


def identity_root() -> Path:
    """results/E147_megamat_input_flies/identity/ (created): the identity's claims, rows, provenance and records."""
    d = atlas().results_dir() / IDENTITY_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def identity_dirs(setting: str, panel: str, wiring: str) -> List[Path]:
    """Where an identity read's run is saved: outputs/E147_megamat_input_flies/identity/<setting>/<panel>/<wiring>/."""
    return [atlas().local_dir(f'{IDENTITY_DIR}/{setting}/{panel}/{wiring}')]


def _table(root: Path) -> pd.DataFrame:
    fs = [pd.read_csv(f) for f in sorted(Path(root).glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def frame_diff(a: pd.DataFrame, b: pd.DataFrame) -> int:
    """How many entries of two frames differ, their odor columns abbreviated (-1 where their shapes, indices or columns
    differ)."""
    from kcq import stage_readings as R
    a, b = R._abbrev(a), R._abbrev(b)
    if a.shape != b.shape or not a.index.equals(b.index) or list(a.columns) != list(b.columns):
        return -1
    return int((a.to_numpy(dtype=float) != b.to_numpy(dtype=float)).sum())


def identity_compare(tag: str, setting: str, panel: str, wiring: str, dirs: Sequence[Path]) -> dict:
    """An identity read's saved run against the setting's reading at the lab's input (`lab_reading`): the differing
    entries of each condition's spike counts and responses, of every stage matrix both runs hold, and the largest
    relative difference of the parameters; exact when every count is 0."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    run = R.load_step3_run(tag, dirs=dirs)
    ref = lab_reading(setting, panel, wiring)
    counts, frames = {}, {}
    for c in R.CONDITIONS:
        for i, k in enumerate(('responses', 'spike_counts')):
            counts[f'{c}/{k}'] = frame_diff(run.conditions[c][k], ref['outputs'][c][i])
        if ref['run'] is not None:
            for k in sorted((set(run.conditions[c]) & set(ref['run'].conditions[c])) - {'responses', 'spike_counts'}):
                frames[f'{c}/{k}'] = frame_diff(run.conditions[c][k], ref['run'].conditions[c][k])
    rel = {k: abs(float(run.info['params'][k]) / float(v) - 1.0) if float(v) != 0.0
           else abs(float(run.info['params'][k])) for k, v in ref['params'].items() if k in run.info['params']}
    rec = dict(tag=tag, setting=setting, panel=panel, wiring=wiring, reference=ref['source'], counts=counts,
               frames=frames, params_max_rel_diff=max(rel.values()) if rel else None,
               input_sha256=run.info.get('input_sha256'))
    if (setting, panel) == ('family', 'megamat'):
        # recorded beside, not gating: E134's hash is of its input as loaded then (same frame_sha256)
        t = B.runs_table()
        t = t[t.tag.astype(str) == B.tag_of(B.Point('megamat', 'held', wiring))]
        rec['e134_input_sha256'] = str(t.input_sha256.iloc[0]) if len(t) else None
        rec['input_as_e134'] = rec['e134_input_sha256'] == rec['input_sha256'] if len(t) else None
    rec['exact'] = bool(all(v == 0 for v in counts.values()) and all(v == 0 for v in frames.values()))
    return rec


def identity_records() -> Dict[str, dict]:
    """Every lane's identity records (identity/identity_*.json), by tag; a tag recorded as not exact by any lane stays
    not exact."""
    out: Dict[str, dict] = {}
    for f in sorted(identity_root().glob('identity_*.json')):
        for tag, rec in _read_json(f).get('records', {}).items():
            if tag not in out or not rec.get('exact'):
                out[tag] = rec
    return out


def require_identity() -> Dict[str, dict]:
    """The identity records; RuntimeError unless the family's megamat identity is exact on at least one wiring and no
    identity has failed."""
    recs = identity_records()
    bad = [t for t, r in recs.items() if not r.get('exact')]
    if bad:
        raise RuntimeError(f'the identity fails on {bad} (results/.../identity/): no read before it is resolved')
    if not any(r.get('setting') == 'family' and r.get('panel') == 'megamat' and r.get('exact') for r in recs.values()):
        raise RuntimeError('no exact identity of the family on megamat yet: run --part identity first')
    return recs


def part_identity(lane: Optional[str] = None, wirings: Sequence[str] = IDENTITY_WIRINGS,
                  settings: Sequence[str] = ('family',), panels: Sequence[str] = ('megamat',),
                  reclaim: bool = False) -> dict:
    """The input checks on each panel, then the zero-deviation input fly read on each wiring for each setting
    (`read_one` under identity/), compared with the setting's reading at the lab's input (`identity_compare`) and
    recorded in identity/identity_<lane>.json (with the input checks). A point already recorded exact is skipped, one
    claimed by another lane is left to it, and one already read but not recorded is compared from its saved run.
    RuntimeError on any difference."""
    M = atlas()
    lane = lane_name(lane)
    lock = M.require_lock()
    root = identity_root()
    path = root / f'identity_{lane}.json'
    rec = _read_json(path) or dict(lane=lane, host=_host(), script=SCRIPT, records={}, input_checks={})
    done = identity_records()
    out: Dict[str, dict] = {}
    for panel in panels:
        rec['input_checks'][panel] = dict(input_checks(panel), time=_now())
        _write_json(path, rec)
        for setting in settings:
            if setting not in M.HELD_SETTINGS:
                raise ValueError(f'the identity reads a re-held setting ({M.HELD_SETTINGS}), not {setting!r}')
            for w in wirings:
                tag = identity_tag(setting, panel, w)
                if done.get(tag, {}).get('exact'):
                    out[tag] = done[tag]
                    continue
                dirs = identity_dirs(setting, panel, w)
                if tag not in set(_table(root).tag.astype(str)):
                    row = read_one(tag, 0, setting, panel, IDENTITY_FLY, w, lambda p=panel: lab_of(p), lane, root,
                                   dirs, reclaim=reclaim, lock=lock, label='E147 identity')
                    gc.collect()
                    if row is None:
                        continue
                r = dict(identity_compare(tag, setting, panel, w, dirs), host=_host(), lane=lane, time=_now(),
                         **commits())
                rec['records'][tag] = r
                _write_json(path, rec)
                out[tag] = r
                print(f"E147 identity {tag}: exact {r['exact']}, counts {r['counts']}", flush=True)
                if not r['exact']:
                    raise RuntimeError(f'E147 identity fails: {r}')
    return out


# ---------------------------------------------------------------------------
# the reads

def planned(step: int, free: bool = False, reverse: bool = False) -> list:
    """atlas's `planned_points(step, free)` in a lane's order (`reverse`: from the other end)."""
    pts = list(atlas().planned_points(step, free))
    return pts[::-1] if reverse else pts


def part_reads(step: int, lane: Optional[str] = None, free: bool = False, reverse: bool = False,
               reclaim: bool = False, max_reads: Optional[int] = None, archive: bool = True) -> dict:
    """One lane of a step's reads (`planned`), each claimed, skipped when done or claimed elsewhere, its input atlas's
    `input_fly`, saved in atlas's `run_dirs` (local only without `archive`), its row in reads_<lane>.csv; after the
    lock, the identity (`require_identity`) and the input checks of the step's panels. `max_reads` stops the lane after
    that many reads of its own. Returns a summary (reads made here, the unheld ones by the lock's rule, those whose
    counts are not exactly at their targets)."""
    M = atlas()
    lane = lane_name(lane)
    lock = M.require_lock()
    require_identity()
    pts = planned(step, free, reverse)
    for panel in sorted({pt.panel for pt in pts}):
        input_checks(panel)
    root = M.results_dir()
    made, unheld, inexact = [], [], []
    for pt in pts:
        tag = M.point_tag(pt.setting, pt.panel, pt.fly, pt.wiring)
        dirs = M.run_dirs(pt.step, pt.setting, pt.wiring, pt.fly)
        row = read_one(tag, pt.step, pt.setting, pt.panel, pt.fly, pt.wiring,
                       lambda f=pt.fly, p=pt.panel: input_of(f, p), lane, root, dirs if archive else dirs[:1],
                       reclaim=reclaim, done=set(M.reads_table().tag.astype(str)), lock=lock)
        gc.collect()
        if row is None:
            continue
        made.append(tag)
        if not free and not row['held']:
            unheld.append(tag)
        if not free and not row['held_exact']:
            inexact.append(tag)
        if max_reads is not None and len(made) >= max_reads:
            break
    tab = M.reads_table()
    want = {M.point_tag(p.setting, p.panel, p.fly, p.wiring) for p in pts}
    summary = dict(step=int(step), free=bool(free), lane=lane, made_here=len(made), unheld_here=unheld,
                   not_exact_here=inexact, done_of_step=int(tab.tag.astype(str).isin(want).sum()), planned=len(pts))
    print(json.dumps(summary, indent=1), flush=True)
    return summary


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('identity', 'reads'))
    ap.add_argument('--step', type=int, choices=(1, 2), help='--part reads: the step')
    ap.add_argument('--free', action='store_true', help="--part reads: the step's points under the free reading")
    ap.add_argument('--lane', help="'<host>' (default) or '<host>-<n>'")
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true', help="retake this lane's own stale claims")
    ap.add_argument('--max-reads', type=int)
    ap.add_argument('--no-archive', action='store_true', help='save the runs locally only')
    ap.add_argument('--wirings', nargs='+', default=list(IDENTITY_WIRINGS), help='--part identity: the wirings')
    ap.add_argument('--settings', nargs='+', default=['family'], choices=('family', 'A1'),
                    help='--part identity: the settings')
    ap.add_argument('--panels', nargs='+', default=['megamat'], choices=('megamat', 'apl16'),
                    help='--part identity: the panels')
    a = ap.parse_args(argv)
    for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT'):
        if k not in os.environ:
            print(f'E147: {k} is not set (CLAUDE.md, "Gotchas": rastermap resets OpenMP to every core without '
                  'OMP_THREAD_LIMIT)', file=sys.stderr, flush=True)
    if a.part == 'identity':
        print(json.dumps(part_identity(a.lane, a.wirings, a.settings, a.panels, a.reclaim), indent=1, default=str))
    else:
        if a.step is None:
            raise SystemExit('--part reads needs --step 1 or 2')
        part_reads(a.step, a.lane, a.free, a.reverse, a.reclaim, a.max_reads, not a.no_archive)


if __name__ == '__main__':
    main()
