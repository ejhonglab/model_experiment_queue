#!/usr/bin/env python3
"""E133: what in a calyx the best family's standing depends on: the male CNS's left calyx against its own re-drawn
wirings, and one difference from the hemibrain moved at a time.

The entry is item E133 of experiment_priorities.md (moved in by viki 2026-11-24 05:51 under the idle-host rule; planned
by a Fable subagent, claude-fable-5-2). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight
wirings (E83's `saved_fit`). "Fixed scales" is E130's step 1: the hemibrain fit's (C2's) factors on the calyx's own base
(E130's bases.json), KC-KC at C2's scale on the calyx's own matrix (the male calyces' from Pratyush's KCclaw-to-KCclaw
table, E130's matrices; the hemibrain's olfsysm's), both response rates re-held (E77_input_screen.hold_rates) at C2's
response counts as fractions of the calyx's KCs x odors. The calyces: the hemibrain under multibouton_claws 'argmax'
and the male CNS's left and right (E10_joint_fit.py's `--connectome`).

This module holds step 0, the identity check, step 1 (the 21 re-draws) and its score, and step 2's runner. Step 2's
transplant builders are experiments/E133_builders.py, which only `part_step2` imports: steps 0 and 1 never load it, so
it can be committed while step-1 lanes run without tripping E108's `_still_clean` guard.

Parts (CLI `--part`):
- step0 (no simulation): (i) from E130's saved runs of the three calyces at fixed scales ('S1_<calyx>') and of their
  counted fits, and of the family's eight wirings: the calibrated pattern part over all 120 pairs and over the pairs that
  hold 2-but, paa or B-myr and the others (`pair_groups`; the module's readings below), against all recordings and each
  of the four halves A to D (E108's `effect_targets`); for the 120 also the amplitude and shape parts, the spread ratio
  and Pearson r (E108's `unit_scores`); the planner's numbers beside. (ii) The three odors' responders with APL on and
  off, by KC type and claw count, and their overlap. (iii) From the model's tables, for each of the three odors'
  glomeruli (those carrying at least a tenth of the odor's summed positive ORN input) on each calyx: the claws by KC
  type, the KCs with two or more such claws, the boutons, and APL's synapses on those claws and boutons. step0.json,
  step0_*.csv, step0_table.md.
- check (simulations): the identity. For each calyx, with no transplant and no re-draw: (i) the runner's tables
  (`calyx_tables`) have E130's wiring and breadth, and `breadth_run` (the path step 2 takes) reproduces that breadth;
  (ii) at E130's held parameters the runner (`build_simulator`) reproduces E130's fixed-scale run S1_<calyx> with 0
  differing spike counts and responses, APL on and off (`replay`); (iii) the step-1 hold from E130's start reaches the
  same held factors and parameters (to a relative 1e-12: the bisection's steps can differ in the last bit between
  CPUs) and the same outputs (0 differing counts and responses). check_<host>.json.
- step1 (simulations, resumable, lanes `--calyces`, `--seeds`, `--reverse`; points claimed by mkdir of
  results/E133_calyx_transplants/claims/<tag>/, `--reclaim` to take over a dead process's): the 21 re-draws at fixed
  scales, seeds 1 to 7 of the within-type claw permutation (pyolfsysm's `shuffle_wPNKC` through E10_joint_fit.bouton_base's
  `wiring='claw_permutation_within_type'`, E72's P draws) of each calyx, E130's step 1 on each (the breadth recomputed on
  the permuted wiring; each claw keeps its KC, its APL weights and its KC's KC-KC partners). A lane first replays the
  identity (ii) of its calyces on its host (replay_<host>.json) and refuses on a difference. Rows in
  runs_1_<host>_<lane>.csv, provenance in provenance_1_<host>_<lane>.json, the runs locally
  (outputs/E133_calyx_transplants/step1/) and in the archive (data/model_outputs/E133_calyx_transplants/step1/), and
  each point's raw-count pair distances in results/.../pairs/<tag>.csv.
- score1 (no simulation): L, H and step 1's class; the right calyx's seven beside (no class); every pattern part with
  its amplitude part, shape part, spread ratio, Pearson r and the grouped errors; the realism columns. readings1.csv,
  runs1.csv, verdict1.json, step1_table.md.
- variants (no simulation): variants.csv and model_pair_distances_long.csv for kcq/apl_effect.py and the archive.
- step2 (simulations; only on step 1's "kept" or "mixed"): the six arms on each host's four wirings, through
  experiments/E133_builders.py.
- commands: the lanes' command lines.

Readings the entry leaves open, fixed here before any run (the reading closest to E130's code):
- The pair groups. The planner's table names "the 39 pairs that hold 2-but, paa or B-myr" and "the other 81", but its
  numbers are those of the 42 pairs that hold any of the three and the 78 that hold none (step 0 reproduces them to
  every decimal given, on all recordings and on the halves; the 3 x 13 = 39 pairs of exactly one of the three, with the
  3 pairs among them counted in the other 81, do not). The groups read are therefore 'p42' / 'p78' (the planner's
  numbers), with 'p39' / 'p81' (the counts as written) and 'p3' (the 3 pairs among the three odors) beside.
- A group's pattern error: the per-pair APL effect's pattern error over all 120 pairs (model minus data, each with its
  mean over the 120 removed, E83's `pattern_error`), root mean square over the group's pairs; so the pattern part
  squared is the pair-count-weighted mean of its two groups' squares.
- The "no class" guard of a step-1 point: its largest spike count per KC and odor above twice its own unpermuted
  calyx's (E130's S1_<calyx>) in either condition; or outside a band (E110's `bands_of`, E130's five); or at a corner (a
  bisection at its bracket's bound, or the held threshold outside 1/4 to 4 times the base's, or APL->KC outside 1/30 to
  30 times the start's: E124's bounds as E130's `bound_flags` reads them).
- Step 1's class is read on all recordings, the halves A to D beside. A class needs every one of the 14 runs of L and
  H to have one (E130's `class_population`); otherwise no class, with the would-be class and the runs named. "The
  hemibrain's median" is the median of the hemibrain's seven re-draws.

Holds (CLAUDE.md, "Connectome data"): the hemibrain and Pratyush's male CNS only (olfsysm's model tables and the shared
directory's data/Male-CNS_Tom/). No olfsysm change. Bouton identities are never inferred from position or tuning. No
recording outside olfdb is read. Nothing is published.

Concentrations: the model's input for paa, euc and -aPine is the ORN data within 1 log10 unit of the flies'.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E133_calyx_transplants.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E133_calyx_transplants.py --part check
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E133_calyx_transplants.py --part step1 --calyces malecns-left hemibrain-argmax malecns-right
    (the lanes: --part commands; a second lane with --reverse takes the seeds from 7 down)
    ./venv/bin/python -u experiments/E133_calyx_transplants.py --part score1 && \\
        ./venv/bin/python -u experiments/E133_calyx_transplants.py --part variants
    ./venv/bin/python -u experiments/E133_calyx_transplants.py --part commands
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import shutil
import socket
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E133_calyx_transplants'
SCRIPT = 'experiments/E133_calyx_transplants.py'
E130_EXP = 'E130_malecns_family'
LEFT, HEMI, RIGHT = 'malecns-left', 'hemibrain-argmax', 'malecns-right'
CALYCES = (HEMI, LEFT, RIGHT)
MALE = (LEFT, RIGHT)
CODE = {HEMI: 'hb', LEFT: 'ml', RIGHT: 'mr'}
SHORT = {HEMI: 'hemibrain (argmax)', LEFT: 'male CNS left', RIGHT: 'male CNS right'}
HOSTS = (HEMI, LEFT)
DONOR = {HEMI: LEFT, LEFT: HEMI}
SEEDS1 = tuple(range(1, 8))
HOST_SEEDS = (0, 1, 2, 3)
WIRING = 'claw_permutation_within_type'
THREE = ('2-but', 'paa', 'B-myr')
GROUPS = ('p42', 'p78', 'p39', 'p81', 'p3')
GLOM_SHARE = 0.1
TARGETS_READ = ('all', 'A', 'B', 'C', 'D')
READOUTS = ('raw', 'cal')
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
K = 2.0
PARAM5 = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')
PARAM_KEYS = PARAM5 + ('wKCKC_scale',)
# step 1's class (the entry's, fixed beforehand)
KEPT_GAP = 0.007
REDEALT_GAP = 0.004
KEPT_ABOVE = 5
MAX_COUNT_FACTOR = 2.0
BAND_COLS = ('band_rate_on', 'band_ratio', 'band_rate_rho_on', 'band_opk_on', 'band_opk_off')
# the planner's read of E130's archive (the entry's table): (120 pairs, the three odors' pairs, the others)
PLANNER = {'S1:hemibrain-argmax': (0.0576, 0.051, 0.061), 'S1:malecns-left': (0.0719, 0.090, 0.060),
           'S1:malecns-right': (0.0625, 0.063, 0.062), 'fit:hemibrain-argmax': (0.0586, 0.048, 0.064),
           'fit:malecns-left': (0.0715, 0.092, 0.057), 'fit:malecns-right': (0.0608, 0.060, 0.061)}
PLANNER_FAMILY = {'p120': (0.047, 0.064), 'three': (0.041, 0.061), 'rest': (0.048, 0.070)}
_CACHE: dict = {}

Point = namedtuple('Point', 'step calyx seed arm')     # arm None in step 1


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E133_calyx_transplants/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


REL_TOL = 1e-12


def _rel(a: float, b: float) -> float:
    """|a - b| relative to the larger magnitude (0 when both are 0)."""
    a, b = float(a), float(b)
    m = max(abs(a), abs(b))
    return 0.0 if m == 0 else abs(a - b) / m


def _nanmean(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(x[np.isfinite(x)].mean()) if np.isfinite(x).any() else float('nan')


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2))) if x.size else float('nan')


def sha256_frame(df: pd.DataFrame) -> str:
    return hashlib.sha256(np.ascontiguousarray(df.to_numpy(dtype=float), dtype='<f8').tobytes()).hexdigest()


def tag_of(pt: Point) -> str:
    """'S1_<calyx>_seed<n>' (step 1); 'S2_<arm>_<host>_seed<n>' (step 2)."""
    if pt.step == 1:
        if pt.arm is not None:
            raise ValueError(f'a step-1 point has no arm: {pt}')
        return f'S1_{pt.calyx}_seed{int(pt.seed)}'
    if pt.step == 2:
        return f'S2_{pt.arm}_{pt.calyx}_seed{int(pt.seed)}'
    raise ValueError(f'step {pt.step}')


def parse_tag(tag: str) -> Point:
    """The point of a tag (`tag_of`'s inverse)."""
    head, seed = tag.rsplit('_seed', 1)
    if head.startswith('S1_'):
        return Point(1, head[3:], int(seed), None)
    if head.startswith('S2_'):
        arm, c = head[3:].split('_', 1)
        return Point(2, c, int(seed), arm)
    raise ValueError(tag)


def lane_name(calyces: Sequence[str], seeds: Sequence[int] = SEEDS1, reverse: bool = False) -> str:
    """'<codes>[-s<seeds>][-rev']: e.g. 'ml-hb-mr', 'ml-s1-2-3-rev'."""
    s = '' if tuple(seeds) == SEEDS1 else '-s' + '-'.join(str(int(x)) for x in seeds)
    return '-'.join(CODE[c] for c in calyces) + s + ('-rev' if reverse else '')


def odor_of(col) -> str:
    """'2h @ -3' -> '2h'."""
    return str(col).split(' @ ')[0]


# ---------------------------------------------------------------------------
# the pair groups and the readings

def pair_groups(pairs: pd.Index, odors: Sequence[str] = THREE) -> Dict[str, np.ndarray]:
    """Boolean masks over pairs (an index of (odor_a, odor_b)): 'p120' (all), 'p42' (any of `odors`), 'p78' (none),
    'p39' (exactly one of `odors`), 'p81' (not exactly one: none, or both), 'p3' (both in `odors`). On the 16 odors
    with the three: 120, 42, 78, 39, 81, 3 pairs."""
    S = set(odors)
    a = np.array([x[0] in S for x in pairs])
    b = np.array([x[1] in S for x in pairs])
    return {'p120': np.ones(len(pairs), dtype=bool), 'p42': a | b, 'p78': ~(a | b), 'p39': a ^ b, 'p81': ~(a ^ b),
            'p3': a & b}


def grouped_errors(err: pd.Series, odors: Sequence[str] = THREE) -> Dict[str, float]:
    """The root mean square of a per-pair pattern error over each group of `pair_groups` (keyed as the groups)."""
    m = pair_groups(err.index, odors)
    v = err.to_numpy(dtype=float)
    return {g: _rms(v[mk]) for g, mk in m.items()}


def group_patterns(d: Dict[str, pd.Series], targets: Optional[Mapping[str, pd.Series]] = None) -> dict:
    """Per effect target t (default E108's 'all' and halves A to D): '<t>_pattern_<group>' for the groups of
    `pair_groups` (the module docstring's reading: the pattern error over the 120 pairs, rms over the group's), and on
    the 3 pairs among the three odors the model's and the data's per-pair APL effect ('p3_effect_model_<a>_<b>',
    '<t>_p3_effect_data_<a>_<b>')."""
    from experiments import E83_cross_wiring_fit as E83
    if targets is None:
        from experiments import E108_paa_extrapolation as E108
        targets = E108.effect_targets()
    me = E83._cs(H.apl_effect(d['on'], d['off']))
    m3 = pair_groups(me.index)['p3']
    row = {}
    for pr, v in me[m3].items():
        row[f'p3_effect_model_{pr[0]}_{pr[1]}'] = float(v)
    for t, de0 in targets.items():
        de = E83._cs(de0).reindex(me.index)
        if de.isna().any():
            raise RuntimeError(f'target {t}: pairs missing')
        err = (me - me.mean()) - (de - de.mean())
        for g, v in grouped_errors(err).items():
            row[f'{t}_pattern_{g}'] = v
        for pr, v in de[m3].items():
            row[f'{t}_p3_effect_data_{pr[0]}_{pr[1]}'] = float(v)
    return row


def unit_readings(d: Dict[str, pd.Series]) -> dict:
    """E110's `scores_of` (E108's `unit_scores` against all recordings and halves A to D: loss, pattern part,
    amplitude and shape parts, spread ratio, ρ and Pearson r on the per-pair APL effect; the two ρ's) and
    `group_patterns`."""
    from experiments import E110_no_kckc as E110
    row = E110.scores_of(d)
    row.update(group_patterns(d))
    return row


def _dist(o: Mapping[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83's `pair_distances`; 'cal' with K re-calibrated on the run's APL-on counts), cached
    in this process and in outputs/.../distances_cache.pkl by the sha256 of the run's spike counts and the readout."""
    import pickle
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist' not in _CACHE:
        try:
            _CACHE['dist'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist'] = {}
    st = _CACHE['dist']
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in st:
        st[key] = E83.pair_distances(dict(o), readout)
        _CACHE['dist_dirty'] = True
    return st[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


# ---------------------------------------------------------------------------
# E130's runs

def e130():
    from experiments import E130_malecns_family as E130
    return E130


def e130_runs_table() -> pd.DataFrame:
    """E130's step-1 rows (runs_1_*.csv), by tag."""
    return e130().runs_table().set_index('tag')


def e130_s1(c: str) -> Dict[str, tuple]:
    o = e130().load_run(f'S1_{c}')
    if o is None:
        raise RuntimeError(f"E130's S1_{c} not found (outputs or archive)")
    return o


def e130_fit(c: str) -> Tuple[Dict[str, tuple], str]:
    """E130's counted fit of a calyx (outputs, tag)."""
    E130 = e130()
    cnt = E130.counted(c)
    if cnt is None:
        raise RuntimeError(f"E130's fit of {c} is not counted yet")
    o = E130.load_fit_outputs(cnt['spec'])
    if o is None:
        raise RuntimeError(f"E130's fit {cnt['spec'].name} has no outputs here or in the archive")
    return o, E130.e10_tag(cnt['spec'])


def reference_max_counts(c: str) -> Tuple[int, int]:
    """The largest spike count per KC and odor of the calyx's own (unpermuted) fixed-scale run, APL on and off."""
    o = e130_s1(c)
    return int(o['on'][1].to_numpy(dtype=float).max()), int(o['off'][1].to_numpy(dtype=float).max())


# ---------------------------------------------------------------------------
# step 0

def step0_units() -> Dict[str, Dict[str, tuple]]:
    """{unit: outputs}: 'S1:<calyx>' (E130's fixed-scale runs), 'fit:<calyx>' (its counted fits), 'family:<w>' (the
    eight wirings' saved fits)."""
    from experiments import E83_cross_wiring_fit as E83
    u = {}
    for c in CALYCES:
        u[f'S1:{c}'] = e130_s1(c)
    for c in CALYCES:
        u[f'fit:{c}'] = e130_fit(c)[0]
    for w in E83.ALL_WIRINGS:
        u[f'family:{w}'] = E83.saved_fit(w)['saved']
    return u


def step0_patterns(units: Mapping[str, Dict[str, tuple]]) -> pd.DataFrame:
    """`unit_readings` per unit and readout."""
    rows = []
    for ro in READOUTS:
        for u, o in units.items():
            rows.append(dict(unit=u, readout=ro, **unit_readings(_dist(o, ro))))
    _save_dist_cache()
    return pd.DataFrame(rows)


def planner_check(rd: pd.DataFrame) -> pd.DataFrame:
    """The planner's numbers against step 0's (calibrated, all recordings): each of the entry's table cells, read with
    the groups 'p42' / 'p78' and, beside, 'p39' / 'p81', rounded to the planner's decimals."""
    cal = rd[rd['readout'] == 'cal'].set_index('unit')
    rows = []
    for u, (v120, v3, vr) in PLANNER.items():
        for (lab, col, want, dec) in (('120 pairs', 'all_pattern', v120, 4), ('three odors', 'all_pattern_p42', v3, 3),
                                      ('the others', 'all_pattern_p78', vr, 3), ('three odors (39)', 'all_pattern_p39', v3, 3),
                                      ('the others (81)', 'all_pattern_p81', vr, 3)):
            got = float(cal.loc[u, col])
            rows.append(dict(unit=u, cell=lab, column=col, planner=want, step0=got, decimals=dec,
                             reproduced=bool(round(got, dec) == round(want, dec))))
    fam = cal.loc[[u for u in cal.index if u.startswith('family:')]]
    for (lab, col, key, dec) in (('family range, 120', 'all_pattern', 'p120', 3), ('family range, three odors', 'all_pattern_p42', 'three', 3),
                                 ('family range, others', 'all_pattern_p78', 'rest', 3),
                                 ('family range, three odors (39)', 'all_pattern_p39', 'three', 3),
                                 ('family range, others (81)', 'all_pattern_p81', 'rest', 3)):
        lo, hi = float(fam[col].min()), float(fam[col].max())
        want = PLANNER_FAMILY[key]
        rows.append(dict(unit='family8', cell=lab, column=col, planner=f'{want[0]} to {want[1]}', step0=f'{lo:.4f} to {hi:.4f}',
                         decimals=dec, reproduced=bool(round(lo, dec) == want[0] and round(hi, dec) == want[1])))
    return pd.DataFrame(rows)


def kc_claws(c: str) -> pd.Series:
    """Claws per KC of a calyx's own wiring (E130's `step0_wiring_only`), indexed by KC ID."""
    W = e130().step0_wiring_only(c)
    return pd.Series(1, index=W.index.get_level_values(H.KC_ID_LEVEL)).groupby(level=0).sum()


def responders(units: Mapping[str, Dict[str, tuple]]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The three odors' responders (two or more spikes) in each calyx's fixed-scale run and fit, APL on and off:
    (by KC type and claw count, a summary by KC type, the overlap of the three odors' responder sets)."""
    det, summ, ov = [], [], []
    for u, o in units.items():
        if not u.startswith(('S1:', 'fit:')):
            continue
        c = u.split(':', 1)[1]
        claws = kc_claws(c)
        for cond in ('on', 'off'):
            r = o[cond][0]
            kc = r.index.get_level_values(H.KC_ID_LEVEL)
            typ = np.asarray(r.index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)
            nc = claws.reindex(kc).to_numpy()
            if np.isnan(nc.astype(float)).any():
                raise RuntimeError(f'{u}: KCs without a claw count')
            nc = nc.astype(int)
            R = pd.DataFrame(r.to_numpy(dtype=float) > 0.5, columns=[odor_of(x) for x in r.columns])
            sets = {od: set(np.flatnonzero(R[od].to_numpy())) for od in THREE}
            for od in THREE:
                x = R[od].to_numpy()
                for t in TYPES + ('all',):
                    m = np.ones(len(x), dtype=bool) if t == 'all' else typ == t
                    summ.append(dict(unit=u, calyx=c, condition=cond, odor=od, kc_type=t, n_kcs=int(m.sum()),
                                     n_responders=int(x[m].sum()), rate=float(x[m].mean()),
                                     share_of_responders=float(x[m].sum() / max(x.sum(), 1)),
                                     claws_per_kc_responders=float(nc[m & x].mean()) if (m & x).any() else np.nan,
                                     claws_per_kc_type=float(nc[m].mean())))
                    if t == 'all':
                        continue
                    for n in np.unique(nc[m]):
                        mm = m & (nc == n)
                        det.append(dict(unit=u, calyx=c, condition=cond, odor=od, kc_type=t, n_claws=int(n),
                                        n_kcs=int(mm.sum()), n_responders=int(x[mm].sum()), rate=float(x[mm].mean())))
            for i, a in enumerate(THREE):
                for b in THREE[i + 1:]:
                    inter, union = len(sets[a] & sets[b]), len(sets[a] | sets[b])
                    ov.append(dict(unit=u, calyx=c, condition=cond, odors=f'{a} & {b}', n_a=len(sets[a]), n_b=len(sets[b]),
                                   n_both=inter, jaccard=inter / union if union else np.nan,
                                   share_of_smaller=inter / min(len(sets[a]), len(sets[b])) if min(len(sets[a]), len(sets[b])) else np.nan))
            all3 = sets[THREE[0]] & sets[THREE[1]] & sets[THREE[2]]
            any3 = sets[THREE[0]] | sets[THREE[1]] | sets[THREE[2]]
            ov.append(dict(unit=u, calyx=c, condition=cond, odors='all three', n_both=len(all3), n_any=len(any3),
                           jaccard=len(all3) / len(any3) if any3 else np.nan))
    return pd.DataFrame(det), pd.DataFrame(summ), pd.DataFrame(ov)


def odor_glomeruli(orn: Optional[pd.DataFrame] = None, share: float = GLOM_SHARE,
                   odors: Sequence[str] = THREE) -> Dict[str, pd.Series]:
    """Each odor's glomeruli carrying at least `share` of its summed positive ORN input (the model's input,
    `H.apl16_orn_deltas()`): {odor: their inputs (descending)}."""
    orn = H.abbrev_columns(H.apl16_orn_deltas() if orn is None else orn)
    out = {}
    for od in odors:
        x = orn[od].astype(float)
        pos = x.clip(lower=0)
        tot = float(pos.sum())
        if not tot > 0:
            raise ValueError(f'{od}: no positive ORN input')
        out[od] = pos[pos >= share * tot].sort_values(ascending=False)
    return out


def glomerulus_tables(c: str, gloms: Mapping[str, pd.Series], W: Optional[pd.DataFrame] = None,
                      counts: Optional[Mapping[str, pd.Series]] = None) -> Tuple[List[dict], List[dict]]:
    """For each odor's glomeruli on a calyx (its own wiring `W` and APL's counted synapses `counts`, E130's
    `step0_wiring_only` and `apl_counts` by default): the claws by KC type, the KCs with two or more such claws, the
    boutons, and APL's synapses on those claws and boutons (sums, shares of the calyx's, per claw / bouton over the
    calyx's mean); (rows per odor, rows per odor and glomerulus)."""
    E130 = e130()
    W = E130.step0_wiring_only(c) if W is None else W
    a = E130.apl_counts(c, W) if counts is None else counts
    v = W.to_numpy(dtype=float)
    claw_glom = np.asarray(W.columns.get_level_values('glomerulus'))[v.argmax(axis=1)]
    has = v.max(axis=1) > 0
    kc = np.asarray(W.index.get_level_values(H.KC_ID_LEVEL))
    typ = np.asarray(W.index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)
    bglom = np.asarray(W.columns.get_level_values('glomerulus'))
    akc, kca = a['apl_kc'].to_numpy(dtype=float), a['kc_apl'].to_numpy(dtype=float)
    apn, pna = a['apl_pn'].to_numpy(dtype=float), a['pn_apl'].to_numpy(dtype=float)
    kc_type = pd.Series(typ, index=kc).groupby(level=0).first()
    rows, by = [], []
    for od, g in gloms.items():
        gs = set(g.index)
        cm = has & np.isin(claw_glom, list(gs))
        bm = np.isin(bglom, list(gs))
        per_kc = pd.Series(cm.astype(int), index=kc).groupby(level=0).sum()
        r = dict(calyx=c, odor=od, glomeruli=' '.join(g.index), n_glomeruli=len(gs),
                 glomeruli_in_wiring=int(len(gs & set(bglom))), n_claws=int(cm.sum()), share_claws=float(cm.mean()),
                 n_boutons=int(bm.sum()), share_boutons=float(bm.mean()),
                 kcs_2plus=int((per_kc >= 2).sum()), kcs_1plus=int((per_kc >= 1).sum()),
                 apl_kc=float(akc[cm].sum()), share_apl_kc=float(akc[cm].sum() / akc.sum()),
                 apl_kc_per_claw_rel=float(akc[cm].mean() / akc.mean()) if cm.any() else np.nan,
                 kc_apl=float(kca[cm].sum()), share_kc_apl=float(kca[cm].sum() / kca.sum()),
                 kc_apl_per_claw_rel=float(kca[cm].mean() / kca.mean()) if cm.any() else np.nan,
                 apl_pn=float(apn[bm].sum()), share_apl_pn=float(apn[bm].sum() / apn.sum()),
                 apl_pn_per_bouton_rel=float(apn[bm].mean() / apn.mean()) if bm.any() else np.nan,
                 pn_apl=float(pna[bm].sum()), share_pn_apl=float(pna[bm].sum() / pna.sum()),
                 pn_apl_per_bouton_rel=float(pna[bm].mean() / pna.mean()) if bm.any() else np.nan)
        for t in TYPES:
            tm = typ == t
            kt = kc_type.reindex(per_kc.index).to_numpy() == t
            r[f'n_claws_{TYPE_COL[t]}'] = int((cm & tm).sum())
            r[f'share_claws_{TYPE_COL[t]}'] = float((cm & tm).sum() / tm.sum())
            r[f'kcs_2plus_{TYPE_COL[t]}'] = int(((per_kc >= 2).to_numpy() & kt).sum())
            r[f'share_kcs_2plus_{TYPE_COL[t]}'] = float(((per_kc >= 2).to_numpy() & kt).sum() / kt.sum())
        rows.append(r)
        for gl, inp in g.items():
            cg = has & (claw_glom == gl)
            bg = bglom == gl
            rr = dict(calyx=c, odor=od, glomerulus=gl, orn_input=float(inp), n_claws=int(cg.sum()), n_boutons=int(bg.sum()),
                      apl_kc=float(akc[cg].sum()), kc_apl=float(kca[cg].sum()), apl_pn=float(apn[bg].sum()),
                      pn_apl=float(pna[bg].sum()))
            for t in TYPES:
                rr[f'n_claws_{TYPE_COL[t]}'] = int((cg & (typ == t)).sum())
            by.append(rr)
    return rows, by


def _f(x, sp: str = '.3f') -> str:
    try:
        return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(float(x), sp)
    except (TypeError, ValueError):
        return str(x)


def step0_table(rd: pd.DataFrame, pc: pd.DataFrame, summ: pd.DataFrame, ov: pd.DataFrame, gl: pd.DataFrame) -> str:
    """Step 0's tables as markdown."""
    cal = rd[rd['readout'] == 'cal'].set_index('unit')
    L = ['# E133 step 0 (generated)', '', '## The pattern part by pair group (calibrated readout, all recordings)', '',
         '| run | 120 pairs | amplitude | shape | spread ratio | Pearson r | 42 pairs with 2-but, paa or B-myr | the other 78 '
         '| 39 (exactly one) | 81 | 3 among the three | halves A / B / C / D, 42 | halves, 78 |',
         '|---|' + '---:|' * 12]
    units = [f'S1:{c}' for c in CALYCES] + [f'fit:{c}' for c in CALYCES] + [u for u in cal.index if u.startswith('family:')]
    for u in units:
        r = cal.loc[u]
        L.append(f"| {u} | {_f(r['all_pattern'], '.4f')} | {_f(r['all_amplitude'], '.4f')} | {_f(r['all_shape'], '.4f')} | "
                 f"{_f(r['all_sd_ratio'], '.2f')} | {_f(r['all_pearson_effect'])} | {_f(r['all_pattern_p42'])} | "
                 f"{_f(r['all_pattern_p78'])} | {_f(r['all_pattern_p39'])} | {_f(r['all_pattern_p81'])} | "
                 f"{_f(r['all_pattern_p3'])} | " + ' / '.join(_f(r[f'{h}_pattern_p42']) for h in TARGETS_READ[1:]) + ' | '
                 + ' / '.join(_f(r[f'{h}_pattern_p78']) for h in TARGETS_READ[1:]) + ' |')
    L += ['', '## The planner\'s numbers', '', '| run | cell | planner | step 0 | reproduced |', '|---|---|---:|---:|---|']
    for r in pc.itertuples():
        L.append(f'| {r.unit} | {r.cell} | {r.planner} | {_f(r.step0, ".4f") if isinstance(r.step0, float) else r.step0} | '
                 f"{'yes' if r.reproduced else 'no'} |")
    L += ['', '## The 3 pairs among the three odors: per-pair APL effect (calibrated)', '']
    p3 = [c for c in cal.columns if c.startswith('p3_effect_model_')]
    L += ['| run | ' + ' | '.join(c.replace('p3_effect_model_', '') for c in p3) + ' |', '|---|' + '---:|' * len(p3)]
    for u in units[:6]:
        L.append(f'| {u} | ' + ' | '.join(_f(cal.loc[u, c], '+.3f') for c in p3) + ' |')
    L.append('| the flies (all recordings) | ' + ' | '.join(_f(cal.iloc[0][c.replace('p3_effect_model_', 'all_p3_effect_data_')], '+.3f')
                                                       for c in p3) + ' |')
    L += ['', '## The three odors\' responders (rate, two or more spikes; all KCs and by type)', '',
          '| run | condition | odor | all | g | a\'b\' | ab | claws per responder, g / a\'b\' / ab |', '|---|---|---|---:|---:|---:|---:|---|']
    s = summ.set_index(['unit', 'condition', 'odor', 'kc_type'])
    for u in units[:6]:
        for cond in ('on', 'off'):
            for od in THREE:
                L.append(f'| {u} | {cond} | {od} | ' + ' | '.join(_f(s.loc[(u, cond, od, t), 'rate'], '.4f') for t in ('all',) + TYPES)
                         + ' | ' + ' / '.join(_f(s.loc[(u, cond, od, t), 'claws_per_kc_responders'], '.2f') for t in TYPES) + ' |')
    L += ['', '## Overlap of the responder sets (KCs responding to both; Jaccard)', '',
          '| run | condition | ' + ' | '.join(f'{a} & {b}' for i, a in enumerate(THREE) for b in THREE[i + 1:]) + ' | all three |',
          '|---|---|' + '---:|' * 4]
    o2 = ov.set_index(['unit', 'condition', 'odors'])
    for u in units[:6]:
        for cond in ('on', 'off'):
            cells = [f"{int(o2.loc[(u, cond, f'{a} & {b}'), 'n_both'])} ({_f(o2.loc[(u, cond, f'{a} & {b}'), 'jaccard'], '.2f')})"
                     for i, a in enumerate(THREE) for b in THREE[i + 1:]]
            cells.append(f"{int(o2.loc[(u, cond, 'all three'), 'n_both'])} ({_f(o2.loc[(u, cond, 'all three'), 'jaccard'], '.2f')})")
            L.append(f'| {u} | {cond} | ' + ' | '.join(cells) + ' |')
    L += ['', '## The three odors\' glomeruli in each calyx', '',
          '| calyx | odor | glomeruli | claws (share) | claws g / a\'b\' / ab | KCs with 2+ such claws g / a\'b\' / ab | boutons | '
          'APL->claw share (per claw x mean) | claw->APL share | APL->bouton share (per bouton x mean) | bouton->APL share |',
          '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in gl.itertuples():
        L.append(f'| {SHORT[r.calyx]} | {r.odor} | {r.glomeruli} | {r.n_claws} ({_f(r.share_claws)}) | '
                 f'{r.n_claws_g} / {r.n_claws_apbp} / {r.n_claws_ab} | {r.kcs_2plus_g} / {r.kcs_2plus_apbp} / {r.kcs_2plus_ab} | '
                 f'{r.n_boutons} | {_f(r.share_apl_kc)} ({_f(r.apl_kc_per_claw_rel, ".2f")}) | {_f(r.share_kc_apl)} | '
                 f'{_f(r.share_apl_pn)} ({_f(r.apl_pn_per_bouton_rel, ".2f")}) | {_f(r.share_pn_apl)} |')
    return '\n'.join(L) + '\n'


def part_step0() -> dict:
    """Step 0 (the module docstring). Writes step0_patterns.csv, step0_planner.csv, step0_responders*.csv,
    step0_overlap.csv, step0_glomeruli*.csv, step0_table.md and step0.json."""
    out = out_dir()
    units = step0_units()
    rd = step0_patterns(units)
    rd.to_csv(out / 'step0_patterns.csv', index=False)
    pc = planner_check(rd)
    pc.to_csv(out / 'step0_planner.csv', index=False)
    det, summ, ov = responders(units)
    det.to_csv(out / 'step0_responders_by_claws.csv', index=False)
    summ.to_csv(out / 'step0_responders.csv', index=False)
    ov.to_csv(out / 'step0_overlap.csv', index=False)
    gloms = odor_glomeruli()
    gl, by = [], []
    for c in CALYCES:
        a, b = glomerulus_tables(c, gloms)
        gl += a
        by += b
    gl, by = pd.DataFrame(gl), pd.DataFrame(by)
    gl.to_csv(out / 'step0_glomeruli.csv', index=False)
    by.to_csv(out / 'step0_glomeruli_by_glomerulus.csv', index=False)
    (out / 'step0_table.md').write_text(step0_table(rd, pc, summ, ov, gl))
    cal = rd[rd['readout'] == 'cal'].set_index('unit')
    holds = {}
    for t in TARGETS_READ:
        holds[t] = dict(left_minus_hemi_p42=float(cal.loc[f'S1:{LEFT}', f'{t}_pattern_p42'] - cal.loc[f'S1:{HEMI}', f'{t}_pattern_p42']),
                        left_minus_hemi_p78=float(cal.loc[f'S1:{LEFT}', f'{t}_pattern_p78'] - cal.loc[f'S1:{HEMI}', f'{t}_pattern_p78']),
                        left_minus_hemi_p39=float(cal.loc[f'S1:{LEFT}', f'{t}_pattern_p39'] - cal.loc[f'S1:{HEMI}', f'{t}_pattern_p39']),
                        left_minus_hemi_p81=float(cal.loc[f'S1:{LEFT}', f'{t}_pattern_p81'] - cal.loc[f'S1:{HEMI}', f'{t}_pattern_p81']))
    res = dict(generated=_now(), host=_host(), script=SCRIPT, commits=H.code_commits(),
               planner_reproduced_42_78=bool(pc[pc['column'].isin(['all_pattern', 'all_pattern_p42', 'all_pattern_p78'])]['reproduced'].all()),
               planner_reproduced_39_81=bool(pc[pc['column'].isin(['all_pattern', 'all_pattern_p39', 'all_pattern_p81'])]['reproduced'].all()),
               n_pairs={g: int(m.sum()) for g, m in pair_groups(_canonical_pairs()).items()},
               fixed_scales_left_minus_hemi=holds, odor_glomeruli={od: list(g.index) for od, g in gloms.items()},
               note=("The planner's '39' and '81' are the 42 pairs that hold any of 2-but, paa and B-myr and the 78 that hold "
                     "none (step0_planner.csv); p39 / p81 are the counts as written (exactly one of the three; the 3 pairs "
                     "among them counted with the others)."))
    _write_json(out / 'step0.json', res)
    print((out / 'step0_table.md').read_text(), flush=True)
    print(json.dumps({k: res[k] for k in ('planner_reproduced_42_78', 'planner_reproduced_39_81', 'n_pairs')}), flush=True)
    return res


def _canonical_pairs() -> pd.Index:
    """The 120 pairs of the 16 odors in E83's canonical order (sorted names)."""
    from experiments import E83_cross_wiring_fit as E83
    od = list(H.APL16_ODORS)
    s = pd.Series(0.0, index=pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(od) for b in od[i + 1:]],
                                                       names=['odor_a', 'odor_b']))
    return E83._cs(s).index


# ---------------------------------------------------------------------------
# the tables, the simulator and the hold

def base_dir() -> Path:
    """This host's cache of the bases' APL-off breadth runs (E10_joint_fit.bouton_base's)."""
    return local_dir('base')


def calyx_tables(c: str, seed: int = 0) -> dict:
    """A calyx's tables at one of its wirings (seed 0 its own; n the within-type claw permutation of seed n), cached in
    the process: 'calyx', 'seed', 'label', 'orn', 'base' (E130's bases.json), 'wnorm' (normalized APL weights: per claw
    'wAPLKC', 'wKCAPL'; per bouton 'wAPLPN', 'wPNAPL'), 'W' (claw x bouton), 'breadth' (per-KC APL-off total spike count
    on that wiring), 'kckc' (the male calyx's matrix, E130's; None for the hemibrain, whose matrix olfsysm loads),
    'meta'. Built by E10_joint_fit.bouton_base with the connectome, E130's base and, for seed n, `wiring=
    'claw_permutation_within_type', wiring_seed=n` (the APL-off breadth run cached under outputs/.../base/)."""
    from experiments import E10_joint_fit as J
    key = (c, int(seed))
    if key in _CACHE.setdefault('tables', {}):
        return _CACHE['tables'][key]
    if c not in CALYCES:
        raise ValueError(f'{c=} not in {CALYCES}')
    E130 = e130()
    orn = H.apl16_orn_deltas()
    wk = dict(wiring=WIRING, wiring_seed=int(seed)) if seed else {}
    initial, wnorm, W, breadth = J.bouton_base(orn, base_dir(), None, connectome=c, base_params=E130.base_params_spec(c), **wk)
    kckc, ksha = male_matrix(c, W) if c in MALE else (None, None)
    T = dict(calyx=c, seed=int(seed), label=f'{c}_seed{int(seed)}', orn=orn, base=dict(initial), wnorm=wnorm, W=W,
             breadth=breadth, kckc=kckc, meta=dict(wiring=WIRING if seed else 'own', wiring_seed=int(seed),
                                                   W_sha256=sha256_frame(W), breadth_sha256=sha256_frame(breadth.to_frame()),
                                                   kckc_sha256=ksha))
    _CACHE['tables'][key] = T
    return T


def male_matrix(c: str, W: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    """A male calyx's KC x KC matrix among the wiring's KCs (E130's `malecns_kckc`; a claw permutation keeps every KC,
    so it is the calyx's own), checked against the sha256 in E130's matrices.json (no file written); (matrix, sha256)."""
    from experiments.E77_input_screen import frame_sha256
    E130 = e130()
    key = ('kckc', c)
    if key not in _CACHE:
        M = E130.malecns_kckc(c, E130.model_kc_ids(W))
        sha = frame_sha256(M)
        want = _read_json(E130.out_dir() / 'matrices.json').get(c, {}).get('sha256')
        if sha != want:
            raise RuntimeError(f"{c}: the KC-KC matrix differs from E130's matrices.json ({sha} vs {want})")
        _CACHE[key] = (M, sha)
    M, sha = _CACHE[key]
    if not M.index.equals(E130.model_kc_ids(W).rename(M.index.name)):
        raise RuntimeError(f"{c}: the wiring's KCs are not the matrix's")
    return M, sha


def drop_tables(c: str, seed: int) -> None:
    _CACHE.setdefault('tables', {}).pop((c, int(seed)), None)
    gc.collect()


def breadth_run(tables: Mapping) -> pd.Series:
    """The per-KC APL-off breadth of `tables`: E10_joint_fit.bouton_base's APL-off run (the base's threshold, APL->KC
    and APL->bouton 0, the base's KC->APL and bouton->APL, no KC-KC) on the tables' wiring passed as `_wPNKC`; each
    KC's total spike count, indexed by KC ID. The path step 2 takes for a changed wiring."""
    from experiments import E10_joint_fit as J
    from pyolfsysm.mb_model import fit_mb_model
    c, b = tables['calyx'], tables['base']
    kws = dict(J.base_kws(None, connectome=c), _wPNKC=tables['W'])
    _, off_spk, _, _ = fit_mb_model(tables['orn'], fixed_thr=b['fixed_thr'], wAPLKC=0.0, wKCAPL=b['wKCAPL_scale'],
                                    wAPLPN=0.0, wPNAPL=b['wPNAPL_scale'], **kws)
    return pd.Series(off_spk.astype(float).sum(axis='columns').to_numpy(), index=off_spk.index.get_level_values(H.KC_ID_LEVEL))


def build_simulator(tables: Mapping, params: Mapping[str, float]) -> Tuple:
    """C2's family on `tables` as a fixed-parameter simulator at `params`: E130's `build_simulator` on given tables
    (the connectome's kwargs, the wiring, the breadth-weighted APL->KC weights (exponent 2), the other normalized APL
    weights, the male calyx's KC-KC matrix as `_wKCKC`); (simulate, info), info's 'kws' the scalar model kwargs and
    'breadth' the breadth behind the APL->KC weights."""
    from experiments import E10_joint_fit as J
    from experiments.E77_input_screen import _scalar_kws
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator
    c = tables['calyx']
    p = dict(params)
    wnorm = tables['wnorm']
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], tables['breadth'], K)
    extra = {} if tables.get('kckc') is None else {'_wKCKC': tables['kckc']}
    kws = J.base_kws(None, connectome=c)
    sim, _ = make_fixed_param_kc_simulator(tables['orn'], initial_params=p, **kws, **extra, _wPNKC=tables['W'],
                                           _wAPLKC=wAPLKC, _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'],
                                           _wPNAPL=wnorm['wPNAPL'], kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)
    return sim, dict(params=p, kws=_scalar_kws(dict(kws, kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)), connectome=c,
                     multibouton_claws=kws['multibouton_claws'], breadth=tables['breadth'], wAPLKC=wAPLKC, W=tables['W'],
                     orn=tables['orn'])


def start_and_targets(tables: Mapping) -> Tuple[Dict[str, float], int, int]:
    """E130's step-1 start on a calyx's tables (C2's five factors on the calyx's base, KC-KC at C2's scale;
    E130.transfer_params) and its held counts (C2's counts as fractions of its KCs x odors, times the calyx's;
    E130.held_targets)."""
    from experiments import E83_cross_wiring_fit as E83
    E130 = e130()
    F = E83.saved_fit('real')
    p = E130.transfer_params(tables['base'], True)
    n_kcs = int(tables['breadth'].shape[0])
    on, off = E130.held_targets(n_kcs, tables['orn'].shape[1], F['target_on'], F['target_off'], F['size'])
    return p, on, off


def hold(tables: Mapping, label: str) -> dict:
    """E130's step 1 on `tables`: from `start_and_targets`, the threshold bisected with APL off, then the tied APL
    scales with APL on (E77_input_screen.hold_rates); its result with 'start'."""
    from experiments.E77_input_screen import hold_rates
    p, on, off = start_and_targets(tables)
    res = hold_rates(lambda: build_simulator(tables, p), target_off=off, target_on=on, params=dict(p), label=label)
    res['start'] = p
    return res


def _ndiff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> Dict[str, object]:
    return dict(counts=int((a[1].to_numpy(dtype=float) != b[1].to_numpy(dtype=float)).sum()),
                responses=int((a[0].to_numpy(dtype=float) != b[0].to_numpy(dtype=float)).sum()),
                same_index=bool(a[1].index.equals(b[1].index) and a[1].columns.equals(b[1].columns)))


def _diff_ok(d: Mapping) -> bool:
    return bool(d['counts'] == 0 and d['responses'] == 0 and d['same_index'])


def e130_held_params(c: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    """E130's S1_<calyx> held parameters, APL on and off (its runs table's params_on / params_off)."""
    r = e130_runs_table().loc[f'S1_{c}']
    return ({k: float(v) for k, v in json.loads(r['params_on']).items()},
            {k: float(v) for k, v in json.loads(r['params_off']).items()})


def replay(c: str) -> dict:
    """Identity (ii): the runner on the calyx's own tables at E130's held parameters against E130's S1_<calyx>, APL on
    and off; with the tables' agreement with E130's wiring (E130's `step0_wiring_only`) and breadth (its archived APL-off
    base run)."""
    E130 = e130()
    t0 = time.time()
    T = calyx_tables(c, 0)
    ref = e130_s1(c)
    p_on, p_off = e130_held_params(c)
    sim, _ = build_simulator(T, p_on)
    diff = {}
    for cnd, q in (('on', p_on), ('off', p_off)):
        r, s = sim(dict(q))
        diff[cnd] = _ndiff((r, s), ref[cnd])
    del sim
    gc.collect()
    W0 = E130.step0_wiring_only(c)
    arch = E130.archive_dir() / 'base' / f'cx-{c}_base_acute_off_spike_counts.parquet'
    b_ok = None
    if arch.exists():
        a = pd.read_parquet(arch)
        b_ref = pd.Series(a.astype(float).sum(axis='columns').to_numpy(), index=a.index.get_level_values(H.KC_ID_LEVEL))
        b_ok = bool(b_ref.equals(T['breadth']))
    rec = dict(calyx=c, differences=diff, wiring_same_as_e130=bool(T['W'].equals(W0)), breadth_same_as_e130=b_ok,
               params_on=p_on, params_off=p_off, seconds=time.time() - t0, host=_host(), checked=_now(),
               olfsysm=olfsysm_build(), commits=H.code_commits())
    rec['passes'] = bool(all(_diff_ok(d) for d in diff.values()) and rec['wiring_same_as_e130'] and b_ok is not False)
    return rec


def olfsysm_build() -> Dict[str, object]:
    import olfsysm
    bi = olfsysm.build_info()
    return dict(olfsysm_git_hash=bi.get('git_hash'), olfsysm_git_dirty=bi.get('git_dirty'),
                olfsysm_cpp_dirty=bi.get('cpp_dirty'), olfsysm_source_sha256=bi.get('source_sha256'))


# ---------------------------------------------------------------------------
# the check

def part_check(calyces: Sequence[str] = CALYCES) -> dict:
    """The identity (the module docstring's part check), each calyx kept once passed; check_<host>.json. Raises if one
    fails."""
    path = out_dir() / f'check_{_host()}.json'
    res = _read_json(path)
    res.setdefault('calyces', {})
    for c in calyces:
        cur = res['calyces'].get(c, {})
        if cur.get('passes'):
            print(f'E133 check {c}: passed before', flush=True)
            continue
        t0 = time.time()
        rec = dict(replay=replay(c))
        print(f"E133 check {c} (ii) replay: {rec['replay']['differences']}, wiring {rec['replay']['wiring_same_as_e130']}, "
              f"breadth {rec['replay']['breadth_same_as_e130']}", flush=True)
        T = calyx_tables(c, 0)
        b = breadth_run(T)
        rec['breadth_run_same'] = bool(b.equals(T['breadth']))
        print(f"E133 check {c} (i) breadth_run: {rec['breadth_run_same']}", flush=True)
        r = e130_runs_table().loc[f'S1_{c}']
        h = hold(T, f'E133 check {c}')
        ref = e130_s1(c)
        p_on, _ = e130_held_params(c)
        # the held factors and parameters are compared to a relative 1e-12: the bisection's geometric steps can differ
        # in the last bit between CPUs (E130's S1 runs were held on atlas; the counts are compared exactly)
        rel = {f'factor_{k}': _rel(h[f'{k}_factor'], float(r[f'{k}_factor'])) for k in ('thr', 'apl')}
        rel.update({k: _rel(h['params'][k], p_on[k]) for k in PARAM_KEYS})
        rec['hold'] = dict(thr_factor=h['thr_factor'], apl_factor=h['apl_factor'], e130_thr_factor=float(r['thr_factor']),
                           e130_apl_factor=float(r['apl_factor']), n_sims=h['n_sims'], seconds=h.get('seconds'),
                           start_same=bool(all(h['start'][k] == json.loads(r['start_params'])[k] for k in PARAM_KEYS)),
                           params_bit_identical=bool(all(float(h['params'][k]) == p_on[k] for k in PARAM_KEYS)),
                           rel_diff=rel, max_rel_diff=float(max(rel.values())), rel_tol=REL_TOL,
                           differences={cnd: _ndiff(h[cnd], ref[cnd]) for cnd in ('on', 'off')})
        rec['hold']['passes'] = bool(rec['hold']['start_same'] and rec['hold']['max_rel_diff'] <= REL_TOL
                                     and all(_diff_ok(d) for d in rec['hold']['differences'].values()))
        print(f"E133 check {c} (iii) hold: {rec['hold']}", flush=True)
        rec.update(passes=bool(rec['replay']['passes'] and rec['breadth_run_same'] and rec['hold']['passes']),
                   seconds=time.time() - t0, checked=_now(), host=_host(), olfsysm=olfsysm_build(), commits=H.code_commits())
        res['calyces'][c] = rec
        _write_json(path, res)
        del h
        gc.collect()
    res['passes'] = bool(all(res['calyces'].get(c, {}).get('passes') for c in calyces))
    _write_json(path, res)
    print('E133 check: ' + ', '.join(f"{c} {'passes' if res['calyces'][c]['passes'] else 'FAILS'}" for c in calyces),
          flush=True)
    if not res['passes']:
        raise RuntimeError(f'E133 check failed ({path})')
    return res


def check_passed(c: str) -> bool:
    """Whether any host's check_<host>.json has passed the full identity on calyx `c`."""
    return any(_read_json(p).get('calyces', {}).get(c, {}).get('passes') for p in out_dir().glob('check_*.json'))


def replay_on_host(calyces: Sequence[str]) -> dict:
    """The replay (identity (ii)) of each calyx on this host, once (replay_<host>.json); RuntimeError on a
    difference."""
    path = out_dir() / f'replay_{_host()}.json'
    res = _read_json(path)
    chk = _read_json(out_dir() / f'check_{_host()}.json').get('calyces', {})
    for c in calyces:
        if res.get(c, {}).get('passes') or chk.get(c, {}).get('replay', {}).get('passes'):
            continue
        res[c] = replay(c)
        _write_json(path, res)
        print(f"E133 replay {c} on {_host()}: {res[c]['differences']} ({res[c]['seconds']:.0f}s)", flush=True)
        if not res[c]['passes']:
            raise RuntimeError(f"E133: the runner does not reproduce E130's S1_{c} on {_host()} ({path})")
    return res


# ---------------------------------------------------------------------------
# step 1: the points

def planned_step1(calyces: Sequence[str] = CALYCES, seeds: Sequence[int] = SEEDS1, reverse: bool = False) -> List[Point]:
    """Step 1's points of a lane: first the left's and the hemibrain's (the 14 runs the class reads), seed by seed in
    the order of `seeds`, each seed's calyces in the order given; then the right's (read beside, no class). `reverse`
    takes the seeds from the other end in both blocks, so a lane run with it meets a forward lane in the middle of each
    block."""
    if any(s not in SEEDS1 for s in seeds) or any(c not in CALYCES for c in calyces):
        raise ValueError(f'calyces {calyces}, seeds {seeds}: outside {CALYCES} x {SEEDS1}')
    ss = [int(s) for s in seeds][::-1] if reverse else [int(s) for s in seeds]
    first = [Point(1, c, s, None) for s in ss for c in calyces if c != RIGHT]
    return first + ([Point(1, RIGHT, s, None) for s in ss] if RIGHT in calyces else [])


def runs_table(step: int = 1) -> pd.DataFrame:
    """Every row of every host and lane of a step (runs_<step>_*.csv), the last row per tag."""
    fs = [pd.read_csv(p) for p in sorted(out_dir().glob(f'runs_{step}_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag', 'calyx', 'seed'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def claim(tag: str, lane: str, reclaim: bool = False, root: Optional[Path] = None) -> bool:
    """Claims a point by mkdir of results/<EXPERIMENT>/claims/<tag>/ (atomic on NFS); False if another process holds
    it (unless `reclaim`)."""
    d = (out_dir() if root is None else root) / 'claims' / tag
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(tag=tag, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(tag: str, root: Optional[Path] = None) -> None:
    shutil.rmtree((out_dir() if root is None else root) / 'claims' / tag, ignore_errors=True)


def save_run(tag: str, o: Mapping[str, tuple], part: str = 'step1', archive: bool = True) -> None:
    """The run's conditions, locally and (with `archive`) in the model-output archive, written atomically."""
    dirs = [local_dir(part)] + ([archive_dir() / part] if archive else [])
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1])):
                p = d / f'{tag}_{cc}_{k}.parquet'
                tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
                x.to_parquet(tmp)
                os.replace(tmp, p)


def load_run(tag: str, part: str = 'step1') -> Optional[Dict[str, tuple]]:
    from experiments import E83_cross_wiring_fit as E83
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, archive_dir() / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def pair_rows(tag: str, o: Mapping[str, tuple]) -> pd.DataFrame:
    """The run's raw-count pair distances (H.PRIMARY_MODEL_RDM, as E130's variants), long: tag, condition ('on' /
    'acute_off'), odor_a, odor_b, distance."""
    return pd.concat({(tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                      (tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                     names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index()


def save_pairs(tag: str, o: Mapping[str, tuple]) -> Path:
    p = out_dir() / 'pairs' / f'{tag}.csv'
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f'.{p.name}.{_host()}.{os.getpid()}.tmp')
    pair_rows(tag, o).to_csv(tmp, index=False)
    os.replace(tmp, p)
    return p


def corner_flags(res: Mapping, base: Mapping[str, float]) -> dict:
    """A held point at a corner: a bisection at its bracket's bound (hold_rates' at_bound_thr / at_bound_apl), or
    E124's bounds as E130's `bound_flags` reads them (the threshold within [1/4, 4] x the base's, APL->KC within
    [1/30, 30] x the start's); APL->bouton at its cap is reported only."""
    bf = e130().bound_flags(res['params'], res['start'], base)
    corner = bool(res['at_bound_thr'] or res['at_bound_apl'] or bf['thr_at_bound'] or bf['aplkc_at_bound'])
    return dict(bf, at_bound_thr=bool(res['at_bound_thr']), at_bound_apl=bool(res['at_bound_apl']), corner=corner)


def guard_flags(max_on: int, max_off: int, ref: Tuple[int, int], factor: float = MAX_COUNT_FACTOR) -> dict:
    """The entry's guard: the largest spike count per KC and odor above `factor` times the reference's (the calyx's own
    unpermuted run in step 1; the host wiring's own run in step 2), in either condition."""
    return dict(ref_max_count_on=int(ref[0]), ref_max_count_off=int(ref[1]),
                guard_max_count=bool(max_on > factor * ref[0] or max_off > factor * ref[1]))


def realism_row(o: Mapping[str, tuple]) -> dict:
    """E83's `realism` (rates, ratio, response rate ρ, odors per KC, largest counts, flag), E130's `type_rates`, E110's
    `run_extras` (silent KCs, KCs with 20 or more spikes to some odor), the five bands (E110's `bands_of`) and
    'in_bands'."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E110_no_kckc as E110
    r = dict(E83.realism(dict(o)))
    r.update(e130().type_rates(dict(o)))
    r.update(E110.run_extras(dict(o)))
    r.update(E110.bands_of(float(r['rate_on']), float(r['rate_off']), float(r['rate_rho_on']), float(r['opk_on_tv']),
                           float(r['opk_off_tv'])))
    r['in_bands'] = bool(all(r[b] for b in BAND_COLS))
    return r


def point_tables(pt: Point) -> Tuple[dict, dict]:
    """(the point's tables, what was changed): step 1, the calyx's tables at the seed; step 2, the arm's builder
    (experiments/E133_builders.py, imported here only) on the host's tables at the seed with the donor's own, the
    breadth recomputed (`breadth_run`) when the wiring changed."""
    if pt.step == 1:
        return calyx_tables(pt.calyx, pt.seed), dict(change='none')
    from experiments import E133_builders as B
    host, donor = calyx_tables(pt.calyx, pt.seed), calyx_tables(DONOR[pt.calyx], 0)
    T = dict(B.build(pt.arm, host, donor, seed=pt.seed))
    if not T['W'].equals(host['W']):
        T['breadth'] = breadth_run(T)
    T['label'] = tag_of(pt)
    return T, dict(change=f'arm {pt.arm} of {DONOR[pt.calyx]} on {pt.calyx} seed {pt.seed}', meta=T.get('meta', {}))


def reference_of(pt: Point) -> Tuple[int, int]:
    """The guard's reference (largest counts, APL on / off): E130's S1_<calyx> in step 1; in step 2 the host wiring's
    own run (E130's for seed 0, step 1's otherwise)."""
    if pt.step == 1 or pt.seed == 0:
        return reference_max_counts(pt.calyx)
    o = load_run(tag_of(Point(1, pt.calyx, pt.seed, None)))
    if o is None:
        raise RuntimeError(f'{tag_of(pt)}: no step-1 run of its host wiring')
    return int(o['on'][1].to_numpy(dtype=float).max()), int(o['off'][1].to_numpy(dtype=float).max())


def run_point(pt: Point, lane: str, table: Path, prov_path: Path, cc: dict, reclaim: bool = False,
              root: Optional[Path] = None, archive: bool = True, part: Optional[str] = None) -> Optional[dict]:
    """One point, claimed (None if done or claimed elsewhere): E130's step-1 hold on the point's tables, the run saved
    (step<n>/, locally and in the archive), its pair distances (pairs/<tag>.csv), its row appended to `table` and its
    provenance to `prov_path`. `root` (default the results folder) holds the claims and pairs (a smoke point elsewhere
    passes a scratch folder with `archive` False). Returns the row."""
    from experiments import E108_paa_extrapolation as E108
    tag = tag_of(pt)
    part = part or f'step{pt.step}'
    root = out_dir() if root is None else root
    done = set(pd.read_csv(table)['tag']) if table.exists() else set()
    if root == out_dir():
        done |= set(runs_table(pt.step)['tag'])
    if tag in done:
        print(f'E133 {tag}: done', flush=True)
        return None
    if not claim(tag, lane, reclaim, root):
        print(f'E133 {tag}: claimed by another process, skipped', flush=True)
        return None
    try:
        E108._still_clean(cc)
        t0 = time.time()
        T, change = point_tables(pt)
        res = hold(T, f'E133 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        if archive:
            save_run(tag, o, part)
        else:
            d = root / part
            d.mkdir(parents=True, exist_ok=True)
            for c, ccn in (('on', 'on'), ('off', 'acute_off')):
                o[c][0].to_parquet(d / f'{tag}_{ccn}_responses.parquet')
                o[c][1].to_parquet(d / f'{tag}_{ccn}_spike_counts.parquet')
        pp = root / 'pairs' / f'{tag}.csv'
        pp.parent.mkdir(parents=True, exist_ok=True)
        pair_rows(tag, o).to_csv(pp, index=False)
        rr = realism_row(o)
        row = dict(part=part, tag=tag, step=pt.step, calyx=pt.calyx, seed=pt.seed, arm=pt.arm or '', wiring=T['meta'].get('wiring', ''),
                   connectome=pt.calyx, multibouton_claws='argmax', start_params=json.dumps(res['start']),
                   thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
                   target_on=res['target_on'], target_off=res['target_off'], n_sims=res['n_sims'],
                   seconds_hold=res.get('seconds'), seconds=time.time() - t0,
                   **{k: float(res['params'][k]) for k in PARAM_KEYS},
                   params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                   params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}),
                   base_params=json.dumps(T['base']), W_sha256=sha256_frame(T['W']),
                   breadth_sha256=res.get('breadth_sha256'), host=_host(), lane=lane)
        row.update(rr)
        row.update(corner_flags(res, T['base']))
        row.update(guard_flags(rr['max_count_on'], rr['max_count_off'], reference_of(pt)))
        row['no_class'] = bool(not row['in_bands'] or row['corner'] or row['guard_max_count'])
        row.update(olfsysm_build())
        row.update(cc)
        from experiments import E83_cross_wiring_fit as E83
        E83._append_row(table, row)
        prov = _read_json(prov_path) or dict(part=part, lane=lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=T['orn']), held=dict(
            calyx=pt.calyx, seed=pt.seed, arm=pt.arm, change=change, wiring=T['meta'],
            base=T['base'], base_source=e130().base_params_spec(pt.calyx), start=res['start'],
            targets=dict(on=res['target_on'], off=res['target_off']),
            hold=("C2's five scales as factors on the calyx's own base, KC-KC at C2's scale; fixed_thr re-held with APL "
                  "off, then the tied APL scales with APL on (E77_input_screen.hold_rates), as E130's step 1"),
            kckc_matrix=(dict(source=f"E130's malecns_kckc (results/{E130_EXP}/matrices.json)", sha256=T['meta'].get('kckc_sha256'))
                         if pt.calyx in MALE else "olfsysm's hemibrain KC->KC counts (connectome_KCKC_weights)"),
            kws=res.get('kws'), breadth_sha256=res.get('breadth_sha256'), W_sha256=row['W_sha256'], history=res.get('history')))
        _write_json(prov_path, prov)
        print(f"E133 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, "
              f"{res['n_off']}/{res['target_off']} off); bands {row['in_bands']}, corner {row['corner']}, guard "
              f"{row['guard_max_count']} (max {row['max_count_on']}/{row['max_count_off']} against "
              f"{row['ref_max_count_on']}/{row['ref_max_count_off']}); {row['seconds']:.0f}s", flush=True)
        del res, o
        gc.collect()
        return row
    except BaseException:
        release(tag, root)
        raise


PRELOAD = ('experiments.E10_joint_fit', 'experiments.E77_input_screen', 'experiments.E83_cross_wiring_fit',
           'experiments.E108_paa_extrapolation', 'experiments.E110_no_kckc', 'experiments.E122_feedback_neuron',
           'experiments.E130_malecns_family', 'experiments.E72_wiring_population', 'experiments.E12_breadth_apl',
           'experiments.E69_kc_firing_ceiling', 'experiments.E77_input_side_common_error', 'kcq.apl_effect',
           'pyolfsysm.mb_model', 'pyolfsysm.tuning', 'pyolfsysm.wpnkc_shuffles')


def _preload() -> None:
    """Imports every module the simulation parts use before E108's commit check (E133_builders only for step 2)."""
    import importlib
    for m in PRELOAD:
        importlib.import_module(m)


def part_step1(calyces: Sequence[str] = CALYCES, seeds: Sequence[int] = SEEDS1, lane: Optional[str] = None,
               reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 1 (the module docstring), resumable, one point claimed at a time; refuses unless the full identity has
    passed on some host for each of the lane's calyces and the replay passes on this host."""
    pts = planned_step1(calyces, seeds, reverse)
    missing = [c for c in calyces if not check_passed(c)]
    if missing:
        raise SystemExit(f'E133 step 1 refused: the identity check has not passed for {missing} (run --part check)')
    _preload()
    from experiments import E108_paa_extrapolation as E108
    cc = E108._commits_clean()
    replay_on_host(list(dict.fromkeys(calyces)))
    lane = lane or lane_name(calyces, seeds, reverse)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    for pt in pts:
        run_point(pt, lane, table, prov_path, cc, reclaim)
        drop_tables(pt.calyx, pt.seed)
    return runs_table(1)


def smoke_point(c: str = LEFT, seed: int = 1, scratch: Optional[Path] = None) -> dict:
    """One step-1 point into a scratch folder (claims, pairs, runs and its row there; nothing in results or the
    archive); returns the row."""
    _preload()
    from experiments import E108_paa_extrapolation as E108
    cc = E108._commits_clean()
    scratch = Path(scratch) if scratch is not None else local_dir('smoke')
    scratch.mkdir(parents=True, exist_ok=True)
    return run_point(Point(1, c, seed, None), 'smoke', scratch / 'runs_smoke.csv', scratch / 'provenance_smoke.json', cc,
                     root=scratch, archive=False, part='step1')


# ---------------------------------------------------------------------------
# step 1's score

def class_step1(left: Sequence[float], hemi: Sequence[float], all_classed: bool = True) -> dict:
    """Step 1's class from the left's and the hemibrain's re-draws' pattern parts (the entry's rule): L and H their
    means; "kept" if L - H >= KEPT_GAP and at least KEPT_ABOVE of the left's are above the hemibrain's median,
    "re-dealt" if L - H <= REDEALT_GAP, otherwise "mixed"; no class (cls None, 'would_be' the class) unless
    `all_classed`; 'incomplete' without seven of each."""
    left, hemi = np.asarray(left, dtype=float), np.asarray(hemi, dtype=float)
    n_l, n_h = int(np.isfinite(left).sum()), int(np.isfinite(hemi).sum())
    if n_l < len(SEEDS1) or n_h < len(SEEDS1):
        return dict(cls='incomplete', n_left=n_l, n_hemi=n_h)
    L, Hm, med = float(left.mean()), float(hemi.mean()), float(np.median(hemi))
    above = int((left > med).sum())
    gap = L - Hm
    if gap >= KEPT_GAP and above >= KEPT_ABOVE:
        cls = 'kept'
    elif gap <= REDEALT_GAP:
        cls = 're-dealt'
    else:
        cls = 'mixed'
    det = dict(L=L, H=Hm, gap=gap, hemi_median=med, left_above_hemi_median=above, n_left=n_l, n_hemi=n_h)
    if not all_classed:
        return dict(det, cls=None, would_be=cls, reason='a run without a class (bands, corner or guard)')
    return dict(det, cls=cls)


def beside_right(right: Sequence[float], hemi: Sequence[float], left: Sequence[float]) -> dict:
    """The right calyx's seven beside, no class: its mean R, R - H, R - L, and how many are above the hemibrain's
    median."""
    right, hemi, left = (np.asarray(x, dtype=float) for x in (right, hemi, left))
    if not (np.isfinite(right).sum() and np.isfinite(hemi).sum()):
        return dict(n_right=int(np.isfinite(right).sum()))
    med = float(np.nanmedian(hemi))
    return dict(R=float(np.nanmean(right)), R_minus_H=float(np.nanmean(right) - np.nanmean(hemi)),
                R_minus_L=float(np.nanmean(right) - np.nanmean(left)) if np.isfinite(left).sum() else np.nan,
                right_above_hemi_median=int((right > med).sum()), n_right=int(np.isfinite(right).sum()))


def collect1() -> Tuple[Dict[str, Dict[str, tuple]], pd.DataFrame]:
    """Step 1's runs and the calyces' own (E130's S1_<calyx>, as seed 0): ({tag: outputs}, runs table)."""
    outs, rows = {}, []
    for c in CALYCES:
        tag = tag_of(Point(1, c, 0, None))
        o = e130_s1(c)
        outs[tag] = o
        r = e130_runs_table().loc[f'S1_{c}']
        row = dict(tag=tag, step=1, calyx=c, seed=0, wiring='own', source=f'{E130_EXP} S1_{c}',
                   thr_factor=float(r['thr_factor']), apl_factor=float(r['apl_factor']),
                   **{k: float(r[k]) for k in PARAM_KEYS})
        row.update(realism_row(o))
        row.update(guard_flags(row['max_count_on'], row['max_count_off'], reference_max_counts(c)))
        row.update(corner_flags(dict(params=e130_held_params(c)[0], start=json.loads(r['start_params']),
                                     at_bound_thr=bool(r['at_bound_thr']), at_bound_apl=bool(r['at_bound_apl'])),
                                e130().calyx_base(c)))
        row['no_class'] = bool(not row['in_bands'] or row['corner'] or row['guard_max_count'])
        rows.append(row)
    for r in runs_table(1).to_dict('records'):
        o = load_run(r['tag'])
        if o is None:
            print(f"E133 score1: {r['tag']} has a row but no outputs here or in the archive", flush=True)
            continue
        outs[r['tag']] = o
        rows.append(dict(r, source=f"runs_1 ({r['host']} {r['lane']})"))
    return outs, pd.DataFrame(rows)


def readings1(outs: Mapping[str, Dict[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per run and readout: `unit_readings`, the calyx and seed."""
    idx = runs.set_index('tag')
    rows = []
    for ro in READOUTS:
        for tag, o in outs.items():
            rows.append(dict(tag=tag, calyx=idx.loc[tag, 'calyx'], seed=int(idx.loc[tag, 'seed']), readout=ro,
                             **unit_readings(_dist(o, ro))))
    _save_dist_cache()
    return pd.DataFrame(rows)


def verdict1(rd: pd.DataFrame, runs: pd.DataFrame) -> dict:
    """Step 1's class on all recordings (calibrated readout), the halves A to D beside, and the right calyx's seven
    beside; the grouped errors' means; the runs without a class."""
    cal = rd[rd['readout'] == 'cal']
    ri = runs.set_index('tag')
    out = dict(rules=dict(kept_gap=KEPT_GAP, redealt_gap=REDEALT_GAP, kept_above=KEPT_ABOVE, max_count_factor=MAX_COUNT_FACTOR,
                          readout='cal', column='<target>_pattern', targets=list(TARGETS_READ)))

    def vals(c, col):
        s = cal[(cal['calyx'] == c) & (cal['seed'] > 0)].set_index('seed')[col]
        return np.array([float(s.get(k, np.nan)) for k in SEEDS1])

    def classed(c):
        tags = [tag_of(Point(1, c, s, None)) for s in SEEDS1]
        return [t for t in tags if t in ri.index and bool(ri.loc[t, 'no_class'])]

    nc = classed(LEFT) + classed(HEMI)
    out['no_class_runs'] = nc
    for t in TARGETS_READ:
        col = f'{t}_pattern'
        out[t] = dict(class_step1(vals(LEFT, col), vals(HEMI, col), all_classed=not nc),
                      right=beside_right(vals(RIGHT, col), vals(HEMI, col), vals(LEFT, col)),
                      own={c: float(cal[(cal['calyx'] == c) & (cal['seed'] == 0)][col].iloc[0]) for c in CALYCES
                           if ((cal['calyx'] == c) & (cal['seed'] == 0)).any()},
                      groups={g: {c: _nanmean(vals(c, f'{t}_pattern_{g}')) for c in CALYCES} for g in GROUPS})
    out['class'] = out['all'].get('cls')
    out['realism_flags'] = [dict(tag=r['tag'], flag=r.get('realism_flag'), bands={b: bool(r[b]) for b in BAND_COLS},
                                 corner=bool(r.get('corner')), guard=bool(r.get('guard_max_count')))
                            for r in runs.to_dict('records') if bool(r.get('no_class'))]
    return out


STEP1_STATS = (('pattern part', 'all_pattern', 'cal', '.4f'), ('halves A / B / C / D', None, 'cal', '.3f'),
               ('amplitude part', 'all_amplitude', 'cal', '.4f'), ('shape part', 'all_shape', 'cal', '.4f'),
               ('spread ratio', 'all_sd_ratio', 'cal', '.2f'), ('Pearson r on the per-pair APL effect', 'all_pearson_effect', 'cal', '.3f'),
               ('ρ on the per-pair APL effect', 'all_rho_effect', 'cal', '.3f'), ('mean APL effect', 'mean_effect', 'cal', '+.3f'),
               ('error on the 42 pairs with 2-but, paa or B-myr', 'all_pattern_p42', 'cal', '.3f'),
               ('error on the other 78', 'all_pattern_p78', 'cal', '.3f'),
               ('error on the 39 (exactly one)', 'all_pattern_p39', 'cal', '.3f'), ('error on the 81', 'all_pattern_p81', 'cal', '.3f'),
               ('APL-effect loss, raw', 'all_loss', 'raw', '.4f'), ('mean APL effect, raw', 'mean_effect', 'raw', '+.3f'))
RUN_STATS = (('response rate, APL on / off', ('rate_on', 'rate_off'), '.4f'), ('APL-off / APL-on', ('ratio',), '.2f'),
             ('response rate ρ, APL on', ('rate_rho_on',), '.3f'), ('odors per KC (TV), on / off', ('opk_on_tv', 'opk_off_tv'), '.3f'),
             ('γ rate on / off', ('rate_on_g', 'rate_off_g'), '.3f'), ('γ off / on', ('ratio_g',), '.2f'),
             ("α′β′ rate on / off", ('rate_on_apbp', 'rate_off_apbp'), '.3f'), ('αβ rate on / off', ('rate_on_ab', 'rate_off_ab'), '.3f'),
             ('largest count, on / off', ('max_count_on', 'max_count_off'), '.0f'),
             ('KCs with 20+ spikes, on / off', ('n_kc_ge20_on', 'n_kc_ge20_off'), '.0f'),
             ('threshold', ('fixed_thr',), '.1f'), ('APL->KC scale', ('wAPLKC_scale',), '.4f'),
             ('threshold / APL factor of the hold', ('thr_factor', 'apl_factor'), '.3f'))


def step1_table(rd: pd.DataFrame, runs: pd.DataFrame, v: Mapping) -> str:
    cal = rd.set_index(['tag', 'readout'])
    ri = runs.set_index('tag')
    a = v['all']
    L = ['# E133 step 1 (generated)', '',
         f"Class: {a.get('cls')}{' (would be ' + str(a.get('would_be')) + ')' if a.get('would_be') else ''}; L - H "
         f"{_f(a.get('gap'), '+.4f')} (L {_f(a.get('L'), '.4f')}, H {_f(a.get('H'), '.4f')}); left above the hemibrain's median: "
         f"{a.get('left_above_hemi_median')} of 7; the right's seven beside: R - H {_f(a.get('right', {}).get('R_minus_H'), '+.4f')}.",
         '', '| target | class | L - H | left above median | R - H |', '|---|---|---:|---:|---:|']
    for t in TARGETS_READ:
        x = v[t]
        L.append(f"| {t} | {x.get('cls') or x.get('would_be')} | {_f(x.get('gap'), '+.4f')} | {x.get('left_above_hemi_median')} | "
                 f"{_f(x.get('right', {}).get('R_minus_H'), '+.4f')} |")
    for c in CALYCES:
        tags = [t for t in [tag_of(Point(1, c, s, None)) for s in (0,) + SEEDS1] if t in ri.index]
        L += ['', f'## {SHORT[c]}', '', '| statistic | ' + ' | '.join(f"seed {int(ri.loc[t, 'seed'])}" for t in tags) + ' |',
              '|---|' + '---:|' * len(tags)]
        for lab, col, ro, sp in STEP1_STATS:
            if col is None:
                cells = [' / '.join(_f(cal.loc[(t, ro), f'{h}_pattern'], sp) for h in TARGETS_READ[1:]) for t in tags]
            else:
                cells = [_f(cal.loc[(t, ro), col], sp) for t in tags]
            L.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
        for lab, cols, sp in RUN_STATS:
            L.append(f'| {lab} | ' + ' | '.join(' / '.join(_f(ri.loc[t, k], sp) for k in cols) for t in tags) + ' |')
        L.append('| in bands / corner / guard | ' + ' | '.join(
            f"{'yes' if ri.loc[t, 'in_bands'] else 'NO'} / {'YES' if ri.loc[t, 'corner'] else 'no'} / "
            f"{'YES' if ri.loc[t, 'guard_max_count'] else 'no'}" for t in tags) + ' |')
    return '\n'.join(L) + '\n'


def part_score1() -> dict:
    """Step 1's score (no simulation): readings1.csv, runs1.csv, verdict1.json, step1_table.md."""
    out = out_dir()
    outs, runs = collect1()
    rd = readings1(outs, runs)
    rd.to_csv(out / 'readings1.csv', index=False)
    runs.to_csv(out / 'runs1.csv', index=False)
    v = verdict1(rd, runs)
    v.update(generated=_now(), host=_host(), script=SCRIPT, commits=H.code_commits(), n_runs=int(len(runs)))
    _write_json(out / 'verdict1.json', v)
    (out / 'step1_table.md').write_text(step1_table(rd, runs, v))
    print(json.dumps({k: v[k] for k in ('class', 'no_class_runs')}, default=str), flush=True)
    print(json.dumps(v['all'], default=str)[:3000], flush=True)
    return v


# ---------------------------------------------------------------------------
# variants, commands

def write_variants() -> pd.DataFrame:
    """variants.csv (every saved run of this experiment, absolute parameters, the calyx and its wiring) and
    model_pair_distances_long.csv (raw spike counts, from pairs/) for the model-output archive and kcq/apl_effect.py;
    provenance.json."""
    from experiments import E10_joint_fit as J
    out = out_dir()
    rows = []
    for step in (1, 2):
        for r in runs_table(step).to_dict('records'):
            c = r['calyx']
            kws = J.base_kws(None, connectome=c)
            rows.append(dict(tag=r['tag'], variant=r['tag'], part=f'step{step}', base=f'bouton_{c}', k=K, connectome=c,
                             pn2kc_connections=kws.get('pn2kc_connections', 'hemibrain'), multibouton_claws=kws['multibouton_claws'],
                             seed=int(r['seed']), wiring=r.get('wiring'), arm=r.get('arm') if isinstance(r.get('arm'), str) else '',
                             kind=(f"C2's family on {SHORT[c]}, within-type claw permutation seed {int(r['seed'])}, "
                                   "fixed scales (E130's step 1)" if step == 1 else f"arm {r.get('arm')} on {SHORT[c]}"),
                             kckc_matrix=(f'{E130_EXP} matrices/kckc_{c}.parquet' if c in MALE else 'hemibrain (olfsysm)'),
                             base_params=r.get('base_params'), **{k: r.get(k) for k in PARAM_KEYS}, kc_kc_inh=r.get('wKCKC_scale'),
                             rate_on=r.get('rate_on'), rate_off=r.get('rate_off'), max_count_on=r.get('max_count_on'),
                             max_count_off=r.get('max_count_off'), realism_flag=r.get('realism_flag'), subdir=f'step{step}'))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    ps = sorted((out / 'pairs').glob('*.csv')) if (out / 'pairs').exists() else []
    keep = set(df['tag']) if len(df) else set()
    parts = [x for x in (pd.read_csv(p) for p in ps) if len(x) and x['tag'].iloc[0] in keep]
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    _write_json(out / 'provenance.json', dict(
        generated=_now(), host=_host(), script=SCRIPT, **H.code_commits(), kc_kc_tau=J.TAU, breadth_k=K,
        bases={f'bouton_{c}': J.base_kws(None, connectome=c) for c in CALYCES},
        base_params={c: e130().calyx_base(c) for c in CALYCES},
        note=("C2's family (E67) on each calyx as E130's step 1 (fixed scales, both response rates re-held), on the "
              "calyx's within-type claw permutations (seeds 1 to 7; E10_joint_fit.bouton_base's wiring "
              "'claw_permutation_within_type'); per-point provenance in provenance_<step>_<host>_<lane>.json")))
    return df


def print_commands() -> None:
    """The lanes' command lines (viki and hal 16 threads, atlas 8): viki and hal take the left's and the hemibrain's
    fourteen from either end and then the right's; atlas starts on the right's."""
    pre = 'OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    lanes = {'viki': (16, (LEFT, HEMI, RIGHT), False), 'hal': (16, (LEFT, HEMI, RIGHT), True),
             'atlas': (8, (RIGHT,), False)}
    for host, (t, cs, rev) in lanes.items():
        print(f"{host}: {pre.format(t=t)} --part step1 --calyces {' '.join(cs)}{' --reverse' if rev else ''}")


# ---------------------------------------------------------------------------
# step 2's runner (the builders are experiments/E133_builders.py)

def planned_step2(hosts: Sequence[str] = HOSTS, arms: Optional[Sequence[str]] = None, seeds: Sequence[int] = HOST_SEEDS,
                  reverse: bool = False) -> List[Point]:
    """Step 2's points: each arm on each host's four wirings (6 arms x 2 hosts x 4 wirings = 48), arm by arm."""
    arms = ('C', 'S', 'A', 'B', 'K', 'ALL') if arms is None else tuple(arms)
    pts = [Point(2, h, int(s), a) for a in arms for h in hosts for s in seeds]
    return pts[::-1] if reverse else pts


def step2_gate() -> Optional[str]:
    """Why step 2 may not run (None when it may): step 1's class must be "kept" or "mixed" (verdict1.json)."""
    cls = _read_json(out_dir() / 'verdict1.json').get('class')
    if cls not in ('kept', 'mixed'):
        return f"step 1's class is {cls!r} (step 2 runs on 'kept' or 'mixed')"
    return None


def part_step2(hosts: Sequence[str] = HOSTS, arms: Optional[Sequence[str]] = None, seeds: Sequence[int] = HOST_SEEDS,
               lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 2's points (the six arms on each host's four wirings), claimed one at a time as step 1's, on the builders
    of experiments/E133_builders.py; refused by `step2_gate`."""
    why = step2_gate()
    if why:
        raise SystemExit(f'E133 step 2 refused: {why}')
    _preload()
    from experiments import E133_builders  # noqa: F401  (loaded before the commit check)
    from experiments import E108_paa_extrapolation as E108
    cc = E108._commits_clean()
    replay_on_host(list(dict.fromkeys(list(hosts) + [DONOR[h] for h in hosts])))
    lane = lane or ('-'.join(CODE[h] for h in hosts) + ('-' + '-'.join(arms) if arms else '') + ('-rev' if reverse else ''))
    table = out_dir() / f'runs_2_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_2_{_host()}_{lane}.json'
    for pt in planned_step2(hosts, arms, seeds, reverse):
        run_point(pt, lane, table, prov_path, cc, reclaim)
    return runs_table(2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('step0', 'check', 'step1', 'smoke', 'score1', 'variants', 'step2', 'commands'))
    ap.add_argument('--calyces', nargs='+', default=list(CALYCES), choices=list(CALYCES),
                    help='step1 / check: the calyces of the lane')
    ap.add_argument('--seeds', nargs='+', type=int, default=None, help='step1: the seeds (default 1 to 7)')
    ap.add_argument('--reverse', action='store_true', help='step1 / step2: run the lane from its other end')
    ap.add_argument('--lane', default=None, help="the lane's name (default from its calyces, seeds and direction)")
    ap.add_argument('--reclaim', action='store_true', help="take over the points' claims (a process that died)")
    ap.add_argument('--hosts', nargs='+', default=list(HOSTS), choices=list(HOSTS), help='step2: the host calyces')
    ap.add_argument('--arms', nargs='+', default=None, help='step2: the arms (default all six)')
    ap.add_argument('--scratch', default=None, help='smoke: the scratch folder (default outputs/.../smoke)')
    a = ap.parse_args()
    if a.part == 'step0':
        part_step0()
    elif a.part == 'check':
        part_check(a.calyces)
    elif a.part == 'step1':
        part_step1(a.calyces, tuple(a.seeds) if a.seeds else SEEDS1, a.lane, a.reverse, a.reclaim)
    elif a.part == 'smoke':
        print(json.dumps(smoke_point(a.calyces[0], (a.seeds or [1])[0], a.scratch), default=str), flush=True)
    elif a.part == 'score1':
        part_score1()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'step2':
        part_step2(a.hosts, a.arms, tuple(a.seeds) if a.seeds else HOST_SEEDS, a.lane, a.reverse, a.reclaim)
    elif a.part == 'commands':
        print_commands()


if __name__ == '__main__':
    main()
