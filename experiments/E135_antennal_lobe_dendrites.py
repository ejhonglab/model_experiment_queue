#!/usr/bin/env python3
"""E135, the antennal lobe held to the recorded PN dendrites: the step from ORNs to PNs read glomerulus by glomerulus,
the measured step and then global parameters and terms per glomerulus fitted to the dendrites and read on held-out
flies and odors, and the recorded dendrites put in the model's PN stage.

The entry is item E135 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04, revised on moving after Tom's words of 09:38, about
09:41 and 09:55, Yang's OK of 09:43 and hal's notes on Zocchi et al. 2022). It builds on E134's shared code and does not
copy it: viki's data_prep/stage_population_data.py and its committed lock (experiments/E134_lock.json), hal's
kcq/population_fit.py (the views and one-fly units), atlas's kcq/stage_readings.py (the model read as each stage is
recorded) and experiments/E134_population_baseline.py (lanes, claims, the lock guard), viki's experiments/E134_flies.py
(the flies' recorded steps and the no-step yardstick).

The readings, fixed before any run (from the entry's Definitions):
- The shared glomeruli are the 28 that the ORN recordings, the PN-dendrite recordings and the model's input all hold
  (E90's `shared_glomeruli`, checked against the lock's populations).
- The ORN input of every arm is the model's own megamat input: the lab's recordings at the ORN axon terminals (olfdb
  dataset 1, 9 flies), converted to rates (`megamat_orn_deltas`). It already holds the lateral signal among ORNs. No arm
  uses Hallem's spiking data. On the APL16 input, paa, euc, -aPine and B-myr come from the validation2 terminal
  recordings.
- The stage's curve is one monotone curve from PN rate to dendrite signal, for all glomeruli and odors (stage_readings'
  `fit_stage_curve`, set on pooled values only).
- Imposed: on the shared glomeruli each PN's odor-evoked rate is set from the recorded dendrites, rank for rank on the
  model's own pooled distribution of PN rates. The arrangement over glomeruli and odors is the flies'; the distribution
  of rates and the time course are the model's. The other glomeruli keep the model's PN stage.
- The global parameters are pn_inhsc, pn_inhadd, ln_inhsc, ln_inhadd and the scale of the PN's saturating curve
  (olfsysm's pn.tanhsc).
- A PN's odor-evoked rate is its mean over the stage window (MODEL_WINDOW, closed) minus its mean over the settled
  pre-odor window (SETTLED_PRE_WINDOW, closed on the left), as stage_readings reads the model's PN dendrites.

The measured step (step 0, recordings only): for each shared glomerulus, the slope through the origin of its dendrite
response on its own ORN-terminal response across odors (fly means), as a multiple of the median glomerulus's (its
measured gain); how closely the dendrites follow the terminals (Pearson r over odors); the regression on its own
response and on the odor's summed terminal response (over the model input's glomeruli); the residual by odor; and the
other input glomeruli whose terminal responses that residual follows best (its candidate partners, ranked by the share
of the residual's variance they explain, with the chance level of the best of them under permuted odor labels).
Intervals from resampled flies (both stages' flies, independently). Each is given on all flies and on both halves of
the PN flies, with their agreement. It is read in dF/F as recorded and in the lock's trial-noise units (STEP_UNITS).
DM1, DL5 and VA2 are shown first. paa and B-myr have no dendrite recording (the dendrites hold the megamat panel only),
so of the three weak-input odors only 2-but is read.

The glomerulus terms: a gain is one factor per shared glomerulus on the drive its PN reads from its own ORN terminals,
between 0.25 and 4; a lateral term adds the terminal response (the ORN delta) of one other glomerulus to that drive,
with a weight between -1 and 1, at most two per glomerulus. olfsysm's LN layer reads only the mean ORN rate over
glomeruli and its PN layer reads its own glomerulus's ORN delta, so the terms act on the PN's drive only and the
normalization keeps the recorded ORN input (`PNStage`, a port of olfsysm's ORN, LN and PN layers checked against them).
What a term stands for: something between the recorded ORN terminals and the PN dendrites; no term is an identified
interaction.

Held to the dendrites: every parameter of this step is fitted to the recorded PN dendrites and to nothing downstream.
No parameter per odor is fitted, and no other stage gains one. The PN flies are split in two halves and the 17 odors
into folds, fixed in this experiment's lock (experiments/E135_lock.json, `--part lock`) before any fit. A fit sees one
half of the flies and the odors of its folds and is read on the other half, on the odors it did not see (each odor
predicted by the fit that left it out, assembled).

The rungs: R0, the model's stage as it is. R1, the global parameters. R2, R1 with the gains: first set from the
measured gains with no gain fitted, then fitted from there. R3, R2 with lateral terms, tried only among each
glomerulus's two best measured candidates and added one at a time for as long as the held-out error falls. A rung is
kept if its held-out error on the matrix itself is lower than the rung below by more than twice the standard error over
resamplings of the flies; where fitted gains are no better than the measured ones by that margin, the measured gains
are the rung. A lateral term must also gain more than the same term fitted to dendrites whose odor labels are permuted.

The gates (the entry's Design 1 and Status): step 0 reads recordings only and runs on E134's committed lock. No fit
(step 1) and no model run (step 2) before E134's measure gate reads "usable" and its megamat runs are scored
(`require_e134_gate`); a planner then re-reads the thresholds.

Step 2 (up to 73 runs at fixed scales, both response rates re-held, on the family's eight wirings): on the megamat
input R1, the highest kept rung (each odor predicted by the fit that left it out), that rung with its terms' glomerulus
labels permuted, the dendrites imposed, imposed with glomerulus labels permuted, imposed with odor labels permuted; on
the APL16 input R1, the highest kept rung, and the dendrites imposed for its 12 megamat odors; one run repeats an
imposed run by the second route. The rung's 24 runs are dropped if none is kept. Nothing of the family is refitted.

The two routes into the model's PN stage: route 1, olfsysm's option for the PN stage's rates given directly (viki's
Work claim, Yang's OK of 09:43; olfsysm b1e89d7's fit_mb_model `pn_rates`, `route1_pn_rates`; the odor's LN stays the
model's, and rates below a PN's floor are refused, see `imposed_rates`); route 2, the ORN
input solved glomerulus by glomerulus so that the unchanged PN stage returns the target rates to within SOLVE_TOL_HZ
(the imposed rates on the shared glomeruli and the model's own on the other input glomeruli; Tom, 09:38: an inversion
fixed by the recordings, not a search; `solve_input`, with its identity check `solve_identity`). A kept rung enters as
its table of PN rates by the same route, with no rank mapping.

The classes (on the eight wirings, against E134's own runs; thresholds reset to three times the standard error over
wirings in E134's table before any arm runs): "upstream", "upstream to the boutons only", "not upstream", "the global
parameters carry", "the glomerulus terms carry"; no class for a run outside a band on the APL16 input or with more than
1 % of bouton responses above 300 Hz (`classify`).

Concentrations: megamat throughout, at 10^-3 on both sides. The APL16 input's paa, euc and -aPine are the ORN data
within 1 log10 unit of the flies' concentrations; the imposed arm does not touch them, a kept rung's terms are applied
to them, and every APL16 number involves them.

Outputs: <shared>/results/E135_antennal_lobe_dendrites/ (compact: step0_*.csv, step0_summary.json, lock.json,
step1_*.csv, runs_2_*.csv, claims/, score_*.csv, provenance_*.json); outputs/E135_antennal_lobe_dendrites/ (local, large:
step 2's runs).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part step0
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part lock --n-folds K
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E135_antennal_lobe_dendrites.py --part step1|step2 [--lane <host>[-<n>]] [--reverse]
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part score|commands
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import itertools
import json
import os
import socket
import subprocess
import time
import warnings
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, stats

from kcq import harness as H
from kcq import population_fit as P

EXPERIMENT = 'E135_antennal_lobe_dendrites'
SCRIPT = 'experiments/E135_antennal_lobe_dendrites.py'
LOCK_COPY = H.HERE / 'experiments' / 'E135_lock.json'
E134 = 'E134_population_baseline'
# E134's verdict files, read before any fit or model run: the measure gate's verdict (written by E134's owner when the
# gate is read; {"verdict": "usable", ...}) and the score of its megamat runs (atlas's `--part score`)
E134_GATE_FILE = 'measure_gate.json'
E134_SCORED_FILE = 'score_classes.json'
E134_STAGE_TABLE = 'score_stage_table.csv'

ORN_KEY = ('orn', 'megamat', 'wild_type')
PN_KEY = ('pn_dendrite', 'megamat', 'wild_type')
MEGAMAT_ODORS: Tuple[str, ...] = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms',
                                  'B-cit', 'Lin', 'aa', 'benz', 'va')
N_SHARED = 28
NAMED_GLOMERULI: Tuple[str, ...] = ('DM1', 'DL5', 'VA2')
WEAK_INPUT_ODORS: Tuple[str, ...] = ('paa', 'B-myr', '2-but')    # only 2-but has a dendrite recording
# the units the measured step is read in: dF/F as recorded (the first look's) and the lock's trial-noise units
STEP_UNITS: Tuple[str, ...] = ('dff', 'noise')
# the units the fits and the imposed table read the dendrites in (the lock's; a planner's question, see the report)
FIT_UNITS = 'noise'

GLOBAL_PARAMS: Tuple[str, ...] = ('pn_inhsc', 'pn_inhadd', 'ln_inhsc', 'ln_inhadd', 'pn_tanhsc')
OLFSYSM_DEFAULTS: Dict[str, float] = {'pn_inhsc': 368.6631, 'pn_inhadd': 31.4088, 'ln_inhsc': 500.0, 'ln_inhadd': 200.0,
                                      'pn_tanhsc': 5.3395}
# R1's bounds, as multiples of the defaults (log-spaced search); pn_inhadd up to the ladder's top (10^6, the gain fixed)
GLOBAL_BOUNDS: Dict[str, Tuple[float, float]] = {'pn_inhsc': (0.05, 20.0), 'pn_inhadd': (0.03, 3.2e4),
                                                 'ln_inhsc': (0.05, 20.0), 'ln_inhadd': (0.05, 20.0),
                                                 'pn_tanhsc': (0.05, 20.0)}
GAIN_BOUNDS: Tuple[float, float] = (0.25, 4.0)
LATERAL_BOUNDS: Tuple[float, float] = (-1.0, 1.0)
MAX_LATERAL = 2
N_CANDIDATES = 2
KEEP_SE_MULT = 2.0
RUNGS: Tuple[str, ...] = ('R0', 'R1', 'R2', 'R3')
SOLVE_TOL_HZ = 0.5
N_BOOT = 1000
N_BOOT_KEEP = 200
N_PERMUTATIONS = 20
N_PERM_PARTNERS = 1000
SEED = 0
CI = 0.9

# the stage windows (s), as kcq/stage_readings.py reads the model's PN dendrites
MODEL_WINDOW: Tuple[float, float] = (0.0, 0.5)
SETTLED_PRE_WINDOW: Tuple[float, float] = (-0.25, 0.0)

# step 2: the arms (input panel, arm); 8 wirings each, plus one repeat of an imposed run by the second route
Arm = namedtuple('Arm', 'panel arm')
MEGAMAT_ARMS: Tuple[str, ...] = ('R1', 'rung', 'rung_glom_perm', 'imposed', 'imposed_glom_perm', 'imposed_odor_perm')
APL16_ARMS: Tuple[str, ...] = ('R1', 'rung', 'imposed')
RUNG_ARMS: Tuple[str, ...] = ('rung', 'rung_glom_perm')
REPEAT_ARM = Arm('megamat', 'imposed_route2_repeat')
WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E134's S then T

# the classes (the entry's; the thresholds are reset from E134's table before any arm runs)
MIN_WIRINGS = 6
SOMA_DROP_UPSTREAM = 0.15
SOMA_WITHIN_BOUTONS_ONLY = 0.05
BOUTON_STRENGTH_HALVED = 0.5
CARRY_SHARE = 0.5
THRESHOLD_SE_MULT = 3.0
BOUTON_RATE_CAP_HZ = 300.0
BOUTON_RATE_CAP_SHARE = 0.01
LANES: Dict[str, Tuple[int, Tuple[Tuple[str, bool], ...]]] = {
    'viki': (16, (('viki', False),)),
    'hal': (16, (('hal-1', False), ('hal-2', True))),
    'atlas': (8, (('atlas', True),)),
}


# ---------------------------------------------------------------------------
# paths, small helpers and guards

def out_dir() -> Path:
    """<shared>/results/E135_antennal_lobe_dendrites/ (created)."""
    return H.results_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _jsonable(x):
    if hasattr(x, 'item'):
        return x.item()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_jsonable) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text()) if Path(path).exists() else {}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def _nanmean(A: np.ndarray, axis: int = 0) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(A, axis=axis)


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git with no uncommitted change."""
    try:
        r = subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def require_e134_lock() -> dict:
    """E134's committed lock (atlas's `E134_population_baseline.require_lock`, which runs viki's check); every part reads
    the recordings through it."""
    from experiments import E134_population_baseline as E134R
    return E134R.require_lock()


def e134_gate_status(results: Optional[Path] = None) -> dict:
    """What E134 has recorded of its measure gate: the verdict file's contents (E134_GATE_FILE) and whether its megamat
    runs are scored (E134_SCORED_FILE), in <shared>/results/E134_population_baseline/ (or `results`)."""
    d = Path(results) if results is not None else H.shared_dir() / 'results' / E134
    gate = _read_json(d / E134_GATE_FILE)
    return dict(dir=str(d), gate_file=E134_GATE_FILE, gate=gate, verdict=gate.get('verdict'),
                scored=(d / E134_SCORED_FILE).exists(), scored_file=E134_SCORED_FILE)


def require_e134_gate(results: Optional[Path] = None) -> dict:
    """E134's measure gate, read before any fit or model run: refuses (PermissionError, naming what it waits for) until
    E134 has recorded that its measure is "usable" and its megamat runs are scored."""
    st = e134_gate_status(results)
    if st['verdict'] != 'usable' or not st['scored']:
        raise PermissionError(
            "E135's fits and model runs wait on E134's measure gate (the entry's Design 1 and Status): they need "
            f"{st['dir']}/{E134_GATE_FILE} with \"verdict\": \"usable\" (now: {st['verdict']!r}), written by E134's "
            f"owner once the gate is read, and {st['dir']}/{E134_SCORED_FILE} (E134's megamat runs scored, its runner's "
            f"--part score; now: {'present' if st['scored'] else 'absent'}). A planner then re-reads E135's thresholds. "
            'Step 0 (recordings only) does not wait.')
    return st


# ---------------------------------------------------------------------------
# step 0: the recorded step (recordings only)

def shared_glomeruli(check: bool = True) -> List[str]:
    """The 28 glomeruli the ORN recordings, the PN dendrites and the model's input share (E90's `shared_glomeruli`),
    checked against the lock's populations (a glomerulus any fly of each stage has) when `check`."""
    from experiments import E90_model_stages as E90
    gl = list(E90.shared_glomeruli())
    if check:
        from data_prep import stage_population_data as D
        require_e134_lock()
        have = []
        for key in (ORN_KEY, PN_KEY):
            pop = D.load_population(*key)
            have.append(set().union(*[set(map(str, t.index)) for t in pop.tables.values()]))
        model = set(H.apl16_orn_deltas().index)
        alt = sorted(have[0] & have[1] & model)
        if alt != sorted(gl):
            raise RuntimeError(f"E90's shared glomeruli {sorted(gl)} differ from the lock's populations' {alt}")
    if len(gl) != N_SHARED:
        raise RuntimeError(f'{len(gl)} shared glomeruli, not {N_SHARED}')
    return sorted(gl)


def recorded_stage(key: Tuple[str, str, str], units: str = 'dff', tol: float = 1e-9) -> Dict[str, pd.DataFrame]:
    """Per fly, units x the 17 megamat odors of trial means, in dF/F (`units='dff'`) or in the lock's trial-noise
    units ('noise'), read through the data module against E134's lock.

    The population file (sha256 checked against the lock) holds the trial-noise values only; dF/F is the same rows
    times their noise SD, from the data module's own reader of the population's source (`orn_long` / `pn_dendrite_long`),
    whose trial means in noise units must equal the population's to `tol`."""
    from data_prep import stage_population_data as D
    if units not in STEP_UNITS:
        raise ValueError(f'{units=}')
    require_e134_lock()
    pop = D.load_population(*key)
    if key == ORN_KEY:
        long = D.orn_long(D.OLFDB_DATASETS['orn_megamat'], 'megamat')
    elif key == PN_KEY:
        long = D.pn_dendrite_long()
    else:
        raise KeyError(key)
    long = long.assign(dff=long.value * long.noise_sd)
    out = {}
    for fly in pop.flies:
        sub = long[long.fly == fly]
        t = pop.tables[fly]
        g = sub.groupby(['unit', 'odor'], sort=False)
        noise = g.value.mean().unstack('odor').reindex(index=t.index, columns=list(pop.odors))
        a, b = noise.to_numpy(dtype=float), t.to_numpy(dtype=float)
        if not (np.array_equal(np.isnan(a), np.isnan(b)) and np.nanmax(np.abs(a - b)) <= tol):
            raise RuntimeError(f'{key} {fly}: the reader differs from the locked population')
        x = noise if units == 'noise' else g.dff.mean().unstack('odor').reindex(index=t.index, columns=list(pop.odors))
        x.index = x.index.astype(str)
        out[str(fly)] = x.astype(float)
    return out


def fly_stack(tables: Mapping[str, pd.DataFrame], glomeruli: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """flies x glomeruli x odors (NaN where a fly lacks a glomerulus or an odor)."""
    return np.stack([t.reindex(index=list(glomeruli), columns=list(odors)).to_numpy(dtype=float)
                     for t in tables.values()])


def _slopes(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per row: the slope through the origin of y on x, Pearson r and the number of jointly finite entries."""
    ok = np.isfinite(x) & np.isfinite(y)
    n = ok.sum(axis=1)
    x0, y0 = np.where(ok, x, 0.0), np.where(ok, y, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        slope = (x0 * y0).sum(axis=1) / (x0 * x0).sum(axis=1)
        mx, my = x0.sum(axis=1) / n, y0.sum(axis=1) / n
        dx, dy = np.where(ok, x - mx[:, None], 0.0), np.where(ok, y - my[:, None], 0.0)
        r = (dx * dy).sum(axis=1) / np.sqrt((dx * dx).sum(axis=1) * (dy * dy).sum(axis=1))
    slope = np.where(n >= 3, slope, np.nan)
    r = np.where(n >= 3, r, np.nan)
    return slope, r, n


def _relative(s: np.ndarray) -> np.ndarray:
    m = np.nanmedian(s[np.isfinite(s) & (s > 0)]) if np.any(np.isfinite(s) & (s > 0)) else np.nan
    return s / m


def measured_step(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                  total: Optional[np.ndarray] = None) -> pd.DataFrame:
    """The measured step on fly means.

    Args:
        orn, pn: flies x glomeruli x odors of the ORN terminals and the PN dendrites (the same glomeruli and odors).
        glomeruli, odors: their labels.
        total: the odor's summed terminal response (odors,), for the regression on own and summed response.

    Returns:
        Per glomerulus: slope (through the origin, PN on own ORN), gain (the slope over the median glomerulus's), r
        (Pearson over odors), n_odors; with `total` also own_coef, total_coef, intercept (least squares of PN on own
        ORN, the summed ORN and 1), r2_own (r squared) and r2_own_total.
    """
    x, y = _nanmean(orn), _nanmean(pn)
    slope, r, n = _slopes(x, y)
    df = pd.DataFrame(dict(slope=slope, gain=_relative(slope), r=r, n_odors=n), index=pd.Index(list(glomeruli),
                      name='glomerulus'))
    if total is not None:
        rows = []
        for i in range(len(glomeruli)):
            ok = np.isfinite(x[i]) & np.isfinite(y[i]) & np.isfinite(total)
            if ok.sum() < 4:
                rows.append((np.nan,) * 4)
                continue
            A = np.column_stack([x[i, ok], total[ok], np.ones(ok.sum())])
            coef, *_ = np.linalg.lstsq(A, y[i, ok], rcond=None)
            res = y[i, ok] - A @ coef
            ss = ((y[i, ok] - y[i, ok].mean()) ** 2).sum()
            rows.append((coef[0], coef[1], coef[2], 1.0 - (res ** 2).sum() / ss if ss > 0 else np.nan))
        R = np.array(rows)
        df['own_coef'], df['total_coef'], df['intercept'], df['r2_own_total'] = R.T
        df['r2_own'] = df.r ** 2
    return df


def residuals(orn: np.ndarray, pn: np.ndarray, step: pd.DataFrame, glomeruli: Sequence[str],
              odors: Sequence[str]) -> pd.DataFrame:
    """The residual by odor of each glomerulus's dendrites (fly mean) about its own slope (glomeruli x odors)."""
    x, y = _nanmean(orn), _nanmean(pn)
    s = step.loc[list(glomeruli), 'slope'].to_numpy(dtype=float)
    return pd.DataFrame(y - s[:, None] * x, index=pd.Index(list(glomeruli), name='glomerulus'), columns=list(odors))


def _corr_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Pearson r between every row of A and every row of B (no NaN), A.shape[0] x B.shape[0]."""
    a = A - A.mean(axis=1, keepdims=True)
    b = B - B.mean(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (a @ b.T) / np.sqrt((a * a).sum(axis=1)[:, None] * (b * b).sum(axis=1)[None, :])


def candidate_partners(resid: pd.DataFrame, orn_mean: pd.DataFrame, own: Optional[pd.DataFrame] = None,
                       dendrites: Optional[pd.DataFrame] = None, n: int = N_CANDIDATES,
                       n_perm: int = N_PERM_PARTNERS, seed: int = SEED) -> pd.DataFrame:
    """Per glomerulus, the other glomeruli whose terminal responses its residual follows best.

    Args:
        resid: glomeruli x odors, the residuals (`residuals`).
        orn_mean: partners x odors, the terminal responses (fly means) of the candidate partners (the model input's
            glomeruli the ORN recordings have); a glomerulus is never its own partner.
        own, dendrites: glomeruli x odors of the own terminal responses and the dendrites (fly means), for each pair's
            joint least-squares weights (dendrites on own and partner, through the origin).
        n: partners kept per glomerulus.
        n_perm: permutations of the residual's odor labels for the chance level of the best partner's r squared.
        seed: of the permutations.

    Returns:
        Rows glomerulus, rank (1 = best), partner, r, r2, chance_r2_95 (the 95th percentile of the best partner's r
        squared under permuted odor labels), and with `own` and `dendrites` own_weight, partner_weight and
        weight_rel (partner over own weight).
    """
    odors = list(resid.columns)
    X = orn_mean.reindex(columns=odors).to_numpy(dtype=float)
    partners = list(orn_mean.index.astype(str))
    rng = np.random.default_rng(seed)
    rows = []
    for g in resid.index:
        e = resid.loc[g].to_numpy(dtype=float)
        ok = np.isfinite(e) & np.isfinite(X).all(axis=0)
        cand = [i for i, p in enumerate(partners) if p != str(g)]
        if ok.sum() < 4 or not cand:
            continue
        Xc = X[np.ix_(cand, np.flatnonzero(ok))]
        r = _corr_rows(e[ok][None, :], Xc)[0]
        if n_perm > 0:
            perm = np.array([rng.permutation(e[ok]) for _ in range(n_perm)])
            chance = float(np.nanpercentile(np.nanmax(_corr_rows(perm, Xc) ** 2, axis=1), 95))
        else:
            chance = np.nan
        order = np.argsort(-(r ** 2))
        for k, j in enumerate(order[:n]):
            row = dict(glomerulus=str(g), rank=k + 1, partner=partners[cand[j]], r=float(r[j]), r2=float(r[j] ** 2),
                       chance_r2_95=chance)
            if own is not None and dendrites is not None:
                xo = own.loc[g, odors].to_numpy(dtype=float)
                y = dendrites.loc[g, odors].to_numpy(dtype=float)
                ok2 = ok & np.isfinite(xo) & np.isfinite(y)
                A = np.column_stack([xo[ok2], X[cand[j], ok2]])
                w, *_ = np.linalg.lstsq(A, y[ok2], rcond=None)
                row.update(own_weight=float(w[0]), partner_weight=float(w[1]),
                           weight_rel=float(w[1] / w[0]) if w[0] != 0 else np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def _total(orn_all: np.ndarray) -> np.ndarray:
    """The odor's summed terminal response over glomeruli (fly means; glomeruli no fly has are left out)."""
    return np.nansum(_nanmean(orn_all), axis=0)


def step_intervals(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                   n_boot: int = N_BOOT, seed: int = SEED, ci: float = CI) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Intervals from resampling the flies of both stages (each with replacement, independently), centered on the
    point values as the analysis of 2026-11-17 centers them (`centered`: the point minus the draws' mean less their
    lower quantile, the point plus their upper quantile less their mean): a resample with repeated flies has noisier fly
    means, so its r sits below the measured one and raw percentiles would not hold the point.

    Returns:
        (per glomerulus gain_lo, gain_hi, r_lo, r_hi, slope_lo, slope_hi, gain_boot_sd, log2_gain_boot_sd;
        per glomerulus x odor resid_lo, resid_hi), at coverage `ci`.
    """
    rng = np.random.default_rng(seed)
    G, O = len(glomeruli), len(odors)
    gains, rs, sl = np.full((n_boot, G), np.nan), np.full((n_boot, G), np.nan), np.full((n_boot, G), np.nan)
    res = np.full((n_boot, G, O), np.nan)
    for b in range(n_boot):
        io = rng.integers(0, orn.shape[0], orn.shape[0])
        ip = rng.integers(0, pn.shape[0], pn.shape[0])
        x, y = _nanmean(orn[io]), _nanmean(pn[ip])
        s, r, _ = _slopes(x, y)
        gains[b], rs[b], sl[b] = _relative(s), r, s
        res[b] = y - s[:, None] * x
    q = ((1 - ci) / 2, 1 - (1 - ci) / 2)
    x0, y0 = _nanmean(orn), _nanmean(pn)
    s0, r0, _ = _slopes(x0, y0)
    point = dict(gain=_relative(s0), r=r0, slope=s0)
    res0 = y0 - s0[:, None] * x0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        cols = {}
        for k, B in (('gain', gains), ('r', rs), ('slope', sl)):
            cols[f'{k}_lo'], cols[f'{k}_hi'] = centered(point[k], B, q)
        cols['r_lo'], cols['r_hi'] = np.clip(cols['r_lo'], -1.0, 1.0), np.clip(cols['r_hi'], -1.0, 1.0)
        g = pd.DataFrame(dict(cols, gain_boot_sd=np.nanstd(gains, axis=0),
                              log2_gain_boot_sd=np.nanstd(np.log2(np.where(gains > 0, gains, np.nan)), axis=0)),
                         index=pd.Index(list(glomeruli), name='glomerulus'))
        lo, hi = centered(res0, res, q)
    e = pd.DataFrame(dict(glomerulus=np.repeat(list(glomeruli), O), odor=np.tile(list(odors), G),
                          resid_lo=lo.ravel(), resid_hi=hi.ravel()))
    return g, e


def centered(point: np.ndarray, draws: np.ndarray, q: Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]:
    """The interval of the draws (axis 0) centered on the point: point - (mean - lower quantile), point + (upper
    quantile - mean)."""
    m = np.nanmean(draws, axis=0)
    return point - (m - np.nanquantile(draws, q[0], axis=0)), point + (np.nanquantile(draws, q[1], axis=0) - m)


def halves(pn_flies: Sequence[str], seed: int = SEED) -> Tuple[List[str], List[str]]:
    """The two halves of the PN flies (sizes floor(n/2) and the rest; hal's `split_halves` with one split), each in the
    flies' order."""
    a, b = P.split_halves(len(pn_flies), n_splits=1, seed=seed)[0]
    return [pn_flies[i] for i in a], [pn_flies[i] for i in b]


def halves_agreement(sa: pd.DataFrame, sb: pd.DataFrame) -> Dict[str, float]:
    """The agreement of two halves' measured steps: Pearson r and Spearman ρ of log2 gains over glomeruli, the RMS
    difference of log2 gains, and Pearson r of the halves' r over glomeruli."""
    ga, gb = np.log2(sa.gain.where(sa.gain > 0)), np.log2(sb.gain.where(sb.gain > 0))
    ok = ga.notna() & gb.notna()
    ra, rb = sa.r, sb.r
    okr = ra.notna() & rb.notna()
    return dict(log2_gain_pearson=float(np.corrcoef(ga[ok], gb[ok])[0, 1]) if ok.sum() > 2 else np.nan,
                log2_gain_spearman=float(stats.spearmanr(ga[ok], gb[ok]).correlation) if ok.sum() > 2 else np.nan,
                log2_gain_rms_diff=float(np.sqrt(np.mean((ga[ok] - gb[ok]) ** 2))) if ok.any() else np.nan,
                r_pearson=float(np.corrcoef(ra[okr], rb[okr])[0, 1]) if okr.sum() > 2 else np.nan,
                n_glomeruli=int(ok.sum()))


def all_splits_agreement(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                         total: Optional[np.ndarray] = None) -> pd.DataFrame:
    """`halves_agreement` over every split of the PN flies into floor(n/2) and the rest (35 for 7 flies)."""
    n = pn.shape[0]
    rows = []
    for a in itertools.combinations(range(n), n // 2):
        b = [i for i in range(n) if i not in a]
        sa = measured_step(orn, pn[list(a)], glomeruli, odors)
        sb = measured_step(orn, pn[b], glomeruli, odors)
        rows.append(dict(half_a=list(a), **halves_agreement(sa, sb)))
    return pd.DataFrame(rows)


def follows_rule(step: pd.DataFrame) -> pd.Series:
    """PROVISIONAL (a planner's question): a glomerulus's dendrites follow its terminals where the 90 % interval of r
    over resampled flies lies above 0. Where they do not, its measured gain is left at 1."""
    return step['r_lo'] > 0


def no_step_rows() -> pd.DataFrame:
    """E134's no-step yardstick and recorded-step summary of the ORN -> PN dendrite step (viki's step 1, the flies
    alone): flies_no_step.csv and flies_steps_summary.csv rows of that step."""
    d = H.shared_dir() / 'results' / E134
    step = f'{"__".join(ORN_KEY)} -> {"__".join(PN_KEY)}'
    out = []
    for name in ('flies_no_step.csv', 'flies_steps_summary.csv'):
        p = d / name
        if not p.exists():
            continue
        t = pd.read_csv(p)
        t = t[t.step == step].copy()
        t.insert(0, 'source', name)
        out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def read_units(units: str, glomeruli: Sequence[str], odors: Sequence[str] = MEGAMAT_ODORS) -> dict:
    """The recordings of both stages in `units`: per-fly tables, the stacks on the shared glomeruli, the ORN stack on
    the model input's glomeruli the ORN recordings have (the partners and the summed response), and the flies."""
    orn_t = recorded_stage(ORN_KEY, units)
    pn_t = recorded_stage(PN_KEY, units)
    model_gl = list(H.apl16_orn_deltas().index)
    orn_have = sorted(set().union(*[set(t.index) for t in orn_t.values()]))
    partners = [g for g in model_gl if g in orn_have]
    return dict(orn_tables=orn_t, pn_tables=pn_t, orn=fly_stack(orn_t, glomeruli, odors),
                pn=fly_stack(pn_t, glomeruli, odors), orn_all=fly_stack(orn_t, partners, odors), partners=partners,
                orn_flies=list(orn_t), pn_flies=list(pn_t))


def step0_units(units: str, glomeruli: Sequence[str], odors: Sequence[str] = MEGAMAT_ODORS, n_boot: int = N_BOOT,
                seed: int = SEED, data: Optional[dict] = None) -> dict:
    """Step 0 in one set of units (`part_step0`'s tables for it)."""
    R = read_units(units, glomeruli, odors) if data is None else data
    orn, pn, orn_all = R['orn'], R['pn'], R['orn_all']
    total = _total(orn_all)
    A, B = halves(R['pn_flies'], seed)
    ia = [R['pn_flies'].index(f) for f in A]
    ib = [R['pn_flies'].index(f) for f in B]
    steps, partners_rows = [], []
    allstep = measured_step(orn, pn, glomeruli, odors, total)
    iv, ev = step_intervals(orn, pn, glomeruli, odors, n_boot, seed)
    allstep = allstep.join(iv)
    allstep['follows_provisional'] = follows_rule(allstep)
    allstep['gain_for_r2_provisional'] = allstep.gain.where(allstep.follows_provisional, 1.0)
    allstep['gain_interval_excludes_1'] = (allstep.gain_lo > 1) | (allstep.gain_hi < 1)
    allstep['rank_gain'] = allstep.gain.rank(ascending=False, method='min')
    allstep['rank_r'] = allstep.r.rank(ascending=False, method='min')
    sets = {'all': (orn, pn, allstep), 'half_a': (orn, pn[ia], None), 'half_b': (orn, pn[ib], None)}
    orn_mean = pd.DataFrame(_nanmean(orn_all), index=R['partners'], columns=list(odors))
    for name, (o, p, st) in sets.items():
        if st is None:
            st = measured_step(o, p, glomeruli, odors, total)
        res = residuals(o, p, st, glomeruli, odors)
        own = pd.DataFrame(_nanmean(o), index=list(glomeruli), columns=list(odors))
        den = pd.DataFrame(_nanmean(p), index=list(glomeruli), columns=list(odors))
        cp = candidate_partners(res, orn_mean, own, den, N_CANDIDATES, N_PERM_PARTNERS, seed)
        cp.insert(0, 'flies', name)
        partners_rows.append(cp)
        st = st.reset_index()
        st.insert(0, 'flies', name)
        steps.append(st)
        if name == 'all':
            resid_all, own_all, den_all = res, own, den
    step_df = pd.concat(steps, ignore_index=True)
    step_df.insert(0, 'units', units)
    part_df = pd.concat(partners_rows, ignore_index=True)
    part_df.insert(0, 'units', units)
    rl = resid_all.stack().rename('resid').reset_index()
    rl.columns = ['glomerulus', 'odor', 'resid']
    rl['orn'] = own_all.stack().to_numpy()
    rl['pn'] = den_all.stack().to_numpy()
    rl['fit'] = rl.pn - rl.resid
    rl = rl.merge(ev, on=['glomerulus', 'odor'], how='left')
    rl['departs'] = (rl.resid_lo > 0) | (rl.resid_hi < 0)
    rl.insert(0, 'units', units)
    sa = step_df[step_df.flies == 'half_a'].set_index('glomerulus')
    sb = step_df[step_df.flies == 'half_b'].set_index('glomerulus')
    agree = halves_agreement(sa, sb)
    spl = all_splits_agreement(orn, pn, glomeruli, odors)
    agree_rows = []
    for k in ('log2_gain_pearson', 'log2_gain_spearman', 'log2_gain_rms_diff', 'r_pearson'):
        agree_rows.append(dict(units=units, measure=k, locked_split=agree[k], mean_all_splits=float(spl[k].mean()),
                               min_all_splits=float(spl[k].min()), max_all_splits=float(spl[k].max()),
                               n_splits=int(len(spl)), half_a=';'.join(A), half_b=';'.join(B)))
    return dict(step=step_df, partners=part_df, residuals=rl, agreement=pd.DataFrame(agree_rows),
                halves=dict(half_a=A, half_b=B), n_orn_flies=int(orn.shape[0]), n_pn_flies=int(pn.shape[0]),
                partners_pool=R['partners'])


def named_rows(step: pd.DataFrame, resid: pd.DataFrame, partners: pd.DataFrame, units: str,
               names: Sequence[str] = NAMED_GLOMERULI, top: int = 4) -> List[dict]:
    """The named glomeruli's readings (all flies): where each stands among the 28, its largest residuals with
    intervals, and its candidate partners."""
    s = step[(step.units == units) & (step.flies == 'all')].set_index('glomerulus')
    out = []
    for g in names:
        if g not in s.index:
            out.append(dict(glomerulus=g, present=False))
            continue
        r = s.loc[g]
        rr = resid[(resid.units == units) & (resid.glomerulus == g)].copy()
        rr['abs'] = rr.resid.abs()
        rr = rr.sort_values('abs', ascending=False).head(top)
        cp = partners[(partners.units == units) & (partners.flies == 'all') & (partners.glomerulus == g)]
        out.append(dict(glomerulus=g, present=True, units=units, gain=float(r.gain), gain_lo=float(r.gain_lo),
                        gain_hi=float(r.gain_hi), r=float(r.r), r_lo=float(r.r_lo), r_hi=float(r.r_hi),
                        rank_gain=int(r.rank_gain), rank_r=int(r.rank_r), n_glomeruli=int(s.gain.notna().sum()),
                        total_coef=float(r.total_coef), r2_own=float(r.r2_own), r2_own_total=float(r.r2_own_total),
                        follows_provisional=bool(r.follows_provisional),
                        largest_residuals=[dict(odor=o, resid=float(v), lo=float(lo), hi=float(hi), departs=bool(d))
                                           for o, v, lo, hi, d in zip(rr.odor, rr.resid, rr.resid_lo, rr.resid_hi,
                                                                      rr.departs)],
                        partners=[dict(partner=p, r=float(x), r2=float(x2), chance_r2_95=float(c))
                                  for p, x, x2, c in zip(cp.partner, cp.r, cp.r2, cp.chance_r2_95)]))
    return out


def part_step0(units: Sequence[str] = STEP_UNITS, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Step 0: the measured step on all flies and on both halves with their agreement, the named glomeruli first, the
    glomeruli that do not follow their terminals, the entries that depart beyond the flies' spread, E134's no-step
    yardstick, and what cannot be read. Writes step0_measured_step.csv, step0_partners.csv, step0_residuals.csv,
    step0_halves.csv, step0_no_step.csv, step0_summary.json and provenance_step0.json.

    The saved PN stages of the earlier experiments against the dendrites (the entry's step 0, item 2) are
    `saved_stage_rows`, a separate part (`--part step0_saved`), since it reads saved model outputs."""
    t0 = time.time()
    lock = require_e134_lock()
    glom = shared_glomeruli()
    res = {u: step0_units(u, glom, MEGAMAT_ODORS, n_boot, seed) for u in units}
    out = out_dir()
    step = pd.concat([r['step'] for r in res.values()], ignore_index=True)
    partners = pd.concat([r['partners'] for r in res.values()], ignore_index=True)
    resid = pd.concat([r['residuals'] for r in res.values()], ignore_index=True)
    agree = pd.concat([r['agreement'] for r in res.values()], ignore_index=True)
    nos = no_step_rows()
    step.to_csv(out / 'step0_measured_step.csv', index=False)
    partners.to_csv(out / 'step0_partners.csv', index=False)
    resid.to_csv(out / 'step0_residuals.csv', index=False)
    agree.to_csv(out / 'step0_halves.csv', index=False)
    nos.to_csv(out / 'step0_no_step.csv', index=False)
    summ = dict(experiment=EXPERIMENT, part='step0', host=_host(), time=_now(), glomeruli=glom,
                odors=list(MEGAMAT_ODORS), units=list(units), n_boot=n_boot, seed=seed,
                halves=res[units[0]]['halves'], n_orn_flies=res[units[0]]['n_orn_flies'],
                n_pn_flies=res[units[0]]['n_pn_flies'], partners_pool=res[units[0]]['partners_pool'],
                not_readable=dict(odors=[o for o in WEAK_INPUT_ODORS if o not in MEGAMAT_ODORS],
                                  read=[o for o in WEAK_INPUT_ODORS if o in MEGAMAT_ODORS],
                                  why='the dendrites hold the megamat panel only'),
                provisional=dict(follows='the 90 % interval of r over resampled flies above 0 (`follows_rule`)',
                                 total='the summed terminal response over the model input glomeruli the ORN '
                                       'recordings have',
                                 flies_noise='fly means (the noise between flies averaged out)'))
    for u in units:
        s = step[(step.units == u) & (step.flies == 'all')].set_index('glomerulus')
        r = resid[resid.units == u]
        summ[u] = dict(
            named=named_rows(step, resid, partners, u),
            not_following=sorted(s.index[~s.follows_provisional.astype(bool)]),
            gain_departs_from_median=sorted(s.index[s.gain_interval_excludes_1.astype(bool)]),
            entries_departing=int(r.departs.sum()), n_entries=int(r.resid.notna().sum()),
            departing_by_glomerulus=r[r.departs].groupby('glomerulus').odor.apply(list).to_dict(),
            gain_range=[float(s.gain.min()), float(s.gain.max())],
            partners_above_chance=sorted(partners[(partners.units == u) & (partners.flies == 'all')
                                                  & (partners['rank'] == 1)
                                                  & (partners.r2 > partners.chance_r2_95)].glomerulus),
            median_r=float(s.r.median()),
            agreement=agree[agree.units == u].set_index('measure').locked_split.to_dict(),
            two_but=r[r.odor == '2-but'][['glomerulus', 'resid', 'resid_lo', 'resid_hi', 'departs']].pipe(
                lambda d: d[d.departs]).to_dict('records'))
    _write_json(out / 'step0_summary.json', summ)
    from data_prep import stage_population_data as D
    pops = lock['populations']
    prov = dict(script=H.script_path(__file__), part='step0', host=_host(), time=_now(), seconds=time.time() - t0,
                **H.code_commits(),
                inputs=dict(lock=dict(path=str(LOCK_E134()), sha256=_sha256(LOCK_E134())),
                            populations={D.pop_name(k): pops[D.pop_name(k)]['file'] for k in (ORN_KEY, PN_KEY)},
                            sources={str(p): _sha256(p) for p in (D.olfdb_path(D.OLFDB_DATASETS['orn_megamat']),
                                                                  D.gh146_copy('GH146'))},
                            e134_flies={n: _sha256(H.shared_dir() / 'results' / E134 / n)
                                        for n in ('flies_no_step.csv', 'flies_steps_summary.csv')
                                        if (H.shared_dir() / 'results' / E134 / n).exists()}),
                n_boot=n_boot, seed=seed)
    _write_json(out / 'provenance_step0.json', prov)
    print(f'step 0: wrote {out} ({time.time() - t0:.0f} s)', flush=True)
    return summ


def LOCK_E134() -> Path:
    return H.HERE / 'experiments' / 'E134_lock.json'


def saved_stage_rows(sources: Optional[Mapping[str, pd.DataFrame]] = None, units: str = FIT_UNITS) -> pd.DataFrame:
    """The entry's step 0, item 2: saved model PN stages (glomeruli x odors of evoked rates, on the odors both have)
    against the dendrites, per source: the model's own gain per glomerulus over the median glomerulus's (its slope on the
    model's ORN input), the measured gain over it (`data_over_model`), and the matrix view's error in one-fly units
    after the stage's curve (the lock's default pn_dendrite reading, no noise). `sources` maps a label to its table
    (default the analysis of 2026-11-17's versions.csv and E90's saved PN stage of the family's real wiring, APL on);
    E102's and E114's saved PN summaries (b_pn.csv, pn_stage.csv) hold no glomerulus table and are quoted as they are
    by `saved_stage_summaries`. Reads saved model outputs: run as its own part (`--part step0_saved`).

    Returns one row per source and glomerulus (model_gain, measured_gain, data_over_model on the source's odors) and,
    per source, a row with glomerulus '_all' holding the matrix view's error and the three views' (`read_prediction`
    against all PN flies on the source's odors and the shared glomeruli)."""
    from data_prep import stage_population_data as D
    glom = shared_glomeruli()
    if sources is None:
        sources = default_saved_sources()
    pop = pn_population(units)
    orn_t = recorded_stage(ORN_KEY, units)
    rows = []
    for label, (pn, orn_in) in sources.items():
        od = [o for o in MEGAMAT_ODORS if o in pn.columns and o in orn_in.columns]
        gl = [g for g in glom if g in pn.index and g in orn_in.index]
        x = orn_in.loc[gl, od].to_numpy(dtype=float)
        y = pn.loc[gl, od].to_numpy(dtype=float)
        sm, _, _ = _slopes(x, y)
        meas = measured_step(fly_stack(orn_t, gl, od), fly_stack(pop.tables, gl, od), gl, od)
        mrel = _relative(sm)
        for i, g in enumerate(gl):
            rows.append(dict(source=label, glomerulus=g, n_odors=len(od), model_gain=float(mrel[i]),
                             measured_gain=float(meas.loc[g, 'gain']),
                             data_over_model=float(meas.loc[g, 'gain'] / mrel[i]) if mrel[i] > 0 else np.nan))
        tgt = Target.make(pop, od, gl)
        sc = read_prediction(pn.loc[gl, od], tgt)
        rows.append(dict(source=label, glomerulus='_all', n_odors=len(od),
                         **{k: v for k, v in sc.as_series().items()}))
    return pd.DataFrame(rows)


def default_saved_sources() -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    """The saved PN stages with glomerulus tables: each version of the analysis of 2026-11-17 (results/
    yang_orn_pn_gh146/versions.csv, in sample, 12 odors; its ORN input `apl16_orn_deltas`) and E90's family on the
    real wiring, APL on (stage_readings' `load_e90_run` and `stage_changes`, the settled pre-odor level; its own input),
    as (PN table, ORN input), odors abbreviated."""
    from kcq import stage_readings as R
    orn = H.apl16_orn_deltas()
    orn.columns = [H.odor_abbrev(c) for c in orn.columns]
    out = {}
    v = pd.read_csv(H.shared_dir() / 'results' / 'yang_orn_pn_gh146' / 'versions.csv')
    v = v[(v['sample'] == 'in_sample') & ~v.version.isin(['orn_input', 'orn_dff', 'gh146'])]
    for name, sub in v.groupby('version'):
        out[f'yang_orn_pn_gh146:{name}'] = (sub.pivot(index='glomerulus', columns='odor', values='value'), orn)
    run = R.load_e90_run('real', 'B')
    od = [o for o in R.run_odors(run) if o in MEGAMAT_ODORS]
    pn = R.stage_changes(run, 'pn_dendrite', 'on', od)
    pn.index = pn.index.astype(str)
    out['E90:B:real:on'] = (pn, R.stage_changes(run, 'orn', 'on', od))
    return out


def saved_stage_summaries() -> pd.DataFrame:
    """E102's and E114's saved PN-stage summaries against the dendrites (results/E102_antennal_lobe_normalization/
    b_pn.csv, d_pn_*.csv, e_pn_*.csv; results/E114_input_compression/pn_stage.csv), quoted with their file."""
    d = H.shared_dir() / 'results'
    out = []
    for f in sorted((d / 'E102_antennal_lobe_normalization').glob('*_pn*.csv')) + [d / 'E114_input_compression' /
                                                                                  'pn_stage.csv']:
        if f.exists():
            t = pd.read_csv(f, comment='#')
            t.insert(0, 'file', str(f.relative_to(d)))
            out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def part_step0_saved(units: str = FIT_UNITS) -> pd.DataFrame:
    """The saved PN stages against the dendrites (`saved_stage_rows`, `saved_stage_summaries`):
    step0_saved_stages.csv, step0_saved_summaries.csv, provenance_step0_saved.json. Reads saved model outputs, no
    simulation."""
    require_e134_lock()
    t0 = time.time()
    rows = saved_stage_rows(units=units)
    rows.to_csv(out_dir() / 'step0_saved_stages.csv', index=False)
    saved_stage_summaries().to_csv(out_dir() / 'step0_saved_summaries.csv', index=False)
    _write_json(out_dir() / 'provenance_step0_saved.json', dict(script=H.script_path(__file__), part='step0_saved',
                host=_host(), time=_now(), seconds=time.time() - t0, units=units, **H.code_commits()))
    return rows


def pn_population(units: str = FIT_UNITS):
    """viki's PN-dendrite population (lock-checked) in `units`: 'noise' as locked; 'dff' the same rows (and single
    trials) times their noise SD (`recorded_stage`'s reader), the criterion then not in the trials' units."""
    from data_prep import stage_population_data as D
    require_e134_lock()
    pop = D.load_population(*PN_KEY, trials=True)
    if units == 'noise':
        return pop
    if units != 'dff':
        raise ValueError(f'{units=}')
    long = D.pn_dendrite_long()
    recorded_stage(PN_KEY, 'noise')                          # the reader against the locked tables
    long = long.assign(value=long.value * long.noise_sd).drop(columns='noise_sd')
    spec = dict(units='dF/F', criterion=pop.criterion, window=pop.window, notes=pop.notes,
                recording_site=pop.recording_site)
    return D.population_from_long(PN_KEY, long, trials=True, spec=spec)


# ---------------------------------------------------------------------------
# this experiment's lock: the halves and the folds, before any fit

def odor_folds(n_folds: int, odors: Sequence[str] = MEGAMAT_ODORS, seed: int = SEED) -> List[List[str]]:
    """The odors split into `n_folds` folds of nearly equal size (a seeded permutation), each odor in one fold."""
    if not 2 <= int(n_folds) <= len(odors):
        raise ValueError(f'{n_folds=} for {len(odors)} odors')
    perm = np.random.default_rng(seed).permutation(len(odors))
    return [[odors[i] for i in sorted(f)] for f in np.array_split(perm, int(n_folds))]


def lock_contents(n_folds: int, seed: int = SEED) -> dict:
    """The halves of the PN flies and the odor folds (each odor left out once), with the commits and E134's lock."""
    from data_prep import stage_population_data as D
    require_e134_lock()
    flies = [str(f) for f in D.load_population(*PN_KEY).flies]
    A, B = halves(flies, seed)
    return dict(experiment=EXPERIMENT, written=_now(), host=_host(), seed=seed, pn_flies=flies,
                halves=dict(half_a=A, half_b=B), n_folds=int(n_folds), folds=odor_folds(n_folds, MEGAMAT_ODORS, seed),
                odors=list(MEGAMAT_ODORS), glomeruli=shared_glomeruli(), fit_units=FIT_UNITS,
                e134_lock_sha256=_sha256(LOCK_E134()), commits=H.code_commits())


def write_lock(n_folds: Optional[int], seed: int = SEED) -> dict:
    """Writes results/.../lock.json and its copy experiments/E135_lock.json (then commit it). The number of odor folds
    is a planner's choice: refuses without it."""
    if n_folds is None:
        raise ValueError('the number of odor folds is not fixed by the entry ("the 17 odors into folds, so that each '
                         'odor is left out once"): give --n-folds once a planner has chosen it (17: leave one out)')
    c = lock_contents(n_folds, seed)
    _write_json(out_dir() / 'lock.json', c)
    _write_json(LOCK_COPY, c)
    print(f'wrote {LOCK_COPY} and {out_dir() / "lock.json"}; commit {LOCK_COPY.name} before step 1', flush=True)
    return c


def read_lock(path: Path = LOCK_COPY) -> dict:
    """The committed lock (PermissionError if absent or uncommitted)."""
    if not path.exists():
        raise PermissionError(f'{path.name}: no lock yet (--part lock --n-folds K, then commit it); no fit before it')
    if path == LOCK_COPY and not _committed(path):
        raise PermissionError(f'{path.name} is not committed (or has uncommitted changes)')
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# step 1: the PN stage alone

@dataclass
class StageParams:
    """The olfsysm parameters of the ORN, LN and PN layers (`from_mp`); the global parameters are R1's."""
    pn_inhsc: float = 368.6631
    pn_inhadd: float = 31.4088
    pn_tanhsc: float = 5.3395
    pn_offset: float = 2.941
    pn_taum: float = 0.01
    ln_inhsc: float = 500.0
    ln_inhadd: float = 200.0
    ln_taum: float = 0.01
    ln_tauGA: float = 0.1
    ln_tauGB: float = 0.4
    ln_thr: float = 1.0
    orn_taum: float = 0.01
    dt: float = 0.0005
    pre_start: float = -2.0
    start: float = -0.5
    end: float = 0.75
    stim_start: float = 0.0
    stim_end: float = 0.5
    n_physical_gloms: int = 51
    noise_mean: float = 0.0
    noise_sd: float = 0.0

    @classmethod
    def from_mp(cls, mp) -> 'StageParams':
        return cls(pn_inhsc=mp.pn.inhsc, pn_inhadd=mp.pn.inhadd, pn_tanhsc=mp.pn.tanhsc, pn_offset=mp.pn.offset,
                   pn_taum=mp.pn.taum, ln_inhsc=mp.ln.inhsc, ln_inhadd=mp.ln.inhadd, ln_taum=mp.ln.taum,
                   ln_tauGA=mp.ln.tauGA, ln_tauGB=mp.ln.tauGB, ln_thr=mp.ln.thr, orn_taum=mp.orn.taum, dt=mp.time_dt,
                   pre_start=mp.time_pre_start, start=mp.time_start, end=mp.time_end, stim_start=mp.time_stim_start,
                   stim_end=mp.time_stim_end, n_physical_gloms=int(mp.orn.n_physical_gloms),
                   noise_mean=mp.pn.noise.mean, noise_sd=mp.pn.noise.sd)

    def with_global(self, params: Optional[Mapping[str, float]]) -> 'StageParams':
        if not params:
            return self
        bad = [k for k in params if k not in GLOBAL_PARAMS]
        if bad:
            raise KeyError(f'not global parameters: {bad}')
        return dataclasses.replace(self, **{k: float(v) for k, v in params.items()})

    def global_values(self) -> Dict[str, float]:
        return {k: float(getattr(self, k)) for k in GLOBAL_PARAMS}


def n_steps(sp: StageParams) -> int:
    """olfsysm's steps_all: unsigned((end - pre_start) / dt)."""
    return int((sp.end - sp.pre_start) / sp.dt)


def step_times(sp: StageParams) -> np.ndarray:
    """The time of every step, pre_start + k dt, rounded to 1e-9 s (pyolfsysm's `get_time_index`)."""
    return np.round(sp.pre_start + np.arange(n_steps(sp)) * sp.dt, 9)


def stim_row(sp: StageParams) -> np.ndarray:
    T = n_steps(sp)
    r = np.zeros(T)
    r[int((sp.stim_start - sp.pre_start) / sp.dt):int((sp.stim_end - sp.pre_start) / sp.dt)] = 1.0
    return r


def windows(sp: StageParams, reading: str = 'stage') -> Tuple[np.ndarray, np.ndarray]:
    """Boolean masks over the steps of a reading's window and pre-odor window.

    Args:
        sp: the layers' parameters.
        reading: 'stage', as kcq/stage_readings.py reads the model's PN dendrites (MODEL_WINDOW closed, the settled
            pre-odor window SETTLED_PRE_WINDOW closed on the left); 'preset', olfsysm's `pn_evoked_rate`, the reading of
            fit_mb_model's `pn_rates` (the steps the odor is on, [stim start, stim end), minus the steps from time_start
            to odor onset, by olfsysm's truncating step arithmetic).
    """
    T = n_steps(sp)
    if reading == 'stage':
        t = step_times(sp)
        win = (t >= MODEL_WINDOW[0] - 1e-9) & (t <= MODEL_WINDOW[1] + 1e-9)
        pre = (t >= SETTLED_PRE_WINDOW[0] - 1e-9) & (t < SETTLED_PRE_WINDOW[1] - 1e-9)
        return win, pre
    if reading == 'preset':
        start = int((sp.start - sp.pre_start) / sp.dt)
        on = int((sp.stim_start - sp.pre_start) / sp.dt)
        off = int((sp.stim_end - sp.pre_start) / sp.dt)
        k = np.arange(T)
        return (k >= on) & (k < off), (k >= start) & (k < on)
    raise ValueError(f'{reading=}')


def orn_layer(sp: StageParams, spont: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """olfsysm's sim_ORN_layer for every odor at once: odors x glomeruli x steps (delta: glomeruli x odors)."""
    T = n_steps(sp)
    s = stim_row(sp)
    odor = spont[None, :, None] + delta.T[:, :, None] * s[None, None, :]
    wsize = 0.02 / sp.dt
    ext = 2.0 / (wsize + 1.0) if wsize > 1.0 else wsize
    for i in range(1, T):
        odor[:, :, i] = ext * odor[:, :, i] + (1 - ext) * odor[:, :, i - 1]
    orn = np.empty_like(odor)
    orn[:, :, 0] = spont[None, :]
    mul = sp.dt / sp.orn_taum
    for t in range(1, T):
        orn[:, :, t] = orn[:, :, t - 1] * (1.0 - mul) + odor[:, :, t] * mul
    return orn


def ln_layer(sp: StageParams, orn: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """olfsysm's sim_LN_layer for every odor at once: (inhA, inhB), odors x steps."""
    return ln_layer_mean(sp, orn.mean(axis=1), orn.shape[1])


def ln_layer_mean(sp: StageParams, m: np.ndarray, n_gloms: int) -> Tuple[np.ndarray, np.ndarray]:
    """olfsysm's sim_LN_layer from the mean ORN rate over the `n_gloms` glomeruli (odors x steps), which is all it
    reads: (inhA, inhB), odors x steps."""
    O, T = m.shape
    scaling = float(n_gloms) / float(sp.n_physical_gloms)
    pot = np.full(O, 300.0)
    resp = np.ones(O)
    inhA = np.full((O, T), 50.0)
    inhB = np.full((O, T), 50.0)
    inh_ln = np.zeros(O)
    for t in range(1, T):
        dA = -inhA[:, t - 1] + resp
        dB = -inhB[:, t - 1] + resp
        dLN = -pot + np.power(m[:, t - 1] * scaling, 3.0) / scaling / 2.0 * inh_ln
        inhA[:, t] = inhA[:, t - 1] + dA * sp.dt / sp.ln_tauGA
        inhB[:, t] = inhB[:, t - 1] + dB * sp.dt / sp.ln_tauGB
        inh_ln = sp.ln_inhsc / (sp.ln_inhadd + inhA[:, t])
        pot = pot + dLN * sp.dt / sp.ln_taum
        resp = (pot - sp.ln_thr) * (pot > sp.ln_thr)
    return inhA, inhB


def pn_layer(sp: StageParams, spont: np.ndarray, orn: np.ndarray, inhA: np.ndarray, inhB: np.ndarray,
             M: Optional[np.ndarray] = None, home: Optional[np.ndarray] = None, full: bool = False,
             reading: str = 'stage'):
    """olfsysm's sim_PN_layer for every odor at once, on virtual PNs: PN n reads the drive orn_delta @ M[:, n] (orn_delta
    the ORN rates minus their spontaneous rates, odors x glomeruli) in place of its own glomerulus's ORN delta, with the
    spontaneous rates of its home glomerulus `home[n]`. M the identity (default) is olfsysm's PN layer.

    Returns the evoked rates (virtual PNs x odors) in `reading` (`windows`), or with `full` the rates over time (odors x
    PNs x steps)."""
    if sp.noise_sd != 0 or sp.noise_mean != 0:
        raise ValueError('the port has no PN noise (olfsysm adds it per step); the model runs with none')
    O, G, T = orn.shape
    if M is None:
        M = np.eye(G)
        home = np.arange(G)
    home = np.asarray(home, dtype=int)
    N = M.shape[1]
    spont_h = spont[home]
    spont_pn = spont_h * sp.pn_inhsc / (spont.sum() + sp.pn_inhadd)
    inh_pn_t = sp.pn_inhsc / (sp.pn_inhadd + 0.25 * inhA + 0.75 * inhB)     # odors x steps
    win, pre = windows(sp, reading)
    pn = np.broadcast_to(spont_h[None, :], (O, N)).copy()
    acc_w = pn * win[0]
    acc_p = pn * pre[0]
    out = np.empty((O, N, T)) if full else None
    if full:
        out[:, :, 0] = pn
    inh = np.zeros((O, 1))
    identity = M is not None and N == G and np.array_equal(M, np.eye(G)) and np.array_equal(home, np.arange(G))
    for t in range(1, T):
        od = orn[:, :, t - 1] - spont[None, :]
        drive = od if identity else od @ M
        d = -pn + spont_pn[None, :] + 200.0 * np.tanh((drive + sp.pn_offset) * sp.pn_tanhsc / 200.0 * inh)
        inh = inh_pn_t[:, t:t + 1]
        pn = pn + d * sp.dt / sp.pn_taum
        pn = np.where(0.0 < pn, pn, 0.0)
        if win[t]:
            acc_w += pn
        if pre[t]:
            acc_p += pn
        if full:
            out[:, :, t] = pn
    if full:
        return out
    return (acc_w / win.sum() - acc_p / pre.sum()).T


@dataclass
class Terms:
    """The glomerulus terms of a rung: gains per shared glomerulus (1 where absent) and lateral terms
    (glomerulus, partner, weight)."""
    gains: Dict[str, float] = field(default_factory=dict)
    laterals: List[Tuple[str, str, float]] = field(default_factory=list)

    def matrix(self, glomeruli: Sequence[str]) -> np.ndarray:
        """The drive matrix M (glomeruli x glomeruli): PN g reads sum_h orn_delta_h M[h, g]."""
        idx = {g: i for i, g in enumerate(glomeruli)}
        M = np.eye(len(glomeruli))
        for g, v in self.gains.items():
            M[idx[g], idx[g]] = float(v)
        n_lat: Dict[str, int] = {}
        for g, p, w in self.laterals:
            if g == p:
                raise ValueError(f'{g}: a lateral term from itself')
            n_lat[g] = n_lat.get(g, 0) + 1
            if n_lat[g] > MAX_LATERAL:
                raise ValueError(f'{g}: more than {MAX_LATERAL} lateral terms')
            M[idx[p], idx[g]] += float(w)
        return M

    def as_dict(self) -> dict:
        return dict(gains={k: float(v) for k, v in self.gains.items()},
                    laterals=[[g, p, float(w)] for g, p, w in self.laterals])

    @classmethod
    def from_dict(cls, d: Mapping) -> 'Terms':
        return cls(gains=dict(d.get('gains', {})), laterals=[tuple(x) for x in d.get('laterals', [])])


class PNStage:
    """olfsysm's ORN, LN and PN layers ported to numpy (`orn_layer`, `ln_layer`, `pn_layer`), on one ORN input, with
    the glomerulus terms on each PN's drive: the LN layer reads the recorded input (the mean ORN rate over glomeruli),
    so the normalization keeps it. Identical to olfsysm's layers with no terms (`check_identity`).

    Args:
        delta: glomeruli x odors, the ORN input's deltas in the model's glomerulus order (all of mp.orn.data.delta's
            rows, those without ORN data at 0).
        spont: the glomeruli's spontaneous ORN rates (mp.orn.data.spont).
        sp: the layers' parameters.
    """

    def __init__(self, delta: pd.DataFrame, spont: pd.Series, sp: Optional[StageParams] = None):
        self.delta = delta.astype(float)
        self.glomeruli = [str(g) for g in delta.index]
        self.odors = [str(o) for o in delta.columns]
        self.spont = spont.reindex(delta.index).to_numpy(dtype=float)
        if np.isnan(self.spont).any():
            raise ValueError('spontaneous rates missing for some glomeruli')
        self.sp = StageParams() if sp is None else sp
        self._orn: Dict[str, np.ndarray] = {}
        self._ln: Dict[tuple, Tuple[np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_mp(cls, mp, glomeruli: Sequence[str], odors: Sequence[str]) -> 'PNStage':
        """From olfsysm's ModelParams (its delta's rows are `glomeruli`, its columns `odors` in olfsysm's order)."""
        d = np.asarray(mp.orn.data.delta, dtype=float)
        s = np.asarray(mp.orn.data.spont, dtype=float).ravel()
        return cls(pd.DataFrame(d, index=list(glomeruli), columns=list(odors)), pd.Series(s, index=list(glomeruli)),
                   StageParams.from_mp(mp))

    def _key(self, delta: np.ndarray) -> str:
        return hashlib.sha1(np.ascontiguousarray(delta).tobytes()).hexdigest()

    def orn(self, delta: Optional[np.ndarray] = None) -> np.ndarray:
        d = self.delta.to_numpy() if delta is None else np.asarray(delta, dtype=float)
        k = self._key(d)
        if k not in self._orn:
            if len(self._orn) > 4:
                self._orn.clear()
            self._orn[k] = orn_layer(self.sp, self.spont, d)
        return self._orn[k]

    def ln(self, sp: StageParams, delta: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        d = self.delta.to_numpy() if delta is None else np.asarray(delta, dtype=float)
        k = (self._key(d), sp.ln_inhsc, sp.ln_inhadd, sp.ln_taum, sp.ln_tauGA, sp.ln_tauGB, sp.ln_thr,
             sp.n_physical_gloms)
        if k not in self._ln:
            if len(self._ln) > 64:
                self._ln.clear()
            self._ln[k] = ln_layer(sp, self.orn(d))
        return self._ln[k]

    def run(self, params: Optional[Mapping[str, float]] = None, terms: Optional[Terms] = None,
            delta: Optional[pd.DataFrame] = None, reading: str = 'stage') -> pd.DataFrame:
        """The PN's odor-evoked rates in `reading` (`windows`), glomeruli x odors (Hz), at the global parameters
        `params` (default the stage's) with the glomerulus terms `terms`, on the stage's input (or `delta`, which the
        LN then reads too)."""
        sp = self.sp.with_global(params)
        d = None if delta is None else delta.reindex(index=self.glomeruli, columns=self.odors).to_numpy(dtype=float)
        orn = self.orn(d)
        inhA, inhB = self.ln(sp, d)
        M = None if terms is None else terms.matrix(self.glomeruli)
        home = None if terms is None else np.arange(len(self.glomeruli))
        x = pn_layer(sp, self.spont, orn, inhA, inhB, M, home, reading=reading)
        return pd.DataFrame(x, index=self.glomeruli, columns=self.odors)

    def run_virtual(self, M: np.ndarray, home: np.ndarray, params: Optional[Mapping[str, float]] = None,
                    reading: str = 'stage') -> np.ndarray:
        """Virtual PNs (`pn_layer`) on the stage's input: evoked rates, virtual PNs x odors."""
        sp = self.sp.with_global(params)
        inhA, inhB = self.ln(sp)
        return pn_layer(sp, self.spont, self.orn(), inhA, inhB, M, home, reading=reading)

    def rate_bounds(self, params: Optional[Mapping[str, float]] = None, reading: str = 'preset',
                    delta: Optional[pd.DataFrame] = None, shift: float = 1e12) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """The evoked rates a PN can reach at the stage's LN (lo: its glomerulus's ORN delta far down, the PN silenced
        from odor onset; hi: far up, its tanh saturated), as olfsysm's `pn_rate_bounds` reads them with the LN kept:
        (lo, hi), glomeruli x odors."""
        sp = self.sp.with_global(params)
        d = None if delta is None else delta.reindex(index=self.glomeruli, columns=self.odors).to_numpy(dtype=float)
        inhA, inhB = self.ln(sp, d)
        k = orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
        rows, cols = np.arange(len(self.glomeruli)), np.arange(len(self.odors))
        one = np.ones((len(rows), len(cols)))
        lo = _rates_at(sp, self.spont, rows, cols, -shift * one, k, inhA, inhB, reading)
        hi = _rates_at(sp, self.spont, rows, cols, shift * one, k, inhA, inhB, reading)
        mk = lambda a: pd.DataFrame(a, index=self.glomeruli, columns=self.odors)
        return mk(lo), mk(hi)

    def check_identity(self, mp, rv, tol: float = 1e-9, params: Optional[Mapping[str, float]] = None) -> dict:
        """The port against olfsysm's layers (`run_ORN_LN_sims`, `run_PN_sims` on `mp`, `rv`, whose input is this
        stage's): the ORN rates, the LN's inhA and inhB and the PN rates over time, at the stage's parameters or
        `params`. Raises RuntimeError beyond `tol`."""
        import olfsysm as osm
        sp = self.sp.with_global(params)
        keep = {k: getattr(mp.pn, k) for k in ('inhsc', 'inhadd', 'tanhsc')}
        keep_ln = {k: getattr(mp.ln, k) for k in ('inhsc', 'inhadd')}
        try:
            mp.pn.inhsc, mp.pn.inhadd, mp.pn.tanhsc = sp.pn_inhsc, sp.pn_inhadd, sp.pn_tanhsc
            mp.ln.inhsc, mp.ln.inhadd = sp.ln_inhsc, sp.ln_inhadd
            osm.run_ORN_LN_sims(mp, rv)
            osm.run_PN_sims(mp, rv)
            o_ref = np.stack([np.asarray(x, dtype=float) for x in rv.orn.sims])
            a_ref = np.stack([np.asarray(x, dtype=float).ravel() for x in rv.ln.ln_inhA.sims])
            b_ref = np.stack([np.asarray(x, dtype=float).ravel() for x in rv.ln.ln_inhB.sims])
            p_ref = np.stack([np.asarray(x, dtype=float) for x in rv.pn.pn_sims])
        finally:
            for k, v in keep.items():
                setattr(mp.pn, k, v)
            for k, v in keep_ln.items():
                setattr(mp.ln, k, v)
        orn = self.orn()
        inhA, inhB = self.ln(sp)
        pn = pn_layer(sp, self.spont, orn, inhA, inhB, full=True)
        if o_ref.shape[-1] != orn.shape[-1]:
            raise RuntimeError(f'olfsysm has {o_ref.shape[-1]} steps, the port {orn.shape[-1]} (pretime removed?)')
        diffs = dict(orn=float(np.abs(o_ref - orn).max()), inhA=float(np.abs(a_ref - inhA).max()),
                     inhB=float(np.abs(b_ref - inhB).max()), pn=float(np.abs(p_ref - pn).max()))
        if hasattr(osm, 'pn_evoked_rates'):
            mp.pn.inhsc, mp.pn.inhadd, mp.pn.tanhsc = sp.pn_inhsc, sp.pn_inhadd, sp.pn_tanhsc
            try:
                e_ref = np.asarray(osm.pn_evoked_rates(mp, rv), dtype=float)
            finally:
                for k, v in keep.items():
                    setattr(mp.pn, k, v)
            e = pn_layer(sp, self.spont, orn, inhA, inhB, reading='preset')
            diffs['evoked_preset'] = float(np.abs(e_ref - e).max())
        res = dict(max_abs_diff=diffs, tol=tol, ok=all(v <= tol for v in diffs.values()), params=sp.global_values())
        if not res['ok']:
            raise RuntimeError(f'the port differs from olfsysm: {res}')
        return res


def flies_subset(pop, flies: Sequence[str]):
    """viki's D.Population restricted to `flies` (tables and trials)."""
    flies = [f for f in pop.flies if f in set(flies)]
    return dataclasses.replace(pop, flies=tuple(flies), tables={f: pop.tables[f] for f in flies},
                               trials=None if pop.trials is None else {f: pop.trials[f] for f in flies
                                                                       if f in pop.trials})


def fit_reading():
    """The reading inside a fit: the lock's pn_dendrite reading with an affine curve (smooth in the parameters), no
    noise, on the 28 shared glomeruli. PROVISIONAL (a planner's question: the fit's loss)."""
    from kcq import stage_readings as R
    return R.Reading('fit_affine_nonoise', 'pn_dendrite', curve='affine', noise='none', glomeruli='shared28')


def read_reading():
    """The reading of a held-out prediction: the lock's default pn_dendrite reading (stage_readings'
    DEFAULT_READINGS['pn_dendrite'])."""
    from kcq import stage_readings as R
    return R.DEFAULT_READINGS['pn_dendrite']


@dataclass
class Target:
    """The dendrites a fit or a reading is held to: viki's population of the flies of one half, the odors, the shared
    glomeruli, and hal's reference of those flies (in their own one-fly units unless `units_from` is given)."""
    flies: object
    odors: List[str]
    glomeruli: List[str]
    ref: P.FlyReference
    _loss: Optional[tuple] = field(default=None, repr=False)

    def loss_terms(self) -> tuple:
        """(the flies' pooled mean and variance, the mean matrix part in glomeruli x odors order, its one-fly unit),
        computed once, for `matrix_loss`."""
        if self._loss is None:
            from kcq import stage_readings as R
            fv, fw = R.pooled_fly_values(self.flies, self.odors)
            w = fw / fw.sum()
            mf = float(np.sum(w * fv))
            vf = float(np.sum(w * (fv - mf) ** 2))
            m = self.ref.mean['matrix/values'].reindex([f'{g}{P.PAIR_SEP}{o}' for g in self.glomeruli
                                                        for o in self.odors]).to_numpy(dtype=float)
            self._loss = (mf, vf, m, float(self.ref.one_fly['matrix/values']))
        return self._loss

    @classmethod
    def make(cls, flies, odors: Sequence[str], glomeruli: Sequence[str], parts: Optional[Sequence[str]] = None,
             units_from: Optional[P.FlyReference] = None) -> 'Target':
        from data_prep import stage_population_data as D
        fl = P.fly_populations(flies)
        fl = [dataclasses.replace(f, values=f.values.reindex(index=list(glomeruli)),
                                  responds=f.responds.reindex(index=list(glomeruli))) for f in fl]
        ref = P.fly_reference(fl, list(odors), D.breadth_bins(len(odors)), parts)
        if units_from is not None:
            ref = dataclasses.replace(ref, one_fly={p: units_from.one_fly[p] for p in ref.one_fly})
        return cls(flies, list(odors), list(glomeruli), ref)


def read_prediction(pred: pd.DataFrame, target: Target, reading=None, seed_index: int = 0) -> P.Score:
    """A predicted table of PN rates (glomeruli x odors) read as the dendrites are recorded (the reading's curve set on
    the prediction's pooled values against the target flies' pooled values, then its noise and criterion) and scored
    against the target in one-fly units (the matrix view and the three views)."""
    from kcq import stage_readings as R
    reading = read_reading() if reading is None else reading
    U = pred.reindex(index=target.glomeruli, columns=target.odors).astype(float)
    curve = R.set_curve(R.pooled_values([U]), target.flies, target.odors, reading, seed_index)
    pop, _, _ = R.finish_population(U, 'pn_dendrite', 'on', reading, curve, target.flies, target.odors, seed_index,
                                    'model')
    return P.score(pop, target.ref)


def matrix_loss(pred: pd.DataFrame, target: Target) -> float:
    """The fit's loss: the matrix view's error in one-fly units after an affine curve set on the pooled values (no
    noise; `fit_reading`; stage_readings' 'affine' curve: the model's pooled mean and SD onto the flies')."""
    return matrix_loss_array(pred.reindex(index=target.glomeruli, columns=target.odors).to_numpy(dtype=float), target)


def matrix_loss_array(U: np.ndarray, target: Target) -> float:
    """`matrix_loss` of an array in the target's glomeruli x odors order."""
    mf, vf, m, unit = target.loss_terms()
    x = U[np.isfinite(U)]
    sm = float(x.std())
    scale = float(np.sqrt(vf) / sm) if sm > 0 else 0.0
    d = (scale * U + (mf - scale * float(x.mean()))).ravel() - m
    ok = np.isfinite(d)
    return float(np.sqrt(np.mean(d[ok] ** 2)) / unit)


def _log_bounds(name: str) -> Tuple[float, float]:
    lo, hi = GLOBAL_BOUNDS[name]
    d = OLFSYSM_DEFAULTS[name]
    return float(np.log(lo * d)), float(np.log(hi * d))


def fit_global(stage: PNStage, target: Target, start: Optional[Mapping[str, float]] = None,
               names: Sequence[str] = GLOBAL_PARAMS, maxiter: int = 400) -> dict:
    """R1: the global parameters fitted to the target (`matrix_loss`), Nelder-Mead on their logs inside GLOBAL_BOUNDS
    (a point outside scores as 10), from `start` (default the stage's)."""
    x0 = stage.sp.global_values() if start is None else dict(start)
    lo = np.array([_log_bounds(n)[0] for n in names])
    hi = np.array([_log_bounds(n)[1] for n in names])

    def f(z):
        if np.any(z < lo) or np.any(z > hi):
            return 10.0
        p = dict(zip(names, np.exp(z)))
        return matrix_loss(stage.run(p).reindex(index=target.glomeruli), target)
    z0 = np.log([x0[n] for n in names])
    r = optimize.minimize(f, z0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-6, maxiter=maxiter))
    p = dict(zip(names, map(float, np.exp(r.x))))
    at = [n for n, z, a, b in zip(names, r.x, lo, hi) if min(z - a, b - z) < 0.02]
    return dict(params=p, loss=float(r.fun), nfev=int(r.nfev), at_bounds=at)


def gain_grid(stage: PNStage, glomeruli: Sequence[str], params: Mapping[str, float], base: Terms,
              grid: np.ndarray) -> np.ndarray:
    """Each shared glomerulus's PN rates with its gain at every grid value and its laterals as in `base` (virtual PNs,
    one batched run): glomeruli x grid x odors. PNs are independent given the LN, so one run gives every glomerulus."""
    G = stage.glomeruli
    idx = {g: i for i, g in enumerate(G)}
    M0 = base.matrix(G)
    cols, home = [], []
    for g in glomeruli:
        for v in grid:
            c = M0[:, idx[g]].copy()
            c[idx[g]] = v
            cols.append(c)
            home.append(idx[g])
    X = stage.run_virtual(np.column_stack(cols), np.array(home), params)
    return X.reshape(len(glomeruli), len(grid), -1)


def fit_gains(stage: PNStage, target: Target, params: Mapping[str, float], start: Mapping[str, float],
              base: Optional[Terms] = None, n_grid: int = 41, maxiter: int = 200) -> dict:
    """R2's fit: the gains of the shared glomeruli (log, inside GAIN_BOUNDS) fitted to the target by L-BFGS-B, with the
    PN rates interpolated from one batched run over a log-spaced grid of gains (`gain_grid`) and re-run exactly at the
    end; from `start` (the measured gains)."""
    base = Terms() if base is None else base
    gl = target.glomeruli
    grid = np.exp(np.linspace(np.log(GAIN_BOUNDS[0]), np.log(GAIN_BOUNDS[1]), n_grid))
    oi = [stage.odors.index(o) for o in target.odors]
    Y = gain_grid(stage, gl, params, base, grid)[:, :, oi]
    lg = np.log(grid)
    step = lg[1] - lg[0]
    ar = np.arange(len(gl))

    def f(z):
        pos = np.clip((z - lg[0]) / step, 0.0, len(lg) - 1.0)
        i0 = np.minimum(np.floor(pos).astype(int), len(lg) - 2)
        w = (pos - i0)[:, None]
        return matrix_loss_array(Y[ar, i0] * (1 - w) + Y[ar, i0 + 1] * w, target)
    z0 = np.log(np.clip([start.get(g, 1.0) for g in gl], *GAIN_BOUNDS))
    r = optimize.minimize(f, z0, method='L-BFGS-B', bounds=[(lg[0], lg[-1])] * len(gl),
                          options=dict(maxiter=maxiter, eps=1e-4))
    gains = {g: float(np.exp(v)) for g, v in zip(gl, r.x)}
    terms = Terms(gains=gains, laterals=list(base.laterals))
    exact = matrix_loss(stage.run(params, terms), target)
    at = [g for g, v in zip(gl, r.x) if min(v - lg[0], lg[-1] - v) < 1e-3]
    return dict(terms=terms, loss=exact, loss_interp=float(r.fun), nit=int(r.nit), at_bounds=at)


def fit_lateral(stage: PNStage, target: Target, params: Mapping[str, float], terms: Terms, g: str, partner: str,
                n_grid: int = 41) -> dict:
    """One lateral term's weight (inside LATERAL_BOUNDS) with everything else held, by a grid over the weight and a
    bounded scalar search on the interpolated rates, re-run exactly at the end."""
    G = stage.glomeruli
    idx = {x: i for i, x in enumerate(G)}
    M0 = terms.matrix(G)
    grid = np.linspace(LATERAL_BOUNDS[0], LATERAL_BOUNDS[1], n_grid)
    cols = []
    for w in grid:
        c = M0[:, idx[g]].copy()
        c[idx[partner]] += w
        cols.append(c)
    Y = stage.run_virtual(np.column_stack(cols), np.full(n_grid, idx[g]), params)
    rest = stage.run(params, terms)

    def f(w):
        X = rest.copy()
        X.loc[g] = [np.interp(w, grid, Y[:, j]) for j in range(Y.shape[1])]
        return matrix_loss(X, target)
    r = optimize.minimize_scalar(f, bounds=LATERAL_BOUNDS, method='bounded', options=dict(xatol=1e-3))
    t = Terms(gains=dict(terms.gains), laterals=list(terms.laterals) + [(g, partner, float(r.x))])
    return dict(terms=t, weight=float(r.x), loss=matrix_loss(stage.run(params, t), target))


def model_gains(stage: PNStage, params: Mapping[str, float], glomeruli: Sequence[str], odors: Sequence[str]) -> pd.Series:
    """The model stage's own gain per glomerulus: the slope through the origin of its evoked PN rate on its ORN input
    over `odors`, over the median glomerulus's."""
    x = stage.delta.loc[list(glomeruli), list(odors)].to_numpy(dtype=float)
    y = stage.run(params).loc[list(glomeruli), list(odors)].to_numpy(dtype=float)
    s, _, _ = _slopes(x, y)
    return pd.Series(_relative(s), index=list(glomeruli))


def measured_gain_start(measured: pd.DataFrame, model_rel: pd.Series) -> Dict[str, float]:
    """R2's starting gains: the measured gain over the model stage's own gain (both over their median glomerulus),
    1 where the dendrites do not follow their terminals (`follows_rule`) or either is undefined, clipped to GAIN_BOUNDS."""
    out = {}
    for g in model_rel.index:
        ok = g in measured.index and bool(measured.loc[g, 'follows_provisional'])
        v = float(measured.loc[g, 'gain'] / model_rel[g]) if ok else 1.0
        out[g] = float(np.clip(v if np.isfinite(v) and v > 0 else 1.0, *GAIN_BOUNDS))
    return out


def keep_rung(err_rung: np.ndarray, err_below: np.ndarray, point_rung: float, point_below: float,
              mult: float = KEEP_SE_MULT) -> Tuple[bool, float]:
    """Kept if the rung's held-out error is lower than the rung below's by more than `mult` standard errors of the
    paired difference over resamplings of the flies. Returns (kept, the standard error)."""
    d = np.asarray(err_below, dtype=float) - np.asarray(err_rung, dtype=float)
    se = float(np.nanstd(d, ddof=1))
    return bool((point_below - point_rung) > mult * se), se


def permute_odor_labels(flies, rng: np.random.Generator, odors: Optional[Sequence[str]] = None):
    """viki's population with its odor labels permuted among `odors` (default all the panel's; one permutation of the
    columns, for every fly and its trials), the other columns kept."""
    odors = list(flies.odors if odors is None else odors)
    perm = [odors[i] for i in rng.permutation(len(odors))]
    ren = dict(zip(perm, odors))
    tables = {}
    for f, t in flies.tables.items():
        x = t.copy()
        x[odors] = t[perm].to_numpy()
        tables[f] = x
    trials = None
    if flies.trials:
        trials = {}
        for f, t in flies.trials.items():
            x = t.copy()
            x.columns = pd.MultiIndex.from_tuples([(ren.get(o, o), r) for o, r in t.columns], names=t.columns.names)
            trials[f] = x.sort_index(axis=1)
    return dataclasses.replace(flies, tables=tables, trials=trials)


def permute_glomerulus_labels(x: pd.DataFrame, glomeruli: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    """A glomeruli x odors table with its shared glomeruli's rows permuted among themselves."""
    g = list(glomeruli)
    perm = [g[i] for i in rng.permutation(len(g))]
    out = x.copy()
    out.loc[g] = x.loc[perm].to_numpy()
    return out


def permute_terms(terms: Terms, glomeruli: Sequence[str], rng: np.random.Generator) -> Terms:
    """A rung's terms with their glomerulus labels permuted among the shared glomeruli (gains and laterals' targets
    move together; a lateral's partner keeps its label)."""
    g = list(glomeruli)
    m = dict(zip(g, [g[i] for i in rng.permutation(len(g))]))
    lat = [(m[a], b, w) for a, b, w in terms.laterals if m[a] != b]
    return Terms(gains={m[k]: v for k, v in terms.gains.items()}, laterals=lat)


def held_out_draws(pred: pd.DataFrame, target: Target, n_boot: int = N_BOOT_KEEP, seed: int = SEED,
                   reading=None) -> np.ndarray:
    """The matrix view's error of a held-out prediction against resampled sets of the target's flies (with replacement)
    in the measured set's one-fly units, one per draw (the same draws for every prediction at the same seed)."""
    from kcq import stage_readings as R
    reading = read_reading() if reading is None else reading
    U = pred.reindex(index=target.glomeruli, columns=target.odors).astype(float)
    curve = R.set_curve(R.pooled_values([U]), target.flies, target.odors, reading, 0)
    pop, _, _ = R.finish_population(U, 'pn_dendrite', 'on', reading, curve, target.flies, target.odors, 0, 'model')
    v = P.population_views(pop, target.odors, target.ref.breadth_edges, ['matrix/values'])['matrix/values']
    A = target.ref.parts['matrix/values'].to_numpy(dtype=float)
    m = v.reindex(target.ref.parts['matrix/values'].columns).to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    out = np.empty(n_boot)
    for b in range(n_boot):
        mb = _nanmean(A[rng.integers(0, n, n)])
        d = m - mb
        out[b] = np.sqrt(np.nanmean(d ** 2)) / target.ref.one_fly['matrix/values']
    return out


def fold_map(table: pd.DataFrame, target: Target) -> pd.DataFrame:
    """A fold's table of PN rates in the fit flies' units: through the fit's own affine curve (`matrix_loss`'s), set on
    the fold's training odors only, applied to every odor. A curve absorbs the rates' scale, so the global parameters
    can change it from fold to fold; the folds' held-out columns are comparable only through these maps, as in the
    analysis of 2026-11-17 (its `loo`)."""
    mf, vf, _, _ = target.loss_terms()
    x = table.reindex(index=target.glomeruli, columns=target.odors).to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    sm = float(x.std())
    scale = float(np.sqrt(vf) / sm) if sm > 0 else 0.0
    return table * scale + (mf - scale * float(x.mean()))


def cross_fit(fit_one: Callable[[Target, List[str]], Tuple[pd.DataFrame, dict]], folds: Sequence[Sequence[str]],
              fit_targets: Mapping[str, Callable[[List[str]], Target]], read_targets: Mapping[str, Target],
              directions: Sequence[Tuple[str, str]] = (('half_a', 'half_b'), ('half_b', 'half_a')),
              n_boot: int = N_BOOT_KEEP, seed: int = SEED, map_folds: bool = True) -> dict:
    """A rung fitted fold by fold on one half and read on the other: per direction (fit half, read half), for each
    fold the fit on the fit half's flies and the other folds' odors (`fit_one(target, train_odors)` -> (the table of PN
    rates on every odor, the fit's record)), the held-out columns assembled into one table (each odor by the fit that
    left it out; with `map_folds` each fold's table first in the fit flies' units, `fold_map`), read on the read half
    (`read_prediction`) with bootstrap draws of its matrix error (`held_out_draws`). Returns per direction the table
    (and 'table_hz', the folds' own rates assembled unmapped), the score, the draws and the folds' records, and the mean
    over directions ('point', 'draws')."""
    out = {}
    for fit_half, read_half in directions:
        cols, cols_hz, recs = {}, {}, []
        for k, fold in enumerate(folds):
            train = [o for f in folds for o in f if o not in set(fold)]
            t = fit_targets[fit_half](train)
            tab, rec = fit_one(t, train)
            mapped = fold_map(tab, t) if map_folds else tab
            for o in fold:
                cols[o] = mapped[o]
                cols_hz[o] = tab[o]
            recs.append(dict(fold=k, held_out=list(fold), **rec))
        tgt = read_targets[read_half]
        table = pd.DataFrame(cols)[tgt.odors]
        sc = read_prediction(table, tgt)
        draws = held_out_draws(table, tgt, n_boot, seed)
        out[f'{fit_half}->{read_half}'] = dict(table=table, table_hz=pd.DataFrame(cols_hz)[tgt.odors], score=sc,
                                                draws=draws, folds=recs, point=float(sc.parts['matrix/values']))
    keys = list(out)
    out['point'] = float(np.mean([out[k]['point'] for k in keys]))
    out['draws'] = np.mean([out[k]['draws'] for k in keys], axis=0)
    return out


N_PERM_CONTROL = 5
N_BOOT_MEASURED = 200


def measured_for(target: Target, orn: np.ndarray, orn_partners: np.ndarray, partner_names: Sequence[str],
                 n_boot: int = N_BOOT_MEASURED, seed: int = SEED, odors_all: Sequence[str] = MEGAMAT_ODORS) -> dict:
    """The measured step on a fit target's flies and odors (step 0's reading): its gains with the follows rule and its
    candidate partners. `orn` and `orn_partners` are the ORN recordings (flies x the target's glomeruli / the
    partner pool x `odors_all`), in the target's units."""
    oi = [list(odors_all).index(o) for o in target.odors]
    pn = fly_stack(target.flies.tables, target.glomeruli, target.odors)
    o = orn[:, :, oi]
    st = measured_step(o, pn, target.glomeruli, target.odors)
    iv, _ = step_intervals(o, pn, target.glomeruli, target.odors, n_boot, seed)
    st = st.join(iv)
    st['follows_provisional'] = follows_rule(st)
    res = residuals(o, pn, st, target.glomeruli, target.odors)
    om = pd.DataFrame(_nanmean(orn_partners[:, :, oi]), index=list(partner_names), columns=target.odors)
    cp = candidate_partners(res, om, n=N_CANDIDATES, n_perm=0)
    return dict(step=st, partners=cp)


def _key(t: Target) -> tuple:
    return (tuple(t.flies.flies), tuple(t.odors))


def run_rungs(stage: PNStage, pop, orn: np.ndarray, orn_partners: np.ndarray, partner_names: Sequence[str],
              lock: Mapping, directions: Sequence[Tuple[str, str]] = (('half_a', 'half_b'), ('half_b', 'half_a')),
              n_boot: int = N_BOOT_KEEP, seed: int = SEED, n_perm: int = N_PERM_CONTROL,
              max_lateral_steps: Optional[int] = None, global_maxiter: int = 400, verbose: bool = True) -> dict:
    """Step 1's rungs on the PN stage alone, in order, each cross-fitted (`cross_fit`: fitted fold by fold on one half
    of the PN flies, read on the other half on the odors it did not see).

    R0 the stage as it is; R1 the global parameters (`fit_global`); R2 with the gains, measured (R2m: the measured gain
    over the R1 stage's own, `measured_gain_start`, no fit) and fitted from there (R2f, `fit_gains`), R2 being R2f if it
    is kept against R2m (`keep_rung`), else R2m; R3 with lateral terms, forward-selected among each glomerulus's
    N_CANDIDATES best measured partners in each fold (a term is a glomerulus and a partner rank, so each fold uses its
    own partners; its weight fitted fold by fold, `fit_lateral`), added one at a time while the held-out error falls and
    the term gains more than the same term fitted to dendrites with odor labels permuted (n_perm permutations among
    the fold's training odors). Each rung is compared with the rung below by `keep_rung`; R1 and R2f are also refitted to
    odor-permuted dendrites (the controls, reported). `pop` is viki's PN population in the fit's units, its odors named
    as the stage's; `lock` gives 'halves', 'folds' and 'glomeruli'.

    Returns dict with 'rungs' (one row per rung: point and per-direction held-out matrix errors, the view errors, the SE
    and kept against the rung below), 'controls' (the permuted refits' points), 'selected' (R3's terms), 'folds' (per
    rung the folds' records), 'tables' (per rung the assembled held-out tables per direction), 'highest_kept'."""
    gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
    folds = [list(f) for f in lock['folds']]
    halves_ = dict(lock['halves'])
    odors = [o for f in folds for o in f]
    fit_t = {h: (lambda train, fl=fl: Target.make(flies_subset(pop, fl), train, gl)) for h, fl in halves_.items()}
    read_t = {h: Target.make(flies_subset(pop, fl), [o for o in pop.odors if o in odors], gl)
              for h, fl in halves_.items()}
    st: Dict[tuple, dict] = {}
    res: Dict[str, dict] = {}

    def say(msg):
        if verbose:
            print(f'{_now()} step 1: {msg}', flush=True)

    def cf(fit_one, targets=None):
        return cross_fit(fit_one, folds, targets or fit_t, read_t, directions, n_boot, seed)

    def f0(t, train):
        return stage.run(), {}

    def fitters(ns: Dict[tuple, dict]):
        """R1, R2m and R2f's fold fits, keeping each fold's state in `ns` (the real fits' or a control's own)."""
        def f1(t, train):
            r = fit_global(stage, t, maxiter=global_maxiter)
            ns.setdefault(_key(t), {})['R1'] = r['params']
            return stage.run(r['params']), dict(params=r['params'], loss=r['loss'], at_bounds=r['at_bounds'])

        def f2m(t, train):
            s_ = ns[_key(t)]
            meas = measured_for(t, orn, orn_partners, partner_names, seed=seed)
            start = measured_gain_start(meas['step'], model_gains(stage, s_['R1'], gl, train))
            s_.update(meas=meas, R2m=Terms(gains=start))
            return stage.run(s_['R1'], s_['R2m']), dict(gains=start)

        def f2f(t, train):
            s_ = ns[_key(t)]
            r = fit_gains(stage, t, s_['R1'], s_['R2m'].gains)
            s_['R2f'] = r['terms']
            return stage.run(s_['R1'], r['terms']), dict(gains=r['terms'].gains, loss=r['loss'],
                                                         at_bounds=r['at_bounds'])
        return f1, f2m, f2f

    f1, f2m, f2f = fitters(st)
    say('R0')
    res['R0'] = cf(f0)
    say('R1')
    res['R1'] = cf(f1)
    say('R2m')
    res['R2m'] = cf(f2m)
    say('R2f')
    res['R2f'] = cf(f2f)
    r2f_kept, _ = keep_rung(res['R2f']['draws'], res['R2m']['draws'], res['R2f']['point'], res['R2m']['point'])
    r2 = 'R2f' if r2f_kept else 'R2m'
    for k in st:
        st[k]['R2'] = st[k][r2]
        st[k]['R3'] = Terms(gains=dict(st[k]['R2'].gains), laterals=[])

    def partner(k, g, rank):
        cp = st[k]['meas']['partners']
        x = cp[(cp.glomerulus == g) & (cp['rank'] == rank)]
        return None if not len(x) else str(x.partner.iloc[0])

    def lateral_fit(term, targets_key_map=None, store=None):
        g, rank = term

        def f(t, train):
            k = _key(t) if targets_key_map is None else targets_key_map[_key(t)]
            s_ = st[k]
            p_ = partner(k, g, rank)
            if p_ is None or p_ not in stage.glomeruli:
                return stage.run(s_['R1'], s_['R3']), dict(partner=None)
            r = fit_lateral(stage, t, s_['R1'], s_['R3'], g, p_)
            if store is not None:
                store[k] = r['terms']
            return stage.run(s_['R1'], r['terms']), dict(partner=p_, weight=r['weight'], loss=r['loss'])
        return f

    current = res[r2]
    selected, excluded, steps = [], set(), 0
    cands = [(g, r) for g in gl for r in range(1, N_CANDIDATES + 1)]
    rng = np.random.default_rng(seed)
    while max_lateral_steps is None or steps < max_lateral_steps:
        n_of = {g: sum(1 for x in selected if x[0] == g) for g in gl}
        trial = [c for c in cands if c not in selected and c not in excluded and n_of[c[0]] < MAX_LATERAL]
        if not trial:
            break
        say(f'R3 step {steps + 1}: {len(trial)} candidate terms')
        best, best_out, best_store = None, None, None
        for c in trial:
            store: Dict[tuple, Terms] = {}
            out = cf(lateral_fit(c, store=store))
            if best_out is None or out['point'] < best_out['point']:
                best, best_out, best_store = c, out, store
        gain = current['point'] - best_out['point']
        if gain <= 0:
            break
        # the same term fitted to odor-permuted dendrites (the fold's training odors permuted)
        perm_gains = []
        for _ in range(n_perm):
            pseed = int(rng.integers(1 << 30))
            kmap: Dict[tuple, tuple] = {}

            def make_perm(h, fl=None):
                def mk(train, fl=halves_[h]):
                    base_t = Target.make(flies_subset(pop, fl), train, gl)
                    pt = Target.make(permute_odor_labels(flies_subset(pop, fl), np.random.default_rng(pseed), train),
                                     train, gl)
                    kmap[_key(pt)] = _key(base_t)
                    return pt
                return mk
            out_p = cf(lateral_fit(best, targets_key_map=kmap), {h: make_perm(h) for h in halves_})
            perm_gains.append(current['point'] - out_p['point'])
        if n_perm and gain <= float(np.mean(perm_gains)):
            excluded.add(best)
            say(f'R3: {best} gains {gain:.4f}, no more than odor-permuted ({np.mean(perm_gains):.4f}); left out')
            continue
        selected.append(best)
        for k, terms in best_store.items():
            st[k]['R3'] = terms
        current = best_out
        steps += 1
        say(f'R3: {best} added, held-out error {current["point"]:.4f} (gain {gain:.4f}; permuted '
            f'{np.mean(perm_gains) if perm_gains else float("nan"):.4f})')
    res['R3'] = current if selected else res[r2]

    # the controls: R1 and R2f refitted to odor-permuted dendrites
    controls = []
    for i in range(n_perm):
        pseed = int(rng.integers(1 << 30))
        pt = {h: (lambda train, fl=fl, pseed=pseed: Target.make(
            permute_odor_labels(flies_subset(pop, fl), np.random.default_rng(pseed), train), train, gl))
            for h, fl in halves_.items()}
        c1, c2m, c2f = fitters({})
        o1 = cf(c1, pt)
        o2m = cf(c2m, pt)
        o2f = cf(c2f, pt)
        controls.append(dict(perm=i, seed=pseed, R1=o1['point'], R2m=o2m['point'], R2f=o2f['point']))

    order = ['R0', 'R1', 'R2m', 'R2f', 'R3']
    below = dict(R1='R0', R2m='R1', R2f='R1', R3=r2)
    rows = []
    for r in order:
        o = res[r]
        row = dict(rung=r, point=o['point'])
        for d in [k for k in o if '->' in k]:
            row[f'point_{d}'] = o[d]['point']
            sc = o[d]['score']
            row.update({f'{d}_{k}': v for k, v in sc.as_series().items() if k.startswith(('error', 'view_'))})
        if r in below:
            kept, se = keep_rung(o['draws'], res[below[r]]['draws'], o['point'], res[below[r]]['point'])
            row.update(below=below[r], se_vs_below=se, kept=kept)
        rows.append(row)
    rungs = pd.DataFrame(rows)
    kept = rungs.set_index('rung')['kept'].to_dict()
    highest = 'R0'
    for r in ('R1', r2, 'R3'):
        if kept.get(r) and (r != 'R3' or selected):
            highest = r
    return dict(rungs=rungs, controls=pd.DataFrame(controls), selected=[list(x) for x in selected],
                excluded=[list(x) for x in excluded], r2=r2, highest_kept=highest,
                folds={r: {d: res[r][d]['folds'] for d in res[r] if '->' in d} for r in order},
                tables={r: {d: res[r][d]['table'] for d in res[r] if '->' in d} for r in order},
                state={str(k): dict(R1=v.get('R1'), R2=v['R2'].as_dict() if 'R2' in v else None,
                                    R3=v['R3'].as_dict() if 'R3' in v else None) for k, v in st.items()})


def part_step1(n_perm: int = N_PERMUTATIONS) -> dict:
    """Step 1: the rungs fitted in order to the dendrites and read on the held-out flies and odors, with the
    odor-permuted controls and the kept rung; refuses before E134's gate and this experiment's committed lock.

    Not run yet: it waits on E134's gate, the lock's folds, and the family's olfsysm parameters for the PN stage
    (E83's `build_simulator(...)` info['cells']['mp']), which a model build reads."""
    require_e134_gate()
    lock = read_lock()
    t0 = time.time()
    stage, info = family_stage('megamat')
    ident = info['identity']
    units = lock.get('fit_units', FIT_UNITS)
    pop = pn_population(units)
    orn_t = recorded_stage(ORN_KEY, units)
    gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
    model_gl = list(H.apl16_orn_deltas().index)
    pool = [g for g in model_gl if g in set().union(*[set(t.index) for t in orn_t.values()]) and g in stage.glomeruli]
    out = run_rungs(stage, pop, fly_stack(orn_t, gl, MEGAMAT_ODORS), fly_stack(orn_t, pool, MEGAMAT_ODORS), pool,
                    lock, n_perm=n_perm)
    d = out_dir()
    out['rungs'].to_csv(d / 'step1_rungs.csv', index=False)
    out['controls'].to_csv(d / 'step1_controls.csv', index=False)
    tabs = []
    for r, by in out['tables'].items():
        for direction, t in by.items():
            x = t.stack().rename('value').reset_index()
            x.columns = ['glomerulus', 'odor', 'value']
            tabs.append(x.assign(rung=r, direction=direction))
    pd.concat(tabs, ignore_index=True).to_csv(d / 'step1_tables.csv', index=False)
    fits = rung_fits_all_flies(stage, pop, fly_stack(orn_t, gl, MEGAMAT_ODORS), fly_stack(orn_t, pool, MEGAMAT_ODORS),
                               pool, lock, out)
    _write_json(d / 'step1_fits.json', dict(folds=out['folds'], state=out['state'], selected=out['selected'],
                                            excluded=out['excluded'], r2=out['r2'], highest_kept=out['highest_kept'],
                                            all_flies=fits))
    _write_json(d / 'provenance_step1.json', dict(script=H.script_path(__file__), part='step1', host=_host(),
                time=_now(), seconds=time.time() - t0, lock_sha256=_sha256(LOCK_COPY), units=units,
                identity=ident, n_perm=n_perm, **H.model_provenance(orn_deltas=stage.delta)))
    return out


def family_stage(panel: str = 'megamat', wiring: str = 'real') -> Tuple[PNStage, dict]:
    """The family's PN stage on a panel's input ('megamat': E107's `panel_input`; 'apl16': the lab's input): E83's
    `build_simulator(wiring, orn=..., breadth=<the lab input's>)`, its olfsysm variables read into a `PNStage` (rows
    the wiring's glomeruli, columns the input's odors abbreviated, in olfsysm's order through the simulator's
    'odor_order'), checked against olfsysm's own layers (`check_identity`) on a copy of the build's variables. A model
    build: after E134's gate only."""
    require_e134_gate()
    from experiments import E83_cross_wiring_fit as X
    orn = panel_orn(panel)
    sim, info = X.build_simulator(wiring, orn=orn, breadth=X.wiring_parts(wiring)['breadth'])
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    names = [None] * len(order)
    for j, i in enumerate(order):
        names[int(i)] = H.odor_abbrev(orn.columns[j])
    gl = [str(g) for g in dict.fromkeys(info['W'].columns.get_level_values('glomerulus'))]
    stage = PNStage.from_mp(mp, gl, names)
    want = orn.copy()
    want.columns = [H.odor_abbrev(c) for c in want.columns]
    want = want.reindex(index=gl, columns=names).fillna(0.0)
    if not np.allclose(stage.delta.to_numpy(), want.to_numpy(), atol=1e-12):
        raise RuntimeError(f"{wiring} {panel}: olfsysm's delta is not the panel's input on the wiring's glomeruli")
    ident = stage.check_identity(mp, rv)
    return stage, dict(sim=sim, info=info, mp=mp, rv=rv, identity=ident, orn=orn)


def panel_orn(panel: str) -> pd.DataFrame:
    """The model's ORN input of a panel: 'megamat' E107's `panel_input('megamat')` (17 odors on the lab input's 36
    glomeruli), 'apl16' the lab's input (`apl16_orn_deltas`)."""
    if panel == 'megamat':
        from experiments import E107_validation2_panel as E107
        return E107.panel_input('megamat')
    if panel == 'apl16':
        return H.apl16_orn_deltas()
    raise ValueError(f'{panel=}')


def rung_fits_all_flies(stage: PNStage, pop, orn: np.ndarray, orn_partners: np.ndarray,
                        partner_names: Sequence[str], lock: Mapping, rungs_out: Mapping) -> dict:
    """The rungs refitted on all PN flies for step 2: fold by fold (each odor's table from the fit that left it out)
    and once on all odors (for the APL16 input's odors outside the folds); R1's parameters, R2's gains (measured or
    fitted, as step 1 chose), R3's selected terms (each fold's own partner of each selected glomerulus and rank).
    PROVISIONAL: which flies step 2's tables are fitted to is a planner's question (see the report)."""
    gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
    folds = [list(f) for f in lock['folds']]
    odors = [o for f in folds for o in f]
    sets = [(k, [o for o in odors if o not in set(f)]) for k, f in enumerate(folds)] + [('all', odors)]
    out = {}
    for k, train in sets:
        t = Target.make(pop, train, gl)
        r1 = fit_global(stage, t)['params']
        meas = measured_for(t, orn, orn_partners, partner_names)
        start = measured_gain_start(meas['step'], model_gains(stage, r1, gl, train))
        terms = Terms(gains=start)
        if rungs_out['r2'] == 'R2f':
            terms = fit_gains(stage, t, r1, start)['terms']
        r2 = terms
        for g, rank in rungs_out['selected']:
            cp = meas['partners']
            x = cp[(cp.glomerulus == g) & (cp['rank'] == rank)]
            if len(x) and str(x.partner.iloc[0]) in stage.glomeruli:
                terms = fit_lateral(stage, t, r1, terms, g, str(x.partner.iloc[0]))['terms']
        out[str(k)] = dict(held_out=[o for o in odors if o not in set(train)], R1=r1, R2=r2.as_dict(),
                           R3=terms.as_dict())
    return out


def rung_table(stage: PNStage, fits: Mapping[str, dict], rung: str, reading: str = 'preset') -> pd.DataFrame:
    """A rung's table of PN rates on a stage's input for step 2 (`rung_fits_all_flies`' fits): each odor a fold held
    out from that fold's fit, every other odor (the APL16 input's paa, euc, -aPine, B-myr) from the fit on all odors;
    rates in `reading`."""
    held = {o: k for k, f in fits.items() if k != 'all' for o in f['held_out']}
    cache: Dict[str, pd.DataFrame] = {}
    cols = {}
    for o in stage.odors:
        k = held.get(o, 'all')
        if k not in fits:
            continue
        if k not in cache:
            f = fits[k]
            cache[k] = stage.run(f['R1'], None if rung == 'R1' else Terms.from_dict(f[rung]), reading=reading)
        cols[o] = cache[k][o]
    return pd.DataFrame(cols)[[o for o in stage.odors if o in cols]]


# ---------------------------------------------------------------------------
# the imposed dendrites, both routes

def rank_assign(values: pd.DataFrame, reference: np.ndarray) -> pd.DataFrame:
    """Rank for rank: where `values` and `reference` have as many finite entries, the entry of rank r gets the
    reference's r-th smallest value (ties in `values` share their values' mean), so the result holds exactly the
    reference's values in `values`' arrangement; otherwise hal's `rank_map` (quantile matching)."""
    ref = np.sort(np.asarray(reference, dtype=float)[np.isfinite(reference)])
    v = values.to_numpy(dtype=float)
    ok = np.isfinite(v)
    if ok.sum() != len(ref):
        return P.rank_map(values, ref)
    flat = v[ok]
    order = np.argsort(flat, kind='stable')
    out_flat = np.empty_like(flat)
    out_flat[order] = ref
    ranks = stats.rankdata(flat)
    for r in np.unique(ranks[ranks != np.round(ranks)]):
        tie = ranks == r
        out_flat[tie] = out_flat[tie].mean()
    out = np.full(v.shape, np.nan)
    out[ok] = out_flat
    return pd.DataFrame(out, index=values.index, columns=values.columns)


OUT_OF_BOUNDS_POLICIES: Tuple[Optional[str], ...] = (None, 'clip')
BOUNDS_MARGIN_HZ = 1e-6


def imposed_rates(model_rates: pd.DataFrame, dendrites: pd.DataFrame, glomeruli: Sequence[str],
                  odors: Optional[Sequence[str]] = None,
                  bounds: Optional[Tuple[pd.DataFrame, pd.DataFrame]] = None,
                  out_of_bounds: Optional[str] = None) -> Tuple[pd.DataFrame, dict]:
    """The imposed table: on the shared glomeruli and `odors` (default the dendrites' columns) the dendrites'
    arrangement mapped rank for rank onto the model's own pooled PN rates there (`rank_assign`: exactly the model's values,
    in the dendrites' order); every other entry keeps the model's.

    A mapped rate can lie outside what that glomerulus's PN can reach for that odor (below its floor, about minus its
    spontaneous rate, since the pooled rates reach the most inhibited glomerulus's floor): with `bounds` (lo, hi, as
    `pn_rate_bounds` or `PNStage.rate_bounds` give them) such entries are counted, and `out_of_bounds` says what is
    done: None (the default) raises ValueError, 'clip' moves them BOUNDS_MARGIN_HZ inside the bounds. Which is right is a
    planner's question (the entry does not settle it).

    Returns (the table, a record: entries mapped, entries outside the bounds, the largest excursion, the policy)."""
    if out_of_bounds not in OUT_OF_BOUNDS_POLICIES:
        raise ValueError(f'{out_of_bounds=}')
    odors = list(dendrites.columns) if odors is None else list(odors)
    gl = list(glomeruli)
    ref = model_rates.loc[gl, odors].to_numpy(dtype=float).ravel()
    mapped = rank_assign(dendrites.loc[gl, odors].astype(float), ref)
    rec = dict(n_mapped=int(np.isfinite(mapped.to_numpy()).sum()), n_outside=0, largest_excursion_hz=0.0,
               policy=out_of_bounds)
    if bounds is not None:
        lo = bounds[0].loc[gl, odors].to_numpy(dtype=float)
        hi = bounds[1].loc[gl, odors].to_numpy(dtype=float)
        m = mapped.to_numpy(dtype=float)
        below, above = m <= lo, m >= hi
        rec['n_outside'] = int(below.sum() + above.sum())
        rec['largest_excursion_hz'] = float(max(np.max(np.where(below, lo - m, 0.0), initial=0.0),
                                                np.max(np.where(above, m - hi, 0.0), initial=0.0)))
        rec['outside'] = [dict(glomerulus=gl[i], odor=odors[j], rate=float(m[i, j]), lo=float(lo[i, j]),
                               hi=float(hi[i, j])) for i, j in zip(*np.nonzero(below | above))]
        if rec['n_outside'] and out_of_bounds is None:
            raise ValueError(f"{rec['n_outside']} imposed rates outside the PN stage's reach (largest by "
                             f"{rec['largest_excursion_hz']:.2f} Hz); choose a policy (out_of_bounds='clip' or "
                             'another map)')
        if rec['n_outside']:
            m = np.clip(m, lo + BOUNDS_MARGIN_HZ, hi - BOUNDS_MARGIN_HZ)
            mapped = pd.DataFrame(m, index=mapped.index, columns=mapped.columns)
    out = model_rates.copy().astype(float)
    out.loc[gl, odors] = mapped.to_numpy()
    return out, rec


def solve_input(stage: PNStage, target: pd.DataFrame, glomeruli: Optional[Sequence[str]] = None,
                params: Optional[Mapping[str, float]] = None, tol: float = SOLVE_TOL_HZ, max_outer: int = 30,
                n_iter: int = 60, bracket: Tuple[float, float] = (-400.0, 2000.0),
                reading: str = 'preset', delta_tol: float = 1e-3) -> Tuple[pd.DataFrame, dict]:
    """Route 2: the ORN input solved glomerulus by glomerulus so that the unchanged PN stage returns `target`.

    The deltas of `glomeruli` (default the target's rows: the input's glomeruli) are solved, every other row of the
    stage's input kept. The ORN layer is linear in the delta, a PN's evoked rate rises monotonically with its own
    glomerulus's delta at a given LN, and the LN reads only the mean ORN rate over glomeruli: so for each odor the LN
    is fixed by one number, the mean delta over glomeruli. The solve finds, odor by odor, the mean at which the deltas
    solved at that mean's LN (every entry at once, by the Illinois method inside `bracket`) have that mean (a secant
    iteration on the mean), and then checks the stage, re-run exactly with its LN on the solved input, against the
    target: to within `tol` Hz on every solved entry. Rates are in `reading` (default 'preset', the reading of route 1's
    `pn_rates`, so that both routes aim at the same table).

    Returns:
        (the solved input, glomeruli x odors as the stage's delta; the solve's record: outer (secant) iterations, the
        largest difference from the target, entries at the bracket's edge (a target outside the stage's reach), the last
        change of the mean delta, and the largest change of the other glomeruli's rates).
    """
    gl = list(target.index if glomeruli is None else glomeruli)
    od = [o for o in stage.odors if o in target.columns]
    if not od:
        raise ValueError("the target has none of the stage's odors")
    sp = stage.sp.with_global(params)
    G = stage.glomeruli
    rows = np.array([G.index(g) for g in gl])
    cols = np.array([stage.odors.index(o) for o in od])
    T = target.loc[gl, od].to_numpy(dtype=float)
    d = stage.delta.to_numpy(dtype=float).copy()
    base = stage.run(params, reading=reading)
    others = [i for i in range(len(G)) if i not in set(rows)]
    k = orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
    fixed = d[others][:, cols].sum(axis=0) if others else np.zeros(len(cols))
    nG = len(G)
    inhA_full, inhB_full = stage.ln(sp)
    inhA_full, inhB_full = inhA_full.copy(), inhB_full.copy()
    edge = dict(n=0)

    def solved_at(m):
        mean_t = stage.spont.mean() + m[:, None] * k[None, :]
        a, b = ln_layer_mean(sp, mean_t, nG)
        inhA_full[cols], inhB_full[cols] = a, b
        x, n_edge = _solve_entries(sp, stage.spont, rows, cols, T, k, inhA_full, inhB_full, reading, bracket, n_iter,
                                   tol / 100)
        edge['n'] = n_edge
        return (x.sum(axis=0) + fixed) / nG - m, x

    # the mean delta of each odor at which the solved deltas have that mean: F(m) = mean(solved(m)) - m, bracketed
    # outward from the input's own mean, then the Illinois method, every odor at once
    thr = delta_tol / nG
    m0 = d[:, cols].mean(axis=0)
    F0, x0 = solved_at(m0)
    it = 1
    done = np.abs(F0) <= thr
    lo_m, flo, hi_m, fhi = m0.copy(), F0.copy(), m0.copy(), F0.copy()
    step = np.where(np.abs(F0) > 1e-3, F0, np.sign(F0) * 1e-3)
    found = done.copy()
    # (outward on both sides: F falls with m where the solved mean rises more slowly than m, and rises where faster)
    for _ in range(max_outer):
        if found.all():
            break
        for sgn in (1.0, -1.0):
            m1 = np.where(found, m0, m0 + sgn * step)
            F1, _ = solved_at(m1)
            it += 1
            cross = ~found & (np.sign(F1) != np.sign(F0))
            lo_m, flo = np.where(cross, m0, lo_m), np.where(cross, F0, flo)
            hi_m, fhi = np.where(cross, m1, hi_m), np.where(cross, F1, fhi)
            found = found | cross
            if found.all():
                break
        step = np.where(found, step, 2.0 * step)
    x, moved = x0, float(np.abs(F0).max())
    if not done.all():
        side = np.zeros(len(cols))
        for _ in range(max_outer):
            den = fhi - flo
            m = np.where(done, m0, np.where(den != 0, hi_m - fhi * (hi_m - lo_m) / np.where(den != 0, den, 1.0),
                                            0.5 * (lo_m + hi_m)))
            Fm, x = solved_at(m)
            it += 1
            moved = float(np.abs(Fm[~done]).max())
            if moved <= thr:
                break
            same = (np.sign(Fm) == np.sign(flo)) & ~done
            other = ~same & ~done
            lo_m, flo = np.where(same, m, lo_m), np.where(same, Fm, flo)
            hi_m, fhi = np.where(other, m, hi_m), np.where(other, Fm, fhi)
            fhi = np.where(same & (side == -1), fhi / 2, fhi)
            flo = np.where(other & (side == 1), flo / 2, flo)
            side = np.where(same, -1.0, np.where(other, 1.0, side))
    def check(xs):
        dd = d.copy()
        dd[np.ix_(rows, cols)] = xs
        got = PNStage(pd.DataFrame(dd, index=G, columns=stage.odors), pd.Series(stage.spont, index=G), sp).run(
            reading=reading)
        return dd, got, float(np.nanmax(np.abs(got.loc[gl, od].to_numpy() - T)))

    d_best, got, diff = check(x)
    edges, how = int(edge['n']), 'bracketed'
    if diff > tol or edges:
        # no self-consistent mean found cleanly (strong normalization can make the rates nearly blind to the mean):
        # plain alternation from the input's own mean, keeping the deltas whose exact rates are nearest the target
        m = m0.copy()
        for _ in range(max_outer):
            Fm, xm = solved_at(m)
            it += 1
            dd, gg, df = check(xm)
            if edge['n'] == 0 and (df < diff or edges):
                d_best, got, diff, edges, how = dd, gg, df, 0, 'alternation'
            if diff <= tol:
                break
            m = m + Fm
    rec = dict(outer=it, max_abs_diff=diff, at_edge=edges, tol=tol, last_mean_move=moved, method=how,
               others_max_change=float(np.abs(got.iloc[others][od].to_numpy() - base.iloc[others][od].to_numpy()
                                              ).max()) if others else 0.0)
    rec['converged'] = bool(rec['max_abs_diff'] <= tol and rec['at_edge'] == 0)
    rec['fixed_point'] = bool(moved <= thr)
    rec['bracketed'] = bool(found.all())
    return pd.DataFrame(d_best, index=G, columns=stage.odors), rec


def _solve_entries(sp: StageParams, spont: np.ndarray, rows: np.ndarray, cols: np.ndarray, T: np.ndarray,
                   kernel: np.ndarray, inhA: np.ndarray, inhB: np.ndarray, reading: str,
                   bracket: Tuple[float, float], n_iter: int, ftol: float) -> Tuple[np.ndarray, int]:
    """Every entry's own delta at which its PN's evoked rate is the target, at a given LN, by the Illinois method
    (regula falsi with the retained end halved) inside `bracket`, all entries at once; entries whose target lies outside
    the rates the bracket reaches are put at its edge and counted."""
    f = lambda x: _rates_at(sp, spont, rows, cols, x, kernel, inhA, inhB, reading) - T
    lo = np.full(T.shape, float(bracket[0]))
    hi = np.full(T.shape, float(bracket[1]))
    flo, fhi = f(lo), f(hi)
    below, above = flo >= 0, fhi <= 0
    inside = ~(below | above)
    side = np.zeros(T.shape)
    x = 0.5 * (lo + hi)
    for _ in range(n_iter):
        den = fhi - flo
        x = np.where(den > 0, hi - fhi * (hi - lo) / np.where(den > 0, den, 1.0), 0.5 * (lo + hi))
        x = np.clip(x, lo, hi)
        fx = f(x)
        neg = fx < 0
        lo, flo = np.where(neg, x, lo), np.where(neg, fx, flo)
        hi, fhi = np.where(neg, hi, x), np.where(neg, fhi, fx)
        fhi = np.where(neg & (side == -1), fhi / 2, fhi)
        flo = np.where(~neg & (side == 1), flo / 2, flo)
        side = np.where(neg, -1.0, 1.0)
        if not inside.any() or float(np.nanmax(np.abs(fx[inside]))) <= ftol:
            break
    x = np.where(below, bracket[0], np.where(above, bracket[1], x))
    return x, int(below.sum() + above.sum())


def _rates_at(sp: StageParams, spont: np.ndarray, rows: np.ndarray, cols: np.ndarray, deltas: np.ndarray,
              kernel: np.ndarray, inhA: np.ndarray, inhB: np.ndarray, reading: str = 'stage') -> np.ndarray:
    """Evoked PN rates of glomeruli `rows` at odors `cols` with their own ORN deltas set to `deltas` (rows x cols), at
    the given LN (odors x steps): each entry a PN reading delta x kernel(t) as its ORN delta over time, its rate in
    `reading`."""
    win, pre = windows(sp, reading)
    R_, C_ = deltas.shape
    s = spont[rows][:, None] * np.ones((1, C_))
    spont_pn = s * sp.pn_inhsc / (spont.sum() + sp.pn_inhadd)
    inh_t = sp.pn_inhsc / (sp.pn_inhadd + 0.25 * inhA[cols] + 0.75 * inhB[cols])        # cols x steps
    pn = s.copy()
    acc_w = pn * win[0]
    acc_p = pn * pre[0]
    inh = np.zeros((1, C_))
    for t in range(1, len(kernel)):
        drive = deltas * kernel[t - 1]
        dd = -pn + spont_pn + 200.0 * np.tanh((drive + sp.pn_offset) * sp.pn_tanhsc / 200.0 * inh)
        inh = inh_t[:, t][None, :]
        pn = pn + dd * sp.dt / sp.pn_taum
        pn = np.where(0.0 < pn, pn, 0.0)
        if win[t]:
            acc_w += pn
        if pre[t]:
            acc_p += pn
    return acc_w / win.sum() - acc_p / pre.sum()


def solve_identity(stage: PNStage, glomeruli: Optional[Sequence[str]] = None, tol: float = SOLVE_TOL_HZ,
                   reading: str = 'preset') -> dict:
    """The solve's identity check: the stage's own rates as the target give back the stage's input (where the rate
    depends on the delta) and the same rates."""
    gl = list(stage.glomeruli if glomeruli is None else glomeruli)
    own = stage.run(reading=reading)
    solved, rec = solve_input(stage, own.loc[gl], gl, tol=tol / 10, reading=reading)
    d0 = stage.delta.loc[gl].to_numpy(dtype=float)
    d1 = solved.loc[gl].to_numpy(dtype=float)
    back = PNStage(solved, pd.Series(stage.spont, index=stage.glomeruli), stage.sp).run(reading=reading)
    rate_diff = float(np.abs(back.loc[gl].to_numpy() - own.loc[gl].to_numpy()).max())
    ok_rows = own.loc[gl].to_numpy() > 0
    res = dict(solve=rec, max_rate_diff=rate_diff, max_delta_diff=float(np.abs(d1 - d0)[ok_rows].max())
               if ok_rows.any() else 0.0, tol=tol)
    res['ok'] = bool(rate_diff <= tol and rec['converged'])
    return res


def route1_pn_rates(table: pd.DataFrame, set_glomeruli: Sequence[str],
                    set_odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Route 1, olfsysm's option for the PN stage's rates given directly (viki's, olfsysm 7528f6d / b1e89d7, Yang's OK
    of 2026-11-24 09:43): fit_mb_model's `pn_rates`, a glomeruli x odors table in the reading of `pn_evoked_rates`
    ('preset' here), NaN keeping the model's PN. Returns `table` on `set_glomeruli` and `set_odors` (default all its
    columns), NaN elsewhere, to pass as `extra_kws=dict(pn_rates=...)` to E83's `build_simulator` (or to
    fit_mb_model). The odor's LN inhibition stays the model's (route 2's solved input changes it)."""
    odors = list(table.columns if set_odors is None else set_odors)
    out = pd.DataFrame(np.nan, index=list(table.index), columns=list(table.columns))
    gl = [g for g in table.index if g in set(set_glomeruli)]
    out.loc[gl, odors] = table.loc[gl, odors].to_numpy(dtype=float)
    if not np.isfinite(out.to_numpy()).any():
        raise ValueError('no rate set')
    if np.isinf(out.to_numpy()).any():
        raise ValueError('infinite rates')
    return out


def route1_check(mp, rv, wPNKC: pd.DataFrame, odors: Sequence, table: pd.DataFrame, tol: float = 1e-6) -> dict:
    """A route-1 run's identity: its PN rates (pyolfsysm's `pn_evoked_rates` of the run's `mp`, `rv`) equal the set
    entries of `table` to `tol` Hz, and the unset entries are the model's own (checked by the caller against the run
    without `pn_rates`). Raises RuntimeError beyond `tol`."""
    from pyolfsysm.mb_model import pn_evoked_rates
    got = pn_evoked_rates(mp, rv, wPNKC=wPNKC, odors=odors)
    t = table.reindex(index=got.index, columns=got.columns)
    m = np.isfinite(t.to_numpy())
    d = float(np.abs(got.to_numpy()[m] - t.to_numpy()[m]).max()) if m.any() else 0.0
    if d > tol:
        raise RuntimeError(f'route 1: the run\'s PN rates differ from the set table by {d} Hz')
    return dict(n_set=int(m.sum()), max_abs_diff=d, tol=tol)


# ---------------------------------------------------------------------------
# step 2: the runs, and the score

def planned_step2(reverse: bool = False, with_rung: bool = True) -> List[Tuple[Arm, str]]:
    """Step 2's points (arm, wiring): 48 on the megamat input, 24 on the APL16 input, one repeat by the second route
    (73); without the rung's 24 if none is kept."""
    pts = [(Arm('megamat', a), w) for a in MEGAMAT_ARMS for w in WIRINGS]
    pts += [(Arm('apl16', a), w) for a in APL16_ARMS for w in WIRINGS]
    pts += [(REPEAT_ARM, WIRINGS[0])]
    if not with_rung:
        pts = [(a, w) for a, w in pts if a.arm not in RUNG_ARMS]
    return pts[::-1] if reverse else pts


def tag_of(arm: Arm, wiring: str) -> str:
    """'<panel>_<arm>_<wiring>', e.g. 'megamat_imposed_glom_perm_P1'."""
    return f'{arm.panel}_{arm.arm}_{wiring}'


def claim(tag: str, lane: str, reclaim: bool = False, root: Optional[Path] = None) -> bool:
    """Claims a point by mkdir of results/E135_antennal_lobe_dendrites/claims/<tag>/ (atomic on NFS), as E134's runner
    does (its `claim`, rooted here)."""
    from experiments import E134_population_baseline as E134R
    return E134R.claim(tag, lane, reclaim=reclaim, root=out_dir() if root is None else root)


def release(tag: str, root: Optional[Path] = None) -> None:
    from experiments import E134_population_baseline as E134R
    E134R.release(tag, root=out_dir() if root is None else root)


def arm_table(arm: Arm, stage: PNStage, dendrites: Optional[pd.DataFrame], glomeruli: Sequence[str],
              fits: Optional[Mapping[str, dict]] = None, rung: Optional[str] = None, seed: int = SEED,
              out_of_bounds: Optional[str] = None) -> Tuple[pd.DataFrame, List[str], dict]:
    """An arm's table of PN rates on a stage's input ('preset' reading) and the glomeruli it sets.

    'R1': the R1 fits' stage (`rung_table`); 'rung' / 'rung_glom_perm': the highest kept rung (its terms' glomerulus
    labels permuted, `permute_terms`, for the control); 'imposed' and its permutations: the dendrites (fly mean of all
    PN flies, glomeruli x odors) on the shared glomeruli, glomerulus or odor labels permuted for the controls, rank for
    rank on the stage's own rates (`imposed_rates`, with the stage's reachable bounds and the `out_of_bounds` policy)."""
    rng = np.random.default_rng(seed)
    gl = [g for g in glomeruli if g in stage.glomeruli]
    rec: dict = dict(arm=arm.arm, panel=arm.panel)
    if arm.arm == 'R1':
        return rung_table(stage, fits, 'R1'), list(stage.glomeruli), rec
    if arm.arm in ('rung', 'rung_glom_perm'):
        if arm.arm == 'rung':
            return rung_table(stage, fits, rung), list(stage.glomeruli), rec
        pf = {k: dict(f, **{rung: permute_terms(Terms.from_dict(f[rung]), gl, np.random.default_rng(seed)).as_dict()})
              for k, f in fits.items()}
        return rung_table(stage, pf, rung), list(stage.glomeruli), rec
    own = stage.run(reading='preset')
    den = dendrites.loc[gl, [o for o in dendrites.columns if o in stage.odors]].astype(float)
    if arm.arm in ('imposed_glom_perm',):
        den = permute_glomerulus_labels(den, gl, rng)
    elif arm.arm == 'imposed_odor_perm':
        cols = list(den.columns)
        perm = [cols[i] for i in rng.permutation(len(cols))]
        den = pd.DataFrame(den[perm].to_numpy(), index=den.index, columns=cols)
    elif arm.arm not in ('imposed', 'imposed_route2_repeat'):
        raise ValueError(arm)
    tab, r = imposed_rates(own, den, gl, bounds=stage.rate_bounds(), out_of_bounds=out_of_bounds)
    rec.update(r)
    return tab, gl, rec


def arm_input(stage: PNStage, table: pd.DataFrame, set_glomeruli: Sequence[str], route: int) -> Tuple[dict, dict]:
    """How an arm's table enters the family's simulator: route 1, extra_kws=dict(pn_rates=<the set entries>) on the
    panel's input (`route1_pn_rates`); route 2, the input solved on the stage (`solve_input`: the set glomeruli at the
    table's rates, the other input glomeruli at their own rates there, the LN on the solved input), as orn=. Returns
    (build_simulator's keyword arguments, the record)."""
    if route == 1:
        return dict(extra_kws=dict(pn_rates=route1_pn_rates(table, set_glomeruli))), dict(route=1)
    if route != 2:
        raise ValueError(f'{route=}')
    own = stage.run(reading='preset')
    target = own.copy()
    target.loc[list(set_glomeruli)] = table.loc[list(set_glomeruli), own.columns].to_numpy()
    rows = [g for g in stage.glomeruli if np.any(stage.delta.loc[g].to_numpy() != 0) or g in set(set_glomeruli)]
    solved, rec = solve_input(stage, target.loc[rows], rows)
    if not rec['converged']:
        raise RuntimeError(f'route 2: the solve did not converge: {rec}')
    return dict(orn=solved), dict(route=2, solve=rec)


def part_step2(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 2, resumable, one point claimed at a time; refuses before E134's gate. A point: its arm's table of PN rates
    (`arm_table`: step 1's kept rung or R1 from `rung_fits_all_flies`; the imposed dendrites and their permutations),
    entered by route 1 (fit_mb_model's pn_rates) or route 2 (the solved input) (`arm_input`), the family's simulator on
    it (`family_stage`'s E83 `build_simulator`, breadth held at the lab input's), both response rates re-held, every
    stage saved and read as E134's step 3 (atlas's runner: its `held_targets`, `run_stages`, `check_identity`,
    `save_run`, still stubs at E135's commit, which take a point of their own panels and need an ORN-input and
    extra-kwargs override, a request to atlas)."""
    require_e134_gate()
    raise NotImplementedError("step 2's runs wait on E134's gate and on atlas's step-3 runner (an ORN-input and "
                              "extra-kwargs override of its build_point and run_stages); E135's points, claims, arm "
                              'tables and both routes are implemented')


def class_thresholds(stage_table: Optional[pd.DataFrame] = None) -> dict:
    """The classes' thresholds reset to THRESHOLD_SE_MULT standard errors over wirings in E134's stage table: per
    megamat stage, 3 x the SD over its eight wirings' population errors / sqrt(8) (the family as fitted), in place of
    the entry's 0.15 (the somata's fall) and 0.05 (the somata within); PROVISIONAL (which thresholds are reset is a
    planner's question). Reads <shared>/results/E134_population_baseline/score_stage_table.csv by default."""
    if stage_table is None:
        p = H.shared_dir() / 'results' / E134 / E134_STAGE_TABLE
        if not p.exists():
            raise FileNotFoundError(f'{p}: E134\'s stage table is not written yet')
        stage_table = pd.read_csv(p)
    t = stage_table
    need = {'stage', 'wiring', 'error'}
    if not need <= set(t.columns):
        raise KeyError(f'the stage table needs columns {sorted(need)}')
    out = {}
    for st, sub in t[t.wiring.isin(WIRINGS)].groupby('stage'):
        e = sub.error.to_numpy(dtype=float)
        se = float(np.nanstd(e, ddof=1) / np.sqrt(np.isfinite(e).sum()))
        out[st] = dict(se=se, threshold=THRESHOLD_SE_MULT * se, n_wirings=int(np.isfinite(e).sum()))
    return out


def classify(scores: pd.DataFrame, thresholds: Optional[Mapping[str, Mapping[str, float]]] = None) -> dict:
    """The entry's classes from step 2's scores against E134's runs.

    Args:
        scores: one row per (arm, wiring, stage) of the megamat input with columns arm, wiring, stage
            ('pn_bouton' / 'kc_soma'), error (the population error), error_base (E134's run of the same wiring),
            strength (the view's error) and strength_base, and excluded (no class: outside a band on the APL16 input,
            or more than 1 % of bouton responses above 300 Hz).
        thresholds: per stage {'threshold': ...} (`class_thresholds`); default the entry's fixed numbers.

    Returns:
        dict with 'upstream', 'boutons_only', 'not_upstream' (bools), 'global_carry', 'terms_carry' and the per-arm
        counts behind them.
    """
    s = scores[~scores.get('excluded', pd.Series(False, index=scores.index)).astype(bool)]
    thr_soma = SOMA_DROP_UPSTREAM if thresholds is None else float(thresholds['kc_soma']['threshold'])
    thr_within = SOMA_WITHIN_BOUTONS_ONLY if thresholds is None else float(thresholds['kc_soma']['threshold'])

    def drop(arm, stage, col='error'):
        x = s[(s.arm == arm) & (s.stage == stage)].set_index('wiring')
        return x[f'{col}_base'] - x[col]

    def n_ok(mask):
        return int(np.sum(np.asarray(mask, dtype=bool)))

    out = {}
    soma = drop('imposed', 'kc_soma')
    soma_g = drop('imposed_glom_perm', 'kc_soma').reindex(soma.index)
    soma_o = drop('imposed_odor_perm', 'kc_soma').reindex(soma.index)
    b = s[(s.arm == 'imposed') & (s.stage == 'pn_bouton')].set_index('wiring')
    halved = (b.strength <= BOUTON_STRENGTH_HALVED * b.strength_base) if len(b) else pd.Series(dtype=bool)
    n_soma = n_ok((soma >= thr_soma) & (soma > soma_g) & (soma > soma_o))
    n_halved = n_ok(halved)
    out['upstream'] = bool(n_soma >= MIN_WIRINGS and n_halved >= MIN_WIRINGS)
    n_within = n_ok(soma.abs() <= thr_within)
    out['boutons_only'] = bool(not out['upstream'] and n_halved >= MIN_WIRINGS and n_within >= MIN_WIRINGS)
    out['not_upstream'] = bool(not out['upstream'] and not out['boutons_only'])
    gain_b = drop('imposed', 'pn_bouton')
    r1_b, r1_s = drop('R1', 'pn_bouton'), drop('R1', 'kc_soma')
    rg_b, rg_s = drop('rung', 'pn_bouton'), drop('rung', 'kc_soma')
    rp_b, rp_s = drop('rung_glom_perm', 'pn_bouton'), drop('rung_glom_perm', 'kc_soma')
    half = lambda x, ref: (x >= CARRY_SHARE * ref.reindex(x.index))
    n_r1 = n_ok(half(r1_b, gain_b) & half(r1_s, soma)) if len(r1_b) else 0
    out['global_carry'] = bool(n_r1 >= MIN_WIRINGS)
    n_rung = (n_ok(half(rg_b, gain_b) & half(rg_s, soma) & (rg_b > rp_b.reindex(rg_b.index))
                   & (rg_s > rp_s.reindex(rg_s.index))) if len(rg_b) else 0)
    out['terms_carry'] = bool(not out['global_carry'] and n_rung >= MIN_WIRINGS)
    out['counts'] = dict(soma_upstream=n_soma, boutons_halved=n_halved, soma_within=n_within, r1_carry=n_r1,
                         rung_carry=n_rung)
    out['thresholds'] = dict(soma_drop=thr_soma, soma_within=thr_within, bouton_halved=BOUTON_STRENGTH_HALVED)
    return out


def part_score() -> dict:
    """Step 2's runs scored and classified; refuses before E134's gate."""
    require_e134_gate()
    raise NotImplementedError("the score reads step 2's runs (atlas's stage readings on them) once they exist")


def print_commands() -> None:
    """The lanes' command lines (step 2): every lane runs the whole list and skips what another has claimed."""
    pre = 'OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    for host, (t, lanes) in LANES.items():
        for name, rev in lanes:
            print(f"{host}: {pre.format(t=t)} --part step2 --lane {name}{' --reverse' if rev else ''}")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('step0', 'step0_saved', 'lock', 'step1', 'step2', 'score', 'commands', 'gate'))
    ap.add_argument('--n-folds', type=int, default=None, help='lock: the number of odor folds (a planner\'s choice)')
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--lane', default=None)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    a = ap.parse_args(argv)
    if a.part == 'step0':
        part_step0(n_boot=a.n_boot, seed=a.seed)
    elif a.part == 'step0_saved':
        part_step0_saved()
    elif a.part == 'lock':
        write_lock(a.n_folds, a.seed)
    elif a.part == 'step1':
        part_step1()
    elif a.part == 'step2':
        part_step2(a.lane, a.reverse, a.reclaim)
    elif a.part == 'score':
        part_score()
    elif a.part == 'commands':
        print_commands()
    elif a.part == 'gate':
        print(json.dumps(e134_gate_status(), indent=1))


if __name__ == '__main__':
    main()
