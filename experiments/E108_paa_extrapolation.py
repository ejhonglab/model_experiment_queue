#!/usr/bin/env python3
"""E108: paa's input at the flies' concentration. paa's ORN input is moved from the recorded 10^-4 to the flies' 10^-3
by rules fixed before any run, and each input is read in the best family at fixed scales with both response rates
re-held.

The entry is item E108 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08 and 16:19; moved into the list 16:34 under the idle-host rule;
claimed by hal 16:37). Tom released olfdb's concentration ramps (dataset 21) at 16:13 ("all that data may be used").
"The family" is C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings
at its own fit (E83's `saved_fit`; S: real, P1, P2, M1; T: P3, P4, M2, M3), read as the mean matrix over the eight
wirings (E72's population). "The ceiling family" is E69's firing ceiling on four wirings (E69's fit on the real
wiring, E72 part G's on P1, P2 and M1; E82's `model_parts`), read beside. "The two ρ's" are ρ on odor-odor distances
with APL on against the control flies and with APL off against APL output blocked (TNT), on raw counts (cosine
distance on spike counts), each KC matrix the mean over recordings (E77's `data_distances('z8')`, the same matrices as
E106's `data_matrices`). ρ over paa's 15 pairs comes from E106's public `rho_pairs` (an explicit list of pairs).

Concentrations. Every output of this experiment carries this: paa's input at 10^-3 is an extrapolation, not a
recording (the lab's ORN recordings have paa at 10^-5 and 10^-4 only; the flies' KCs got 10^-3), and euc and -aPine
stay at the recorded 10^-1.75 and 10^-1.5 against the flies' 10^-2 (B-myr is at the flies' 10^-2).

The inputs (fixed before any run; glomeruli x the 16 odors, the lab's 36 rows; only paa's column differs from the
lab's input, `H.apl16_orn_deltas()`; inputs/<name>.csv with sha256):
- lo, hi: paa at 10^-5 (the validation2 ORN file's diagnostic column, `validation2_orn_deltas(drop_diags=False)`) and
  at 10^-4 (the panel's, the lab input's paa column), from the same six flies, in the model's units (estimated spike
  deltas). C (CAP): the largest value in the lab's input.
- X1, linear: hi + (hi - lo), capped at C (`rule_x1`).
- X2, saturating: hi + f (C - hi) with f = (hi - lo) / (C - lo), each glomerulus moving toward the common ceiling C (X1's
  cap) by the fraction of the remaining distance it covered from 10^-5 to 10^-4 (`rule_x2`).
- X3, by the lab's concentration ramps: top = a + b mid + c low, the response a decade up as a function of the two
  below, fitted by least squares on the ramps' triples (`fit_x3`, `rule_x3`).
- I0, the lab's input (the reference: the saved fits).
- Controls: W5, paa at 10^-5 in place of 10^-4 (the wrong way); P1 to P4, X1's changes (X1 - hi) permuted across the 36
  glomeruli (a derangement drawn with seed 1 to 4) and capped at C as X1 is; S1, paa scaled by X1's total over hi's
  (sum over the 36 glomeruli), its pattern unchanged.
- Step 4 (beside, no gate): B0 to B5, each validation2 ORN fly's own B-myr vector (olfdb dataset 4, the six flies), as
  E85's 'measured' model flies: the lab's B-myr plus the fly's deviation from the mean over the flies that have the
  glomerulus, in the model's units (each fly scaled to-avg-max, then one slope through the origin of the lab's four
  validation odors on the flies' mean; its R² is reported); paa's column is the lab's.

The ramps (part inputs; olfdb dataset 21, http://hal:8080/api/download?ids=21&format=long; the compact copy of what is
read goes to <shared>/data/E108_paa_extrapolation/ with a README):
- The seven single-odor ramps of panel mono_nat_ramp_cond (RAMP_ODORS: 2h, p-cre, IaA, t2h, menth, 2-but listed as
  2-but_pfo, 1p3ol; five concentrations a decade apart, 6 to 8 flies), on the model's 36 glomeruli where measured.
- The unit: pyolfsysm's dF/F-to-spike-delta fit (the one that made the model's input files; `dff_fit`): each fly scaled
  'to-avg-max' (its largest |mean over repeats| over all its panels, odors and ROIs brought to the saved avg_flymax),
  then the fit's slope through the origin (`load_data_and_refit_dff2spiking_model` on pyolfsysm's data/internal). The
  value of an entry is the mean over the flies that have it. `unit_check.csv`: the diagnostic odors that the ramp
  flies and the validation2 file share, the ramp flies' values against the file's (a check of the unit, which chooses
  nothing).
- Three triples of adjacent decades per odor and glomerulus (low, mid, top).
- The score that chooses the rule (`rule_scores`): the RMSE of the predicted top concentration over all triples (odors x
  the model's glomeruli, pooled). X1 and X2 fit nothing and are scored on every triple with C as above; X3 is scored
  with one odor's ramp left out at a time (fitted on the other six odors' triples). The rule with the lowest RMSE
  carries the verdict (`carrying_rule`); all three are run and reported. The six natural extracts (EXTRACTS) are scored
  beside, X3 fitted on all seven single odors, and choose nothing.

Parts (CLI `--part`):
- inputs (no simulations): the ramps' tables, scores.csv, scores_by_odor.csv, scores_extracts.csv, unit_check.csv,
  rule.json (the carrying rule and X3's coefficients), paa_inputs.csv (lo, hi and every rule's and control's paa column),
  bmyr_flies.csv, inputs/ with sha256.json, and the compact data copy.
- 0, the identity: on each of the eight wirings, the simulator built with the lab's input and its breadth recomputed
  (E85's `build`, with this experiment's cache) reproduces the saved APL-on and APL-off runs with 0 count and response
  differences, and its breadth equals the wiring's; on each of the ceiling family's four, E82's simulator with the
  breadth recomputed for the lab's input (`E77_input_screen.handle_breadth`) the same (identity.json). If anything
  differs the experiment stops there.
- B (simulations; resumable; `--lane a|b` takes every other point of the plan, so two processes run side by side):
  each planned point (`b_plan`: the rules, then the controls, on the eight wirings; then the same inputs on the
  ceiling family's four; then step 4's six flies on the eight wirings) at fixed scales with both response rates
  re-held (the breadth recomputed for the input, the somatic threshold bisected to the fit's APL-off response count,
  then the tied APL scale to its APL-on count; E85's `held_reading` for the family, `E77_input_screen.hold_rates` on
  E82's simulator for the ceiling family). Rows in runs_B_<lane>.csv (absolute parameters, the factors, realism, paa's
  and B-myr's response fractions, `H.code_commits()` of the process, which must be clean); runs under
  outputs/E108_paa_extrapolation/B/.
- score (no simulations): readings.csv (per input, unit and readout: the two ρ's against all recordings and halves A
  to D over the 120 pairs and over paa's 15; the APL-effect loss, its pattern part, amplitude part, shape part, spread
  ratio, ρ and Pearson r on the per-pair APL effect and the mean APL effect against all recordings and halves A to D;
  the block, the natural-source alignment on and off and the structure vs wild-type somata; the mean over the unit's
  runs of the response rates, response rate ρ on and off, and paa's and B-myr's response fractions; the runs outside
  the bands), runs.csv, gate.csv and gate.json (`gate`, per rule and control), classes.csv (`class_of`, with the column
  "provisional, E106 pending"), bmyr.csv (step 4), comparison_table.md.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.

The gate (the entry's, fixed before any run; `gate`; means over the eight wirings at fixed scales, each wiring paired
with itself at the lab's input; the two ρ's on raw counts, the APL-effect statistics under the calibrated readout
unless named):
- Gain: one ρ higher by GAIN (0.010) or more for the population against all recordings, on at least GAIN_WIRINGS (6) of
  the 8 wirings and for the population against each of halves A to D; the other ρ not lower by more than OTHER_SLACK
  (0.005) for the population against all recordings.
- The APL effect (E94's equal): the pattern part not higher by more than EQUAL_PATTERN (0.003) against all recordings
  and against each of halves A to D; the raw-count loss not higher by more than EQUAL_RAW (0.005); the mean APL effect
  not lower by more than EQUAL_MEAN (0.01). A lower loss does not fail.
- Realism: every one of the input's eight runs inside the four bands.
- Response rates: response rate ρ with APL on and with APL off (the bands' units) not lower by more than RATE_RHO_SLACK
  (0.05) as the mean over the eight wirings. paa's own response fractions are reported against the flies' (Yang's z,
  0.065 / 0.137, as E64 read them; and in baseline-noise units beside).
- Pending Tom's answers of 2026-11-17 14:44 (applied until he decides; the class is also given without them): the block
  not lower by more than BLOCK_SLACK (0.010); the natural-source alignment over the 66 megamat pairs with APL on and with
  APL off, and the structure vs wild-type somata (APL on), each not lower by more than STRUCTURE_SLACK (0.05); all on the
  population's raw-count matrices.
- Controls: W5's and the mean over P1 to P4's gains, on each ρ that gains, less than CONTROL_RATIO (0.5) times the
  candidate's.
- Classes (`class_of`): 'better' (all hold), 'a price' (the gain holds and a constraint fails; each failure named with
  its size; a failure of the controls alone is named the same way), 'no gain'. The outcome letter (a) to (d) of the
  entry beside: (c) / (d) by whether W5 lowers either ρ by GAIN or more.
- The class is written after E106's note (its yardsticks); until then classes.csv says "provisional, E106 pending".
  No threshold here moves after a run is read.

Choices the entry left open (fixed here before any run, with the reason):
- X3's form: ordinary least squares top = a + b mid + c low, pooled over the training odors' triples and glomeruli, no
  cap (the simplest function of the two below; it holds X1 as a special case, a=0, b=2, c=-1).
- X2 where the ceiling is not above the low value (C <= low, which no paa glomerulus meets): f = 0, the prediction is
  the mid value (no remaining distance).
- X1 and X2 have no floor: the entry caps X1 at the largest rate and names no floor, so a falling glomerulus keeps
  falling (DC1 and VA3 fall from 10^-5 to 10^-4). paa_inputs.csv gives the smallest values beside the lab input's.
- The permuted controls are re-capped at C after the permutation (X1's rule caps every value); S1 is not capped
  (capping would change its pattern).
- The controls of the gate are X1's, as the entry defines them, whichever rule carries the verdict.
- The ceiling family is read for every input and control; step 4 on the eight wirings of the family only.

Departures from the entry: none so far.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no connectome table is read beyond the model's own
hemibrain wiring and E72's re-draws of it (the eight wirings, the four of the ceiling family); no raw-only recording.
Nothing is published to claude.ai.

Outputs: results/E108_paa_extrapolation/; spike counts and responses under outputs/E108_paa_extrapolation/B/;
<shared>/data/E108_paa_extrapolation/ (the compact data).

Usage (from the repository root; two lanes at most on hal, under systemd-run --user):
    ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part B --lane a
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part B --lane b
    ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part score
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import io
import json
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E108_paa_extrapolation'
OLFDB = 'http://hal:8080'
OLFDB_RAMPS, OLFDB_VALIDATION = 21, 4
RAMP_PANEL = 'mono_nat_ramp_cond'
RAMP_ODORS = ('2h', 'p-cre', 'IaA', 't2h', 'menth', '2-but_pfo', '1p3ol')
EXTRACTS = ('ACV', 'ban', 'basil', 'mint', 'pearb', 'wine')
DIAG_PANEL = 'diagnostics'                      # olfdb's name; the validation2 file's is 'glomeruli_diagnostics'
PAA_LO_COL = ('glomeruli_diagnostics', 'paa @ -5')
PAA_HI_COL = ('validation2', 'paa @ -4')
BMYR = 'B-myr'
RULES = ('X1', 'X2', 'X3')
PERM_SEEDS = (1, 2, 3, 4)
PERMS = tuple(f'P{s}' for s in PERM_SEEDS)
CONTROLS = ('W5',) + PERMS + ('S1',)
INPUTS = ('I0',) + RULES + CONTROLS
N_BMYR_FLIES = 6
BMYR_INPUTS = tuple(f'B{f}' for f in range(N_BMYR_FLIES))
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
FC_MODELS = {'fc_real': 'fit1', 'fc_P1': 'P1_fc', 'fc_P2': 'P2_fc', 'fc_M1': 'M1_fc'}
FC_WIRINGS = tuple(FC_MODELS)
UNITS = (('all8', ALL_WIRINGS), ('S', S_WIRINGS), ('T', T_WIRINGS), ('fc4', FC_WIRINGS))
READOUTS = ('raw', 'cal')
HALVES = ('A', 'B', 'C', 'D')
TARGETS = ('all',) + HALVES                    # 'all' = all recordings (E83's 'z8')
PAA15 = tuple(('paa', o) for o in H.APL16_ODORS if o != 'paa')
FLIES_PAA_FRAC = {'on': 0.065, 'off': 0.137}    # the entry's, Yang's z (as E64 read them; `A.frac_data(c, 'z')['paa']`)

# the gate (fixed in the entry)
GAIN = 0.010
GAIN_WIRINGS = 6
OTHER_SLACK = 0.005
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
RATE_RHO_SLACK = 0.05
BLOCK_SLACK = 0.010
STRUCTURE_SLACK = 0.05
CONTROL_RATIO = 0.5
E106_PENDING = 'provisional, E106 pending'

_CACHE: Dict[object, object] = {}


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def data_dir() -> Path:
    """<shared>/data/E108_paa_extrapolation/: the compact copy of what this experiment reads from olfdb."""
    d = H.shared_dir() / 'data' / EXPERIMENT
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def _flagged(x) -> bool:
    return isinstance(x, str) and x.strip() not in ('', 'nan')


def _commits_clean() -> dict:
    """`H.code_commits()` at the start of a simulation part, which must be clean (E89's)."""
    from experiments import E89_apl_breadth_weighting as E89
    return E89._commits_clean()


# ---------------------------------------------------------------------------
# the rules (made-up data in the tests)

def rule_x1(lo: np.ndarray, hi: np.ndarray, cap: float) -> np.ndarray:
    """X1, linear: `hi` plus its change from `lo`, capped at `cap` (no floor)."""
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    return np.minimum(hi + (hi - lo), float(cap))


def rule_x2(lo: np.ndarray, hi: np.ndarray, ceiling: float) -> np.ndarray:
    """X2, saturating: `hi` moved toward `ceiling` by the fraction of the remaining distance that it covered from `lo`
    to `hi`, hi + f (ceiling - hi) with f = (hi - lo) / (ceiling - lo); f = 0 (the prediction `hi`) where the ceiling is
    not above `lo`. A falling glomerulus (f < 0) moves away from the ceiling by the same rule."""
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    room = float(ceiling) - lo
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(room > 0, (hi - lo) / np.where(room > 0, room, 1.0), 0.0)
    return hi + f * (float(ceiling) - hi)


def fit_x3(triples: pd.DataFrame) -> np.ndarray:
    """X3's coefficients (a, b, c) of top = a + b mid + c low: ordinary least squares over the rows of `triples`
    (columns 'low', 'mid', 'top'; rows with any NaN dropped)."""
    t = triples[['low', 'mid', 'top']].astype(float).dropna()
    if len(t) < 3:
        raise ValueError(f'X3 needs at least 3 triples, got {len(t)}')
    X = np.column_stack([np.ones(len(t)), t['mid'].to_numpy(), t['low'].to_numpy()])
    coef, *_ = np.linalg.lstsq(X, t['top'].to_numpy(), rcond=None)
    return coef


def rule_x3(lo: np.ndarray, hi: np.ndarray, coef: np.ndarray) -> np.ndarray:
    """X3: a + b hi + c lo (`fit_x3`'s coefficients; `lo` the low value, `hi` the mid one)."""
    a, b, c = (float(x) for x in coef)
    return a + b * np.asarray(hi, dtype=float) + c * np.asarray(lo, dtype=float)


def predictions(triples: pd.DataFrame, cap: float) -> pd.DataFrame:
    """Every triple's predicted top by X1 and X2 (no fit) and by X3 left one odor out ('X3', fitted on the other odors'
    triples): `triples` with columns 'pred_X1', 'pred_X2', 'pred_X3' added."""
    t = triples.copy()
    t['pred_X1'] = rule_x1(t['low'], t['mid'], cap)
    t['pred_X2'] = rule_x2(t['low'], t['mid'], cap)
    t['pred_X3'] = np.nan
    for o in t['odor'].unique():
        held = t['odor'] == o
        coef = fit_x3(t[~held])
        t.loc[held, 'pred_X3'] = rule_x3(t.loc[held, 'low'], t.loc[held, 'mid'], coef)
    return t


def rule_scores(triples: pd.DataFrame, cap: float) -> pd.DataFrame:
    """The score that chooses the rule: per rule (index RULES), the RMSE of the predicted top concentration over all the
    triples (odors x glomeruli, pooled; X1 and X2 on every triple, X3 with one odor's ramp left out at a time), the
    number of triples, the bias (mean predicted minus measured) and the RMSE of a prediction of no change (top = mid)
    beside."""
    p = predictions(triples, cap)
    rows = {}
    for r in RULES:
        e = (p[f'pred_{r}'] - p['top']).to_numpy(dtype=float)
        rows[r] = dict(rmse=_rms(e), bias=float(np.mean(e)), n_triples=int(len(e)), n_odors=int(p['odor'].nunique()),
                       n_capped=int((p['mid'] + (p['mid'] - p['low']) > cap).sum()) if r == 'X1' else 0)
    rows['none (top = mid)'] = dict(rmse=_rms((p['mid'] - p['top']).to_numpy(dtype=float)),
                                    bias=float(np.mean(p['mid'] - p['top'])), n_triples=int(len(p)),
                                    n_odors=int(p['odor'].nunique()), n_capped=0)
    out = pd.DataFrame.from_dict(rows, orient='index')
    out.index.name = 'rule'
    return out


def scores_by_odor(triples: pd.DataFrame, cap: float) -> pd.DataFrame:
    """`rule_scores`' RMSE per held-out odor (X3 left that odor out)."""
    p = predictions(triples, cap)
    rows = []
    for o, g in p.groupby('odor', sort=False):
        rows.append(dict(odor=o, n_triples=len(g), **{r: _rms((g[f'pred_{r}'] - g['top']).to_numpy(dtype=float)) for r in RULES},
                         none=_rms((g['mid'] - g['top']).to_numpy(dtype=float))))
    return pd.DataFrame(rows)


def carrying_rule(scores: pd.DataFrame) -> str:
    """The rule of RULES with the lowest RMSE (`rule_scores`'s table)."""
    s = scores.loc[[r for r in RULES if r in scores.index], 'rmse'].astype(float)
    if s.isna().any():
        raise ValueError(f'a rule has no score: {s.to_dict()}')
    return str(s.idxmin())


def derangement(n: int, seed: int) -> np.ndarray:
    """The first permutation of range(n) drawn with `default_rng(seed)` that leaves no element in place."""
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


def paa_columns(lo: pd.Series, hi: pd.Series, cap: float, coef: np.ndarray) -> pd.DataFrame:
    """paa's column (index the glomeruli) of every rule and control: 'lo', 'hi', X1, X2, X3, W5, P1 to P4 (X1's
    changes, X1 - hi, permuted across the glomeruli by `derangement(n, seed)` and capped at `cap`), S1 (hi times X1's
    total over hi's)."""
    if not lo.index.equals(hi.index):
        raise ValueError('lo and hi have different glomeruli')
    x1 = rule_x1(lo, hi, cap)
    out = pd.DataFrame(dict(lo=lo.to_numpy(float), hi=hi.to_numpy(float), X1=x1, X2=rule_x2(lo, hi, cap),
                            X3=rule_x3(lo, hi, coef), W5=lo.to_numpy(float)), index=hi.index)
    d = x1 - hi.to_numpy(float)
    for s, name in zip(PERM_SEEDS, PERMS):
        out[name] = np.minimum(hi.to_numpy(float) + d[derangement(len(d), s)], cap)
    out['S1'] = hi.to_numpy(float) * (float(x1.sum()) / float(hi.sum()))
    return out


# ---------------------------------------------------------------------------
# the data

def olfdb_long(ds: int, refresh: bool = False) -> Tuple[pd.DataFrame, dict]:
    """olfdb dataset `ds` in its long format and its metadata, cached under outputs/E108_paa_extrapolation/olfdb/
    (fetched from OLFDB the first time)."""
    from urllib.request import urlopen
    d = local_dir('olfdb')
    p, pm = d / f'dataset{ds}.parquet', d / f'dataset{ds}.json'
    if refresh or not (p.exists() and pm.exists()):
        with urlopen(f'{OLFDB}/api/download?ids={ds}&format=long') as r:
            raw = r.read()
        df = pd.read_csv(io.BytesIO(raw))
        with urlopen(f'{OLFDB}/api/datasets/{ds}') as r:
            meta = json.loads(r.read())
        df.to_parquet(p)
        keep = {k: meta.get(k) for k in ('id', 'title', 'uploader_name', 'created_at', 'panel_names', 'indicator', 'driver')}
        keep.update(download_sha256=hashlib.sha256(raw).hexdigest(), n_rows=int(len(df)),
                    url=f'{OLFDB}/api/download?ids={ds}&format=long')
        _write_json(pm, keep)
    return pd.read_parquet(p), _read_json(pm)


def dff_fit() -> Dict[str, float]:
    """pyolfsysm's dF/F -> spike-delta fit, the one that made the model's input files: 'slope' (refitted from its saved
    data, `load_data_and_refit_dff2spiking_model`), 'avg_flymax' (its saved choice) and 'scaling_method' (must be
    'to-avg-max'; no constant, no separate inhibitory fit)."""
    if 'dff' not in _CACHE:
        import contextlib
        from pyolfsysm import mb_model as M
        md = Path(M.data_root) / 'internal'
        ch = M.read_dff_to_spiking_model_choices(md)
        if str(ch['scaling_method']) != 'to-avg-max' or bool(ch['add_constant']) or bool(ch['separate_inh_model']):
            raise RuntimeError(f'unexpected dF/F -> spike choices {ch.to_dict()}')
        with contextlib.redirect_stdout(io.StringIO()):
            model, inh = M.load_data_and_refit_dff2spiking_model(md, verbose=False)
        if inh is not None or len(model.params) != 1:
            raise RuntimeError(f'unexpected dF/F -> spike model {model.params.to_dict()}')
        _CACHE['dff'] = dict(slope=float(model.params.iloc[0]), avg_flymax=float(ch['avg_flymax']), scaling_method='to-avg-max')
    return dict(_CACHE['dff'])


def fly_means(df: pd.DataFrame) -> pd.DataFrame:
    """Per fly, panel, odor and ROI the mean over repeats (columns date, fly_num, panel, odor1, roi, value)."""
    return df.groupby(['date', 'fly_num', 'panel', 'odor1', 'roi'], as_index=False)['value'].mean()


def fly_maxima(fm: pd.DataFrame) -> pd.DataFrame:
    """Each fly's largest |mean over repeats| over all its panels, odors and ROIs ('to-avg-max''s divisor)."""
    m = fm.assign(a=fm['value'].abs()).groupby(['date', 'fly_num'], as_index=False)['a'].max()
    return m.rename(columns={'a': 'fly_max'})


def to_spike_deltas(fm: pd.DataFrame, maxima: pd.DataFrame, fit: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """`fm` (`fly_means`' layout) with 'rate': value / fly_max x avg_flymax x slope (pyolfsysm's fit)."""
    fit = dff_fit() if fit is None else fit
    m = fm.merge(maxima, on=['date', 'fly_num'], how='left', validate='many_to_one')
    if m['fly_max'].isna().any():
        raise RuntimeError('a fly has no maximum')
    return m.assign(rate=m['value'] / m['fly_max'] * fit['avg_flymax'] * fit['slope'])


def compact_tables(refresh: bool = False) -> Dict[str, pd.DataFrame]:
    """What this experiment reads from olfdb, compact (also written to `data_dir()` by `write_compact`): 'ramps' (dataset
    21's ramp panel and diagnostics, mean over repeats per fly, odor and ROI), 'ramp_maxima' (its flies' maxima),
    'validation' (dataset 4's validation2 panel and diagnostics, the same layout), 'validation_maxima'; from the shared
    copy when olfdb has not been read on this host (`refresh` reads olfdb again)."""
    if 'compact' in _CACHE and not refresh:
        return _CACHE['compact']
    names = ('ramps', 'ramp_maxima', 'validation', 'validation_maxima')
    local = {n: local_dir('olfdb') / f'{n}.csv' for n in names}
    shared = {n: H.shared_dir() / 'data' / EXPERIMENT / f'{n}.csv' for n in names}
    if not refresh and all(p.exists() for p in local.values()):
        out = {n: pd.read_csv(local[n]) for n in names}
    elif not refresh and all(p.exists() for p in shared.values()):
        out = {n: pd.read_csv(shared[n]) for n in names}
    else:
        r, _ = olfdb_long(OLFDB_RAMPS, refresh)
        v, _ = olfdb_long(OLFDB_VALIDATION, refresh)
        fr, fv = fly_means(r), fly_means(v)
        out = dict(ramps=fr, ramp_maxima=fly_maxima(fr), validation=fv, validation_maxima=fly_maxima(fv))
        for n in names:
            out[n].to_csv(local[n], index=False)
    for n in ('ramps', 'validation'):
        out[n]['date'] = out[n]['date'].astype(str)
        out[n]['fly_num'] = out[n]['fly_num'].astype(int)
    for n in ('ramp_maxima', 'validation_maxima'):
        out[n]['date'] = out[n]['date'].astype(str)
        out[n]['fly_num'] = out[n]['fly_num'].astype(int)
    _CACHE['compact'] = out
    return out


def write_compact() -> dict:
    """The compact copy and its README in `data_dir()` (`compact_tables`), with olfdb's metadata of both datasets."""
    t = compact_tables()
    d = data_dir()
    rec = {}
    for n, df in t.items():
        p = d / f'{n}.csv'
        p.write_text(df.to_csv(index=False))
        rec[n] = dict(file=p.name, rows=int(len(df)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    metas = {ds: _read_json(local_dir('olfdb') / f'dataset{ds}.json') for ds in (OLFDB_RAMPS, OLFDB_VALIDATION)}
    _write_json(d / 'olfdb_metadata.json', {str(k): v for k, v in metas.items()})
    fit = dff_fit()
    lines = [
        '# E108: what paa\'s extrapolation reads from olfdb (compact copy)', '',
        'Written by experiments/E108_paa_extrapolation.py (`write_compact`, part inputs) in the model_experiment_queue '
        'repository, so that the inputs can be rebuilt without olfdb. Each file is the mean over the three repeats per fly, '
        'panel, odor and ROI (dF/F as olfdb serves it), or each fly\'s largest |mean| over all its panels, odors and ROIs '
        '(the divisor of pyolfsysm\'s to-avg-max scaling).', '',
        f'- ramps.csv: olfdb dataset {OLFDB_RAMPS} ({metas[OLFDB_RAMPS].get("title")}; uploaded by '
        f'{metas[OLFDB_RAMPS].get("uploader_name")}), panels {RAMP_PANEL} and diagnostics, all ROIs; download sha256 '
        f'{metas[OLFDB_RAMPS].get("download_sha256")}. Tom released it for queue work 2026-11-17 16:13.',
        '- ramp_maxima.csv: its flies\' maxima.',
        f'- validation.csv: olfdb dataset {OLFDB_VALIDATION} ({metas[OLFDB_VALIDATION].get("title")}), panels validation2 '
        f'and diagnostics; download sha256 {metas[OLFDB_VALIDATION].get("download_sha256")}.',
        '- validation_maxima.csv: its flies\' maxima.',
        '- olfdb_metadata.json: both datasets\' metadata as olfdb gave them.', '',
        f'The unit (estimated spike deltas): value / fly_max x {fit["avg_flymax"]:.6f} (avg_flymax) x {fit["slope"]:.6f} '
        '(the slope of pyolfsysm\'s dF/F -> spike-delta fit, refitted from pyolfsysm/data/internal).', '',
        '| file | rows | sha256 |', '|---|---|---|'] + [f'| {r["file"]} | {r["rows"]} | {r["sha256"][:16]} |' for r in rec.values()]
    (d / 'README.md').write_text('\n'.join(lines) + '\n')
    return rec


def ramp_rates() -> pd.DataFrame:
    """The ramp panel in estimated spike deltas, the mean over the flies that have each entry: (odor, conc, glomerulus,
    rate, n_flies), for the single-odor ramps and the extracts, on every ROI."""
    t = compact_tables()
    fm = t['ramps'][t['ramps']['panel'] == RAMP_PANEL]
    s = to_spike_deltas(fm, t['ramp_maxima'])
    s = s.assign(odor=[o.split(' @ ')[0] for o in s['odor1']], conc=[float(o.split(' @ ')[1]) for o in s['odor1']])
    g = s.groupby(['odor', 'conc', 'roi'])
    out = g['rate'].mean().rename('rate').to_frame().join(g['rate'].size().rename('n_flies')).reset_index()
    return out.rename(columns={'roi': 'glomerulus'})


def ramp_triples(rates: Optional[pd.DataFrame] = None, odors: Sequence[str] = RAMP_ODORS,
                 glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The triples of adjacent decades per ramp and glomerulus of the model's 36 (default; `glomeruli`): columns odor,
    glomerulus, conc_low, low, mid, top (three per odor with five concentrations); only where all three are measured."""
    rates = ramp_rates() if rates is None else rates
    gl = list(lab_input().index) if glomeruli is None else list(glomeruli)
    rows = []
    for o in odors:
        r = rates[(rates['odor'] == o) & rates['glomerulus'].isin(gl)]
        cs = sorted(r['conc'].unique())
        if len(cs) < 3:
            raise RuntimeError(f'{o}: {len(cs)} concentrations')
        if not np.allclose(np.diff(cs), 1.0):
            raise RuntimeError(f'{o}: concentrations {cs} are not a decade apart')
        wide = r.pivot_table(index='glomerulus', columns='conc', values='rate')
        for i in range(len(cs) - 2):
            c = cs[i:i + 3]
            for g in wide.index:
                v = wide.loc[g, c].to_numpy(dtype=float)
                if np.isfinite(v).all():
                    rows.append(dict(odor=o, glomerulus=g, conc_low=c[0], low=v[0], mid=v[1], top=v[2]))
    return pd.DataFrame(rows)


def lab_input() -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()`: 36 glomeruli x (panel, odor))."""
    if 'lab' not in _CACHE:
        _CACHE['lab'] = H.apl16_orn_deltas().astype(float)
    return _CACHE['lab'].copy()


def cap() -> float:
    """C: the largest value in the lab's input (X1's cap and X2's ceiling)."""
    return float(lab_input().to_numpy(dtype=float).max())


def paa_rates() -> pd.DataFrame:
    """paa at 10^-5 ('lo', the validation2 file's diagnostic column) and 10^-4 ('hi', the lab input's paa column, the
    file's panel column), the lab input's 36 glomeruli."""
    from pyolfsysm.mb_model import validation2_orn_deltas
    lab = lab_input()
    v = validation2_orn_deltas(drop_diags=False)
    hi_col = [c for c in lab.columns if H.odor_abbrev(c) == 'paa']
    if hi_col != [PAA_HI_COL]:
        raise RuntimeError(f"the lab input's paa column is {hi_col}, not {PAA_HI_COL}")
    if not np.array_equal(v.loc[lab.index, PAA_HI_COL].to_numpy(float), lab[PAA_HI_COL].to_numpy(float)):
        raise RuntimeError("the validation2 file's paa @ -4 is not the lab input's")
    out = pd.DataFrame(dict(lo=v.loc[lab.index, PAA_LO_COL].to_numpy(float), hi=lab[PAA_HI_COL].to_numpy(float)), index=lab.index)
    if out.isna().any().any():
        raise RuntimeError('paa is missing in some glomerulus')
    return out


def triples() -> pd.DataFrame:
    if 'triples' not in _CACHE:
        _CACHE['triples'] = ramp_triples()
    return _CACHE['triples'].copy()


def x3_coef() -> np.ndarray:
    """X3's coefficients fitted on all seven single-odor ramps' triples (applied to paa)."""
    return fit_x3(triples())


def paa_table() -> pd.DataFrame:
    """`paa_columns` for paa (lo, hi, every rule and control)."""
    if 'paa_table' not in _CACHE:
        p = paa_rates()
        _CACHE['paa_table'] = paa_columns(p['lo'], p['hi'], cap(), x3_coef())
    return _CACHE['paa_table'].copy()


def with_paa(column: np.ndarray) -> pd.DataFrame:
    """The lab's input with paa's column replaced."""
    orn = lab_input()
    orn[PAA_HI_COL] = np.asarray(column, dtype=float)
    return orn


# step 4: B-myr per validation2 fly

def bmyr_flies() -> dict:
    """The six validation2 ORN flies' B-myr in the model's units (E85's 'measured' model flies): 'keys' ((date,
    fly_num)), 'dev' (flies x the 36 glomeruli: the fly's value minus the mean over the flies that have the glomerulus;
    0 where the fly lacks it), 'slope' (least squares through the origin of the lab's four validation odors on the
    flies' mean, after to-avg-max), 'r2' / 'r2_origin', 'n_measured'."""
    if 'bmyr' in _CACHE:
        return _CACHE['bmyr']
    t = compact_tables()
    lab = H.abbrev_columns(lab_input())
    gl = list(lab.index)
    fm = t['validation'][t['validation']['panel'] == 'validation2']
    m = fm.merge(t['validation_maxima'], on=['date', 'fly_num'], validate='many_to_one')
    A = dff_fit()['avg_flymax']
    m = m.assign(scaled=m['value'] / m['fly_max'] * A, odor=[o.split(' @ ')[0] for o in m['odor1']])
    want = {o: f'{o} @ {H.ORN_INPUT_CONC[o]:g}' for o in H.VALIDATION2_4}
    m = m[m['odor1'].isin(set(want.values())) & m['roi'].isin(gl)]
    keys = sorted({(str(d), int(f)) for d, f in zip(m['date'], m['fly_num'])})
    if len(keys) != N_BMYR_FLIES:
        raise RuntimeError(f'expected {N_BMYR_FLIES} validation2 flies, found {len(keys)}')
    mean = m.groupby(['odor', 'roi'])['scaled'].mean()
    x = np.array([mean.get((o, g), np.nan) for o in H.VALIDATION2_4 for g in gl])
    y = np.array([lab.loc[g, o] for o in H.VALIDATION2_4 for g in gl], dtype=float)
    ok = np.isfinite(x)
    slope = float(x[ok] @ y[ok] / (x[ok] @ x[ok]))
    res = y[ok] - slope * x[ok]
    r2 = float(1 - res @ res / ((y[ok] - y[ok].mean()) @ (y[ok] - y[ok].mean())))
    r2o = float(1 - res @ res / (y[ok] @ y[ok]))
    b = m[m['odor'] == BMYR]
    V = np.full((len(keys), len(gl)), np.nan)
    for (d, f, g), v in b.groupby(['date', 'fly_num', 'roi'])['scaled'].mean().items():
        V[keys.index((str(d), int(f))), gl.index(g)] = v * slope
    with np.errstate(invalid='ignore'):
        mu = np.nanmean(V, axis=0)
    dev = np.where(np.isfinite(V), V - mu[None], 0.0)
    _CACHE['bmyr'] = dict(keys=keys, dev=dev, slope=slope, r2=r2, r2_origin=r2o, n_measured=np.isfinite(V).sum(axis=1),
                          values=V, glomeruli=gl)
    return _CACHE['bmyr']


def bmyr_input(f: int) -> pd.DataFrame:
    """Model fly `f`'s input: the lab's with B-myr's column plus the fly's deviation (`bmyr_flies`)."""
    B = bmyr_flies()
    orn = lab_input()
    col = [c for c in orn.columns if H.odor_abbrev(c) == BMYR]
    if len(col) != 1 or list(orn.index) != B['glomeruli']:
        raise RuntimeError("the lab input's B-myr column or glomeruli are not the flies'")
    orn[col[0]] = orn[col[0]].to_numpy(float) + B['dev'][f]
    return orn


def input_of(name: str) -> pd.DataFrame:
    """An input by name: 'I0' (the lab's), a rule (RULES), a control (CONTROLS) or a B-myr fly (BMYR_INPUTS); the lab
    input's layout."""
    c = _CACHE.setdefault('inputs', {})
    if name not in c:
        if name == 'I0':
            c[name] = lab_input()
        elif name in RULES or name in CONTROLS:
            c[name] = with_paa(paa_table()[name].to_numpy(float))
        elif name in BMYR_INPUTS:
            c[name] = bmyr_input(int(name[1:]))
        else:
            raise ValueError(f'unknown input {name!r}')
    return c[name].copy()


def input_path(name: str) -> Path:
    return out_dir() / 'inputs' / f'{name}.csv'


def write_input(name: str) -> dict:
    """inputs/<name>.csv (glomeruli x abbreviated odors) and its record in inputs/sha256.json."""
    from experiments.E77_input_screen import frame_sha256
    orn = input_of(name)
    lay = H.abbrev_columns(orn)
    lay.index.name = 'glomerulus'
    p = input_path(name)
    p.parent.mkdir(exist_ok=True)
    p.write_text(lay.to_csv())
    rec = dict(file=str(p.relative_to(H.shared_dir())), file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
               frame_sha256=frame_sha256(orn))
    shas = _read_json(p.parent / 'sha256.json')
    shas[name] = rec
    _write_json(p.parent / 'sha256.json', shas)
    return rec


def unit_check() -> pd.DataFrame:
    """The diagnostic odors (at the same concentration) that the ramp flies (dataset 21) and the validation2 file share:
    per odor, the slope through the origin of the file's values on the ramp flies' (both in the model's units), Pearson r
    and the number of glomeruli; and 'pooled'. Chooses nothing."""
    from pyolfsysm.mb_model import validation2_orn_deltas
    t = compact_tables()
    v = validation2_orn_deltas(drop_diags=False)
    vd = {c[1]: v[c] for c in v.columns if c[0] == 'glomeruli_diagnostics'}
    fm = t['ramps'][t['ramps']['panel'] == DIAG_PANEL]
    s = to_spike_deltas(fm, t['ramp_maxima']).groupby(['odor1', 'roi'])['rate'].mean()
    rows, xs, ys = [], [], []
    for o in sorted(set(s.index.get_level_values(0)) & set(vd)):
        g = [x for x in s.loc[o].index if x in vd[o].index and np.isfinite(vd[o][x])]
        if len(g) < 3:
            continue
        x, y = s.loc[o].reindex(g).to_numpy(float), vd[o].reindex(g).to_numpy(float)
        xs.append(x)
        ys.append(y)
        rows.append(dict(odor=o, n_glomeruli=len(g), slope=float(x @ y / (x @ x)), r=float(np.corrcoef(x, y)[0, 1]),
                         ramp_flies_max=float(x.max()), file_max=float(y.max())))
    if xs:
        x, y = np.concatenate(xs), np.concatenate(ys)
        rows.append(dict(odor='pooled', n_glomeruli=len(x), slope=float(x @ y / (x @ x)), r=float(np.corrcoef(x, y)[0, 1]),
                         ramp_flies_max=float(x.max()), file_max=float(y.max())))
    return pd.DataFrame(rows)


def write_inputs() -> dict:
    """Part inputs (no simulations): the compact data copy, the ramps' triples and the scores (scores.csv,
    scores_by_odor.csv, scores_extracts.csv), unit_check.csv, rule.json, paa_inputs.csv, bmyr_flies.csv and inputs/."""
    out = out_dir()
    rec = dict(compact=write_compact(), dff_fit=dff_fit(), cap=cap())
    rates = ramp_rates()
    rates.to_csv(out / 'ramp_rates.csv', index=False)
    tr = triples()
    tr.to_csv(out / 'ramp_triples.csv', index=False)
    sc = rule_scores(tr, cap())
    sc.to_csv(out / 'scores.csv')
    scores_by_odor(tr, cap()).to_csv(out / 'scores_by_odor.csv', index=False)
    # the extracts, beside (X3 fitted on all seven single odors)
    ext = ramp_triples(rates, EXTRACTS)
    coef = x3_coef()
    rows = []
    for o, g in ext.groupby('odor', sort=False):
        rows.append(dict(odor=o, n_triples=len(g), X1=_rms(rule_x1(g['low'], g['mid'], cap()) - g['top']),
                         X2=_rms(rule_x2(g['low'], g['mid'], cap()) - g['top']),
                         X3=_rms(rule_x3(g['low'], g['mid'], coef) - g['top']), none=_rms(g['mid'] - g['top'])))
    allx = dict(odor='all extracts', n_triples=len(ext), X1=_rms(rule_x1(ext['low'], ext['mid'], cap()) - ext['top']),
                X2=_rms(rule_x2(ext['low'], ext['mid'], cap()) - ext['top']), X3=_rms(rule_x3(ext['low'], ext['mid'], coef) - ext['top']),
                none=_rms(ext['mid'] - ext['top']))
    pd.DataFrame(rows + [allx]).to_csv(out / 'scores_extracts.csv', index=False)
    unit_check().to_csv(out / 'unit_check.csv', index=False)
    rule = carrying_rule(sc)
    pt = paa_table()
    lab = lab_input()
    pt.to_csv(out / 'paa_inputs.csv')
    B = bmyr_flies()
    bm = pd.DataFrame(B['values'], index=[f'B{i}' for i in range(len(B['keys']))], columns=B['glomeruli'])
    bm.insert(0, 'fly', [f'{d} fly {f}' for d, f in B['keys']])
    bm.to_csv(out / 'bmyr_flies.csv')
    rec.update(carrying_rule=rule, x3_coef=dict(zip(('a', 'b_mid', 'c_low'), (float(x) for x in coef))),
               scores=sc.reset_index().to_dict('records'), n_triples=int(len(tr)), ramp_odors=list(RAMP_ODORS),
               paa=dict(lab_min=float(lab.to_numpy().min()), lab_max=float(lab.to_numpy().max()),
                        **{f'{k}_min': float(pt[k].min()) for k in pt.columns},
                        **{f'{k}_max': float(pt[k].max()) for k in pt.columns},
                        **{f'{k}_sum': float(pt[k].sum()) for k in pt.columns}),
               bmyr=dict(flies=[f'{d} fly {f}' for d, f in B['keys']], slope=B['slope'], r2=B['r2'], r2_origin=B['r2_origin'],
                         n_measured=[int(x) for x in B['n_measured']]))
    rec['inputs'] = {n: write_input(n) for n in INPUTS + BMYR_INPUTS}
    rec.update(H.code_commits())
    _write_json(out / 'rule.json', rec)
    print(f'E108 inputs: {len(tr)} triples; RMSE ' + ', '.join(f'{r} {sc.loc[r, "rmse"]:.2f}' for r in sc.index)
          + f'; carrying rule {rule}; X3 coef {np.round(coef, 4).tolist()}', flush=True)
    return rec


# ---------------------------------------------------------------------------
# the models: the family's eight wirings and the ceiling family's four

def is_fc(unit: str) -> bool:
    return unit in FC_MODELS


def fc_parts(unit: str) -> dict:
    """E82's `model_parts` of a ceiling-family wiring, with 'W' (its claw x bouton wiring; None for the hemibrain's)."""
    from experiments import E82_protocol_runs as R
    m = dict(R.model_parts(FC_MODELS[unit]))
    m['W'] = None if m['wiring'] is None else R.wiring_parts(m['wiring'])['W']
    return m


def fc_build(unit: str, orn: pd.DataFrame, breadth: Optional[pd.Series] = None):
    """The ceiling-family simulator of `unit` on the input `orn` (E82's `build_simulator`: its refractory period and
    adaptation, C2's normalized APL weights, breadth exponent 2 with the breadth recomputed for `orn` on its wiring,
    `E77_input_screen.handle_breadth`). Returns (simulate, info)."""
    from experiments import E77_input_screen as S
    from experiments import E82_protocol_runs as R
    m = fc_parts(unit)
    b = S.handle_breadth(None, orn=orn, wPNKC=m['W']) if breadth is None else breadth
    return R.build_simulator(dict(m['kc']), wPNKC=m['W'], breadth=b, params=m['p'], orn=orn)


def fit_of(unit: str) -> dict:
    """A unit's fit: 'saved', 'p', 'target_on', 'target_off', 'folder', 'tag'."""
    if is_fc(unit):
        m = fc_parts(unit)
        return dict(saved=m['saved'], p=dict(m['p']), target_on=m['target_on'], target_off=m['target_off'],
                    folder=m['source'][0], tag=m['source'][1])
    P = E83.saved_fit(unit)
    return dict(saved=P['saved'], p=dict(P['p']), target_on=P['target_on'], target_off=P['target_off'], folder=P['folder'],
                tag=P['tag'])


def cache_tag(orn: pd.DataFrame, name: str) -> str:
    return E83.input_tag_for(orn, f'E108_{name}')


def identity(units: Sequence[str] = ALL_WIRINGS + FC_WIRINGS) -> dict:
    """Part 0 (the module docstring); resumable per unit; identity.json; raises if anything differs."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E82_protocol_runs as R
    path = out_dir() / 'identity.json'
    res = _read_json(path)
    res.update(_commits_clean())
    orn = lab_input()
    res.setdefault('units', {})
    for u in units:
        if u in res['units'] and 'count_diffs' in res['units'][u]:
            continue
        t0 = time.time()
        F = fit_of(u)
        try:
            if is_fc(u):
                m = fc_parts(u)
                sim, info = fc_build(u, orn)
                own = R.own_breadth(m)
            else:
                sim, info = E85.build(u, orn, cache_tag(orn, 'I0'), cache_dir=local_dir())
                own = E83.wiring_parts(u)['breadth']
        except RuntimeError as e:
            res['units'][u] = dict(error=str(e), count_diffs=-1, response_diffs=-1, breadth_equal=False)
            _write_json(path, res)
            continue
        b = info['breadth']
        beq = bool(b.index.equals(own.index) and np.array_equal(b.to_numpy(dtype=float), own.to_numpy(dtype=float)))
        on = sim(dict(F['p']))
        off = sim(dict(F['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
        ns, nr = E83._same({'on': on, 'off': off}, F['saved'])
        res['units'][u] = dict(fit=f"{F['folder']}:{F['tag']}", breadth_equal=beq, count_diffs=ns, response_diffs=nr,
                               seconds=round(time.time() - t0, 1))
        _write_json(path, res)
        print(f'E108 0 {u}: {res["units"][u]}', flush=True)
        del sim, info, on, off
        gc.collect()
    want = list(ALL_WIRINGS + FC_WIRINGS)
    done = [res['units'][u] for u in want if u in res['units']]
    res['all_exact'] = bool(len(done) == len(want) and all(d.get('breadth_equal') and d['count_diffs'] == 0
                                                           and d['response_diffs'] == 0 for d in done))
    _write_json(path, res)
    bad = {u: d for u, d in res['units'].items() if not (d.get('breadth_equal') and d['count_diffs'] == 0 and d['response_diffs'] == 0)}
    if bad:
        raise RuntimeError(f'E108 part 0: the lab input does not reproduce the saved fits: {bad} (the experiment stops here)')
    return res


def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('all_exact'))


# ---------------------------------------------------------------------------
# part B

def run_tag(name: str, unit: str) -> str:
    return f'B_{name}_{unit}'


def b_plan() -> List[Tuple[str, str]]:
    """[(input, unit)] of part B in order: the rules, then the controls, on the eight wirings; the same inputs on the
    ceiling family's four; step 4's flies on the eight."""
    plan = [(n, w) for n in RULES for w in ALL_WIRINGS]
    plan += [(n, w) for n in CONTROLS for w in ALL_WIRINGS]
    plan += [(n, u) for n in RULES + CONTROLS for u in FC_WIRINGS]
    plan += [(n, w) for n in BMYR_INPUTS for w in ALL_WIRINGS]
    return plan


def lane_plan(lane: Optional[str]) -> List[Tuple[str, str]]:
    """The plan's points of a lane: 'a' the even positions, 'b' the odd ones, None all."""
    plan = b_plan()
    if lane is None:
        return plan
    i = {'a': 0, 'b': 1}[lane]
    return plan[i::2]


def held(unit: str, orn: pd.DataFrame, name: str) -> dict:
    """One held reading: hold_rates' dict ('on', 'off', 'params', 'params_off', 'thr_factor', 'apl_factor', ...)."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E77_input_screen as S
    tag = cache_tag(orn, name)
    if not is_fc(unit):
        return E85.held_reading(unit, orn, tag, cache_dir=local_dir())
    F = fit_of(unit)
    res = S.hold_rates(lambda: fc_build(unit, orn), target_off=F['target_off'], target_on=F['target_on'], params=F['p'],
                       label=f'{unit} {tag}')
    res['input_tag'] = tag
    return res


def save_run(tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir('B')
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}): this host's outputs, else the model-output archive; or None."""
    for d in (H.outputs_dir(EXPERIMENT) / 'B', H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'B'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_outputs(name: str, unit: str) -> Dict[str, tuple]:
    """A run of the readings: I0 is the unit's saved fit; the others part B's."""
    if name == 'I0':
        return fit_of(unit)['saved']
    o = load_run(run_tag(name, unit))
    if o is None:
        raise RuntimeError(f'no saved run {run_tag(name, unit)}')
    return o


def run_extras(o: Dict[str, tuple]) -> dict:
    """paa's and B-myr's response fractions (two or more spikes), APL on and off."""
    row = {}
    for c in ('on', 'off'):
        r = H.abbrev_columns(o[c][0].astype(float) > 0.5)
        for od in ('paa', BMYR):
            row[f'frac_{od}_{c}'] = float(r[od].mean())
    return row


def run_row(name: str, unit: str, tag: str, res: dict, o: Dict[str, tuple], rec: dict, cc: dict) -> dict:
    F = fit_of(unit)
    row = dict(part='B', input=name, unit=unit, family='fc' if is_fc(unit) else 'C2', tag=tag, fit_folder=F['folder'],
               fit_tag=F['tag'], input_file_sha256=rec['file_sha256'], input_frame_sha256=rec['frame_sha256'],
               input_tag=res['input_tag'], orn_input=rec['file'], thr_factor=res['thr_factor'], apl_factor=res['apl_factor'],
               n_on=res['n_on'], n_off=res['n_off'], target_on=res['target_on'], target_off=res['target_off'],
               at_bound_thr=res['at_bound_thr'], at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'],
               seconds=res.get('seconds'), **{k: res['params'][k] for k in E83.PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}),
               K_run=E83.readout_K(o), breadth_sha256=res.get('breadth_sha256'))
    row.update(E83.realism(o))
    row.update(run_extras(o))
    row.update(cc)
    return row


def _read_rows(p: Path, tries: int = 5) -> pd.DataFrame:
    """A lane's table, read again after a moment if another lane is appending to it."""
    for i in range(tries):
        try:
            return pd.read_csv(p)
        except (pd.errors.ParserError, pd.errors.EmptyDataError):
            if i == tries - 1:
                raise
            time.sleep(2.0)
    raise AssertionError


def runs_b() -> pd.DataFrame:
    """Every lane's rows (runs_B_<lane>.csv)."""
    ps = sorted(out_dir().glob('runs_B_*.csv'))
    if not ps:
        return pd.DataFrame(columns=['tag', 'input', 'unit'])
    return pd.concat([_read_rows(p) for p in ps], ignore_index=True, sort=False).drop_duplicates('tag', keep='last')


def part_b(lane: Optional[str] = None, max_points: Optional[int] = None) -> pd.DataFrame:
    """Part B (the module docstring): the lane's points in order, resumable (a point is done when its row is in any
    lane's table). Needs part 0 to have passed and clean checkouts."""
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identity must hold before any reading)')
    cc = _commits_clean()
    path = out_dir() / f'runs_B_{lane or "all"}.csv'
    prov_path = out_dir() / f'provenance_B_{lane or "all"}.json'
    n = 0
    for name, unit in lane_plan(lane):
        tag = run_tag(name, unit)
        if tag in set(runs_b()['tag']) or (max_points is not None and n >= max_points):
            continue
        orn = input_of(name)
        rec = write_input(name)
        res = held(unit, orn, name)
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = run_row(name, unit, tag, res, o, rec, cc)
        E83._append_row(path, row)
        prov = _read_json(prov_path) or dict(part='B', lane=lane, script='experiments/E108_paa_extrapolation.py', points={})
        prov['points'][tag] = H.model_provenance(orn_deltas=orn)
        _write_json(prov_path, prov)
        print(f"E108 B {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, "
              f"{res['n_off']}/{res['target_off']} off); paa {row['frac_paa_on']:.3f} / {row['frac_paa_off']:.3f}; "
              f"flag '{row.get('realism_flag') or ''}'; {res['seconds']:.0f}s", flush=True)
        n += 1
        del res, o
        gc.collect()
    return runs_b()


# ---------------------------------------------------------------------------
# scoring

def dist_targets() -> Dict[str, Dict[str, pd.Series]]:
    """The KC data's mean distances per condition ({'on', 'off'}, canonical order): 'all' (all recordings, E77's 'z8')
    and the halves A to D (E77's `recording_halves(0)` A / B, `recording_halves(1)` A / B as C / D)."""
    if 'dist_targets' not in _CACHE:
        from experiments import E77_input_side_common_error as E77
        ctrl = H.yang_recording_rdms_values('TNT3c', 'graded')
        tnt = H.yang_recording_rdms_values('TNT', 'graded')
        out = {'all': {c: E83._cs(v) for c, v in E77.data_distances('z8').items()}}
        for seed, names in ((0, ('A', 'B')), (1, ('C', 'D'))):
            hs = E77.recording_halves(seed)
            for k, h in zip(('A', 'B'), names):
                out[h] = {'on': E83._cs(H.canonical_pairs(ctrl.loc[hs[k]['ctrl']].mean())),
                          'off': E83._cs(H.canonical_pairs(tnt.loc[hs[k]['tnt']].mean()))}
        _CACHE['dist_targets'] = out
    return _CACHE['dist_targets']


def effect_targets() -> Dict[str, pd.Series]:
    """The per-pair APL effect targets: 'all' (E83's 'z8') and the halves A to D (E89's `targets89`)."""
    if 'effect_targets' not in _CACHE:
        from experiments import E89_apl_breadth_weighting as E89
        tg = E89.targets89()
        _CACHE['effect_targets'] = {'all': tg['z8'], **{h: tg[f'half{h}'] for h in HALVES}}
    return _CACHE['effect_targets']


def rho_pairs_with(model: pd.Series, data: pd.Series, pairs: Optional[Sequence[Tuple[str, str]]] = None) -> float:
    """ρ on odor-odor distances of `model` against `data` over `pairs` (all both have when None; E106's `rho_pairs`)."""
    from experiments import E106_rho_gap_parts as E106
    return float(E106.rho_pairs(model, data, pairs=pairs))


def unit_scores(d: Dict[str, pd.Series]) -> dict:
    """The readings of one unit's distances ({'on', 'off'}; a run's or a population's): the two ρ's against every
    distance target over the 120 pairs ('rho_<on|off>_<target>') and paa's 15 ('rho15_...'); against every effect
    target the loss, the pattern part, its amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair
    APL effect ('<target>_<stat>'); the mean APL effect and the mean distances; the block, the natural-source alignment
    on and off and the structure vs wild-type somata (on and off)."""
    on, off = E83._cs(d['on']), E83._cs(d['off'])
    row = dict(mean_on=float(on.mean()), mean_off=float(off.mean()))
    for t, dd in dist_targets().items():
        for c, m in (('on', on), ('off', off)):
            row[f'rho_{c}_{t}'] = rho_pairs_with(m, dd[c])
            row[f'rho15_{c}_{t}'] = rho_pairs_with(m, dd[c], PAA15)
    me = E83._cs(H.apl_effect(on, off))
    row['mean_effect'] = float(me.mean())
    for t, de0 in effect_targets().items():
        de = E83._cs(de0).reindex(me.index)
        err = (me - me.mean()) - (de - de.mean())
        pattern = _rms(err.to_numpy())
        sp = H.apl_effect_spread_shape(on, off, H.canonical_pairs(de))
        if not np.isclose(pattern ** 2, sp['apl_effect_amplitude'] ** 2 + sp['apl_effect_shape'] ** 2, rtol=1e-9, atol=1e-12):
            raise RuntimeError(f'{t}: the amplitude and shape parts do not add up to the pattern part (pairs misaligned?)')
        row.update({f'{t}_loss': _rms((me - de).to_numpy()), f'{t}_pattern': pattern,
                    f'{t}_amplitude': sp['apl_effect_amplitude'], f'{t}_shape': sp['apl_effect_shape'],
                    f'{t}_sd_ratio': sp['apl_effect_sd_ratio'], f'{t}_rho_effect': H.spearman(me, de),
                    f'{t}_pearson_effect': float(np.corrcoef(me.to_numpy(), de.to_numpy())[0, 1])})
    st, _ = H.compare_on_off(on, off)
    row['block'] = float(st['block_fruity_small_iso'])
    for c, m in (('on', on), ('off', off)):
        r = H.compare_to_references(m)
        row[f'vcf_{c}'] = float(r['spearman_resid_vs_vcf_hyp'])
        row[f'struct_wt_{c}'] = float(r['spearman_resid_vs_remap_wt'])
    return row


def _dist(name: str, unit: str, readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83's `pair_distances`; 'cal' with K set on the run's APL-on counts), cached in this
    process and in outputs/.../distances_cache.pkl by the sha256 of the run's spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    mem = _CACHE.setdefault('dist_mem', {})
    key = (name, unit, readout)
    if key in mem:
        return mem[key]
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
    store = _CACHE['dist_store']
    o = run_outputs(name, unit)
    skey = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if skey not in store:
        store[skey] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    mem[key] = store[skey]
    return mem[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist_store']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def runs_table() -> pd.DataFrame:
    """Every run: the saved fits as I0 (eight wirings and the ceiling family's four, with their realism) and part B's,
    with 'flagged'."""
    rows = []
    for u in ALL_WIRINGS + FC_WIRINGS:
        F = fit_of(u)
        rows.append(dict(tag=f'I0_{u}', input='I0', unit=u, family='fc' if is_fc(u) else 'C2', thr_factor=1.0, apl_factor=1.0,
                         at_bound_thr=False, at_bound_apl=False, fit_folder=F['folder'], fit_tag=F['tag'],
                         **E83.realism(F['saved']), **run_extras(F['saved'])))
    rb = runs_b()
    df = pd.concat([pd.DataFrame(rows), rb], ignore_index=True, sort=False) if len(rb) else pd.DataFrame(rows)
    df = df.drop_duplicates(['input', 'unit'], keep='last').reset_index(drop=True)
    df['flagged'] = [_flagged(x) for x in df['realism_flag']]
    df['at_bound'] = df['at_bound_thr'].fillna(False).astype(bool) | df['at_bound_apl'].fillna(False).astype(bool)
    return df


RUN_MEANS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'rate_rho_on_zcrit', 'rate_rho_off_zcrit',
             'opk_on_tv', 'opk_off_tv', 'frac_paa_on', 'frac_paa_off', 'frac_B-myr_on', 'frac_B-myr_off', 'max_count_on',
             'max_count_off', 'thr_factor', 'apl_factor')


def readings() -> pd.DataFrame:
    """readings.csv (the module docstring), with 'd_<column>' the change from I0's same unit and readout; and runs.csv."""
    out = out_dir()
    runs = runs_table()
    runs.to_csv(out / 'runs.csv', index=False)
    rows = []
    for name, g in runs.groupby('input', sort=False):
        have = set(g['unit'])
        meta = {r.unit: r for r in g.itertuples()}
        units = [(u, (u,)) for u in g['unit']] + [(u, us) for u, us in UNITS if all(x in have for x in us)]
        for ro in READOUTS:
            dd = {u: _dist(name, u, ro) for u in g['unit']}
            for unit, us in units:
                single = len(us) == 1 and unit in dd
                d = dd[unit] if single else E83.population([dd[u] for u in us])
                sub = g[g['unit'].isin(us)]
                rows.append(dict(input=name, unit=unit, kind='wiring' if single else 'population', n_wirings=len(us),
                                 family='fc' if all(is_fc(u) for u in us) else 'C2', readout=ro,
                                 n_flagged=int(sum(meta[u].flagged for u in us)),
                                 n_at_bound=int(sum(meta[u].at_bound for u in us)),
                                 **{f'run_{k}': float(sub[k].astype(float).mean()) for k in RUN_MEANS if k in sub},
                                 **unit_scores(d)))
    _save_dist_cache()
    rd = pd.DataFrame(rows)
    base = rd[rd['input'] == 'I0'].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c not in ('input', 'unit', 'kind', 'n_wirings', 'family', 'readout')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, r), np.nan) for u, r, v in zip(rd['unit'], rd['readout'], rd[c])]
    rd.to_csv(out / 'readings.csv', index=False)
    return rd


def _v(rd: pd.DataFrame, name: str, unit: str, col: str, readout: str) -> float:
    r = rd[(rd['input'] == name) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _have(runs: pd.DataFrame, name: str, units: Sequence[str]) -> bool:
    got = set(zip(runs['input'], runs['unit']))
    return name == 'I0' or all((name, u) in got for u in units)


def wirings_needed(n: int) -> int:
    """The wirings that must gain: GAIN_WIRINGS of the family's 8; the same share, rounded up, of another population
    (3 of the ceiling family's 4, which is read beside)."""
    return GAIN_WIRINGS if n == len(ALL_WIRINGS) else int(np.ceil(GAIN_WIRINGS / len(ALL_WIRINGS) * n))


def gate(rd: pd.DataFrame, runs: pd.DataFrame, name: str, ref: str = 'I0', units: Sequence[str] = ALL_WIRINGS,
         pop: str = 'all8') -> dict:
    """The gate of the entry (the module docstring) for input `name` against `ref` on the population `pop` of `units`
    (default the family's eight). `rd`: `readings`' layout; `runs`: `runs_table`'s ('input', 'unit', 'flagged')."""
    out = dict(input=name, reference=ref, population=pop, complete=bool(_have(runs, name, units)))
    if not out['complete']:
        out.update(gain=False, note='runs missing')
        return out
    D = lambda col, ro='raw', u=pop: _v(rd, name, u, col, ro) - _v(rd, ref, u, col, ro)
    other = {'on': 'off', 'off': 'on'}
    gains = {}
    for c in ('on', 'off'):
        d_all = D(f'rho_{c}_all')
        n_w = int(sum(D(f'rho_{c}_all', 'raw', w) >= GAIN for w in units))
        d_h = {h: D(f'rho_{c}_{h}') for h in HALVES}
        ok = bool(d_all >= GAIN and n_w >= wirings_needed(len(units)) and all(v >= GAIN for v in d_h.values()))
        out.update({f'd_rho_{c}_all': d_all, f'n_wirings_gain_{c}': n_w, **{f'd_rho_{c}_{h}': v for h, v in d_h.items()},
                    f'd_rho15_{c}_all': D(f'rho15_{c}_all')})
        gains[c] = ok
    for c in ('on', 'off'):
        out[f'gain_{c}'] = bool(gains[c] and out[f'd_rho_{other[c]}_all'] >= -OTHER_SLACK)
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    # the APL effect (E94's equal), under the calibrated readout; the raw-count loss
    fails = {}
    for t in TARGETS:
        out[f'd_pattern_{t}'] = D(f'{t}_pattern', 'cal')
        if out[f'd_pattern_{t}'] > EQUAL_PATTERN:
            fails[f'pattern part ({t})'] = out[f'd_pattern_{t}']
    out['d_raw_loss'] = D('all_loss', 'raw')
    if out['d_raw_loss'] > EQUAL_RAW:
        fails['raw-count loss'] = out['d_raw_loss']
    out['d_mean_effect'] = D('mean_effect', 'cal')
    if out['d_mean_effect'] < -EQUAL_MEAN:
        fails['mean APL effect'] = out['d_mean_effect']
    for k in ('amplitude', 'shape', 'sd_ratio', 'rho_effect', 'pearson_effect', 'loss'):
        out[f'd_{k}_all_cal'] = D(f'all_{k}', 'cal')
    out['apl_effect_ok'] = not any(k.startswith(('pattern', 'raw', 'mean')) for k in fails)
    # realism
    r = runs[(runs['input'] == name) & runs['unit'].isin(units)]
    out['n_flagged'] = int(r['flagged'].sum())
    if out['n_flagged']:
        fails['realism (runs outside the bands)'] = out['n_flagged']
    out['realism_ok'] = out['n_flagged'] == 0
    # response rates (the mean over the wirings)
    rr = runs[runs['input'] == ref]
    rr = rr[rr['unit'].isin(units)]
    for c in ('on', 'off'):
        v = float(r[f'rate_rho_{c}'].astype(float).mean()) - float(rr[f'rate_rho_{c}'].astype(float).mean())
        out[f'd_rate_rho_{c}'] = v
        if v < -RATE_RHO_SLACK:
            fails[f'response rate ρ, APL {c}'] = v
        out[f'frac_paa_{c}'] = float(r[f'frac_paa_{c}'].astype(float).mean())
        out[f'frac_paa_{c}_ref'] = float(rr[f'frac_paa_{c}'].astype(float).mean())
        out[f'frac_paa_{c}_flies'] = FLIES_PAA_FRAC[c]
    out['rates_ok'] = not any(k.startswith('response rate') for k in fails)
    # pending Tom's answers of 2026-11-17 14:44
    pend = {}
    out['d_block'] = D('block')
    if out['d_block'] < -BLOCK_SLACK:
        pend['the block'] = out['d_block']
    for c in ('on', 'off'):
        out[f'd_vcf_{c}'] = D(f'vcf_{c}')
        if out[f'd_vcf_{c}'] < -STRUCTURE_SLACK:
            pend[f'natural-source alignment, APL {c}'] = out[f'd_vcf_{c}']
    out['d_struct_wt_on'] = D('struct_wt_on')
    out['d_struct_wt_off'] = D('struct_wt_off')
    if out['d_struct_wt_on'] < -STRUCTURE_SLACK:
        pend['structure vs wild-type somata'] = out['d_struct_wt_on']
    out['pending_ok'] = not pend
    # controls (X1's, the entry's), on each ρ that gains
    ctl = {}
    for c in ('on', 'off'):
        if not out[f'gain_{c}']:
            continue
        g0 = out[f'd_rho_{c}_all']
        w5 = _v(rd, 'W5', pop, f'rho_{c}_all', 'raw') - _v(rd, ref, pop, f'rho_{c}_all', 'raw')
        pg = [_v(rd, p, pop, f'rho_{c}_all', 'raw') - _v(rd, ref, pop, f'rho_{c}_all', 'raw') for p in PERMS if _have(runs, p, units)]
        pm = float(np.mean(pg)) if len(pg) == len(PERMS) else np.nan
        ctl[c] = dict(gain=g0, W5=w5, perm_mean=pm, perm=pg,
                      ok=bool(np.isfinite(w5) and np.isfinite(pm) and w5 < CONTROL_RATIO * g0 and pm < CONTROL_RATIO * g0))
        if not ctl[c]['ok']:
            fails[f'controls, APL {c} (W5 {w5:+.3f}, permuted mean {pm:+.3f} against {g0:+.3f})'] = g0
    out['controls'] = ctl
    out['controls_ok'] = all(v['ok'] for v in ctl.values()) if ctl else None
    out['fails'] = fails
    out['pending_fails'] = pend
    w5_lose = [_v(rd, 'W5', pop, f'rho_{c}_all', 'raw') - _v(rd, ref, pop, f'rho_{c}_all', 'raw') for c in ('on', 'off')]
    out['W5_d_rho'] = w5_lose
    out['W5_loses'] = bool(np.nanmin(w5_lose) <= -GAIN) if np.isfinite(w5_lose).any() else None
    return out


def class_of(g: dict, with_pending: bool = True) -> str:
    """'better' (the gain and every constraint hold), 'a price' (the gain holds and a constraint fails; the pending
    items count unless `with_pending` is False) or 'no gain' (`gate`'s dict)."""
    if not g.get('gain'):
        return 'no gain'
    failed = bool(g['fails']) or (with_pending and bool(g['pending_fails']))
    return 'a price' if failed else 'better'


def outcome_of(g: dict) -> str:
    """The entry's outcome letter: (a) better, (b) a price, (c) no gain and W5 loses, (d) no gain and W5 loses
    nothing."""
    c = class_of(g)
    if c == 'better':
        return 'a'
    if c == 'a price':
        return 'b'
    return 'c' if g.get('W5_loses') else 'd'


def gate_tables(rd: pd.DataFrame, runs: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """gate.csv (one row per rule and control, the family's eight; the ceiling family's four beside, 'fc4') and
    gate.json (with the carrying rule's verdict and class)."""
    rule = _read_json(out_dir() / 'rule.json').get('carrying_rule')
    rows, full = [], {}
    for pop, units in (('all8', ALL_WIRINGS), ('fc4', FC_WIRINGS)):
        for name in RULES + CONTROLS:
            g = gate(rd, runs, name, units=units, pop=pop)
            full[f'{pop}:{name}'] = g
            flat = {k: v for k, v in g.items() if not isinstance(v, (dict, list))}
            flat['fails'] = '; '.join(f'{k} {v:+.4f}' if isinstance(v, float) else f'{k} {v}' for k, v in g.get('fails', {}).items())
            flat['pending_fails'] = '; '.join(f'{k} {v:+.4f}' for k, v in g.get('pending_fails', {}).items())
            if g.get('complete'):
                flat['class'] = class_of(g)
                flat['class_without_pending'] = class_of(g, with_pending=False)
                flat['outcome'] = outcome_of(g)
            flat['carrying'] = name == rule
            rows.append(flat)
    gt = pd.DataFrame(rows)
    gt.to_csv(out_dir() / 'gate.csv', index=False)
    v = full.get(f'all8:{rule}', {})
    verdict = dict(carrying_rule=rule, complete=bool(v.get('complete')), status=E106_PENDING)
    if v.get('complete'):
        verdict.update(gain=v['gain'], cls=class_of(v), class_without_pending=class_of(v, False), outcome=outcome_of(v),
                       fails=v['fails'], pending_fails=v['pending_fails'])
    _write_json(out_dir() / 'gate.json', dict(verdict=verdict, gates=full, **H.code_commits()))
    return gt, verdict


def class_table(gt: pd.DataFrame) -> pd.DataFrame:
    """classes.csv: per rule (the family's eight), its class, the class without the pending items, the outcome and the
    column "provisional, E106 pending" (the class is written after E106's note)."""
    g = gt[(gt['population'] == 'all8') & gt['input'].isin(RULES)].copy()
    cols = [c for c in ('input', 'carrying', 'complete', 'gain', 'class', 'class_without_pending', 'outcome', 'fails',
                        'pending_fails') if c in g.columns]
    ct = g[cols].copy()
    ct[E106_PENDING] = True
    ct.to_csv(out_dir() / 'classes.csv', index=False)
    return ct


def bmyr_table(rd: pd.DataFrame) -> pd.DataFrame:
    """Step 4 (beside, no gate): per B-myr fly, the family's two ρ's (population of eight, raw counts) over the 120 pairs
    and over B-myr's 15, against all recordings, with I0's beside and the range over the flies."""
    from experiments import E106_rho_gap_parts as E106
    rows = []
    bm15 = tuple((BMYR, o) for o in H.APL16_ODORS if o != BMYR)
    dd = dist_targets()['all']
    for name in ('I0',) + BMYR_INPUTS:
        r = rd[(rd['input'] == name) & (rd['unit'] == 'all8') & (rd['readout'] == 'raw')]
        if not len(r):
            continue
        d = E83.population([_dist(name, w, 'raw') for w in ALL_WIRINGS])
        rows.append(dict(input=name, rho_on=float(r['rho_on_all'].iloc[0]), rho_off=float(r['rho_off_all'].iloc[0]),
                         rho15_bmyr_on=float(E106.rho_pairs(d['on'], dd['on'], pairs=bm15)),
                         rho15_bmyr_off=float(E106.rho_pairs(d['off'], dd['off'], pairs=bm15)),
                         pattern_cal=_v(rd, name, 'all8', 'all_pattern', 'cal'), mean_effect_cal=_v(rd, name, 'all8', 'mean_effect', 'cal'),
                         n_flagged=int(r['n_flagged'].iloc[0])))
    bt = pd.DataFrame(rows)
    if len(bt) > 1:
        f = bt[bt['input'] != 'I0']
        bt = pd.concat([bt, pd.DataFrame([dict(input='range over flies', **{c: f'{f[c].min():.3f} to {f[c].max():.3f}'
                                                                             for c in ('rho_on', 'rho_off', 'rho15_bmyr_on', 'rho15_bmyr_off')})])],
                       ignore_index=True)
    bt.to_csv(out_dir() / 'bmyr.csv', index=False)
    return bt


def _f(x, spec='.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


def comparison_table(rd: pd.DataFrame, gt: pd.DataFrame, verdict: dict) -> None:
    """comparison_table.md: the lab's input, each rule and each control on the family's eight (population), the
    statistics of the gate, and the gate's verdict per input."""
    sc = pd.read_csv(out_dir() / 'scores.csv', index_col=0) if (out_dir() / 'scores.csv').exists() else None
    L = ['# E108: paa\'s input at the flies\' concentration', '',
         'paa\'s input at 10^-3 is an extrapolation, not a recording; euc and -aPine stay at the recorded 10^-1.75 and 10^-1.5 '
         '(the flies\' 10^-2). Status: ' + E106_PENDING + '.', '']
    if sc is not None:
        L += ['## The rule (RMSE of the predicted top concentration on the ramps\' triples, estimated spike deltas)', '',
              '| rule | RMSE | bias | triples |', '|---|---|---|---|']
        L += [f'| {r} | {_f(sc.loc[r, "rmse"])} | {_f(sc.loc[r, "bias"], "+.3f")} | {int(sc.loc[r, "n_triples"])} |' for r in sc.index]
        L += ['', f'Carrying rule: {verdict.get("carrying_rule")}.', '']
    cols = [('ρ on odor-odor distances, APL on (raw)', 'rho_on_all', 'raw'), ('ρ on odor-odor distances, APL off (raw)', 'rho_off_all', 'raw'),
            ('ρ, paa\'s 15 pairs, APL on (raw)', 'rho15_on_all', 'raw'), ('ρ, paa\'s 15 pairs, APL off (raw)', 'rho15_off_all', 'raw'),
            ('APL-effect loss (readout)', 'all_loss', 'cal'), ('pattern part (readout)', 'all_pattern', 'cal'),
            ('amplitude part (readout)', 'all_amplitude', 'cal'), ('shape part (readout)', 'all_shape', 'cal'),
            ('spread ratio (readout)', 'all_sd_ratio', 'cal'), ('mean APL effect (readout)', 'mean_effect', 'cal'),
            ('APL-effect loss (raw)', 'all_loss', 'raw'), ('pattern part (raw)', 'all_pattern', 'raw'),
            ('mean APL effect (raw)', 'mean_effect', 'raw'), ('ρ on the per-pair APL effect (readout)', 'all_rho_effect', 'cal'),
            ('response rate ρ, APL on', 'run_rate_rho_on', 'raw'), ('response rate ρ, APL off', 'run_rate_rho_off', 'raw'),
            ('paa\'s response fraction, APL on', 'run_frac_paa_on', 'raw'), ('paa\'s response fraction, APL off', 'run_frac_paa_off', 'raw'),
            ('the block (raw)', 'block', 'raw'), ('natural-source alignment, APL on (raw)', 'vcf_on', 'raw'),
            ('natural-source alignment, APL off (raw)', 'vcf_off', 'raw'), ('structure vs wild-type somata (raw)', 'struct_wt_on', 'raw'),
            ('runs outside the bands', 'n_flagged', 'raw')]
    for pop, title in (('all8', 'the family, mean matrix over its eight wirings'), ('fc4', 'the ceiling family, four wirings (beside)')):
        names = [n for n in INPUTS if len(rd[(rd['input'] == n) & (rd['unit'] == pop)])]
        if not names:
            continue
        L += [f'## {title}', '', '| statistic | ' + ' | '.join(names) + ' |', '|---|' + '---|' * len(names)]
        for lab, col, ro in cols:
            L.append(f'| {lab} | ' + ' | '.join(_f(_v(rd, n, pop, col, ro), '.0f' if col == 'n_flagged' else '.3f') for n in names) + ' |')
        L.append('')
    g = gt[gt['population'] == 'all8']
    if len(g):
        L += ['## The gate (the family\'s eight, against the lab\'s input)', '',
              '| input | gain on / off | wirings ≥ 0.010 on / off | APL effect | realism | response rates | pending | controls | class | outcome |',
              '|---|---|---|---|---|---|---|---|---|---|']
        ok = lambda x: '-' if x is None or (isinstance(x, float) and np.isnan(x)) else ('ok' if x else 'fails')
        for r in g.to_dict('records'):
            if not r.get('complete'):
                L.append(f'| {r["input"]} | runs missing | | | | | | | | |')
                continue
            L.append(f'| {r["input"]}{" (carrying)" if r["carrying"] else ""} | {_f(r["d_rho_on_all"], "+.3f")} / '
                     f'{_f(r["d_rho_off_all"], "+.3f")} | {int(r["n_wirings_gain_on"])} / {int(r["n_wirings_gain_off"])} | '
                     f'{ok(r["apl_effect_ok"])} | {ok(r["realism_ok"])} | {ok(r["rates_ok"])} | {ok(r["pending_ok"])} | '
                     f'{ok(r.get("controls_ok"))} | {r.get("class", "")} | {r.get("outcome", "")} |')
        L += ['', f'Verdict ({E106_PENDING}): ' + json.dumps({k: v for k, v in verdict.items() if k not in ('fails', 'pending_fails')}), '']
    (out_dir() / 'comparison_table.md').write_text('\n'.join(L) + '\n')


def score() -> dict:
    """Part score (no simulations)."""
    rd = readings()
    runs = runs_table()
    gt, verdict = gate_tables(rd, runs)
    class_table(gt)
    bmyr_table(rd)
    comparison_table(rd, gt, verdict)
    _write_json(out_dir() / 'provenance_score.json', dict(part='score', script='experiments/E108_paa_extrapolation.py',
                                                          **H.code_commits()))
    print(json.dumps(verdict, indent=1, default=str), flush=True)
    return verdict


def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM), one row per part-B run (E93's layout)."""
    out = out_dir()
    rows, parts = [], []
    for _, r in runs_b().iterrows():
        fam = ("the firing-ceiling family's" if r.family == 'fc' else "C2's family's") + f' {r.unit}'
        rows.append(dict(tag=r.tag, variant=r.tag, part='B', arm=r.input, base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         kind=(f'{fam} fit with the ORN input {r.input} (results/{EXPERIMENT}/inputs/{r.input}.csv: the lab\'s input '
                               + ('with paa\'s column extrapolated or a control' if r.input in RULES + CONTROLS
                                  else 'with one validation2 fly\'s B-myr') +
                               '); breadth recomputed for the input; the somatic threshold and the tied APL scale re-held at the '
                               'fit\'s APL-off and APL-on counts'),
                         config=r.input, wiring=r.unit, family=r.family, ref_experiment=r.fit_folder, ref_fit=r.fit_tag,
                         orn_input=r.orn_input, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, thr_factor=r.thr_factor, apl_factor=r.apl_factor, rate_on=r.rate_on,
                         rate_off=r.rate_off, max_count_on=r.max_count_on, max_count_off=r.max_count_off,
                         realism_flag=r.realism_flag, subdir='B'))
        o = load_run(r.tag)
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('inputs', '0', 'B', 'score', 'variants'))
    ap.add_argument('--lane', default=None, choices=('a', 'b'))
    ap.add_argument('--max-points', type=int, default=None)
    a = ap.parse_args()
    if a.part == 'inputs':
        write_inputs()
    elif a.part == '0':
        identity()
    elif a.part == 'B':
        part_b(a.lane, a.max_points)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()


if __name__ == '__main__':
    main()
