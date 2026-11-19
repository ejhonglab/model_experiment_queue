#!/usr/bin/env python3
"""E112: the two ρ's on odor-odor distances read as the recordings are made: the best family's saved runs through the
calibrated readout, the per-KC session z-score and the trial noise; the gap, its bound and the gates of the input
experiments read again (no new simulations in part A).

The entry is item E112 of experiment_priorities.md (proposed by viki 2026-11-19 13:58; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 13:50; claimed by hal 14:01). "The family" is the breadth-2 APL-effect fit (C2, E67)
on its eight wirings (E72's part B; E83's `saved_fit`: S = real, P1, P2, M1; T = P3, P4, M2, M3), "the two ρ's" ρ on
odor-odor distances with APL on against the control flies and with APL off against APL output blocked (TNT), each KC
matrix the mean over recordings of the cosine distance (E106's `data_matrices`), through E106's `rho_pairs`.

The readings (`reading_distances`), each applied to saved spike counts with existing functions, both conditions read
the same way (calibrated on APL on only, so the APL effect stays a prediction):
- R0, raw counts (cosine on spike counts; E83's `pair_distances(o, 'raw')`, the present reading);
- R1, E65's calibrated saturating readout, K set per run on its APL-on counts (E83's `pair_distances(o, 'cal')`);
- R2, E73's per-KC session z-score with the saturating map, graded, in expectation over the measurement noise,
  calibrated per run on its APL-on counts as E74 / E93 do (`zcell_readings`; equal to E93's `zcell_distances`);
- R3, R2 through sampled measurement noise (`pseudo_recording_distances`): pseudo-recordings of three trials, 11 with APL
  on and 8 with APL off as in the flies (`N_PSEUDO`), each trial the reading's noise-free value plus Gaussian noise of the
  KC's trial sd, the cosine distance of each pseudo-recording's trial means averaged over the pseudo-recordings as the
  data's are; seeded per wiring (`wiring_seed`: the same draws for every run on a wiring, so changes are paired);
- R4, the threshold-jitter trial noise (E43's part B pipeline at E72's part E settings: seed 43, 10 pseudo-recordings
  of three trials, the noise matched to the flies' trial-to-trial variability), raw counts: the mean over
  pseudo-recordings of the distances of their three-trial means, APL on at the level matched to the control flies and
  APL off at the level matched to APL output blocked (TNT) (`trial_noise_distances`, from partB/pairs_noise_long.csv);
- beside, R4z: R2 applied to the same jittered trials (E73's `_trials` reading: calibrated on the APL-on trials), where
  the trials are saved (`trial_zcell_distances`).

Part A (`--part A`, no simulations):
1. per reading, per wiring and per population (the eight, S, T; R4 on the four wirings E72's part E ran, which are S),
   the scores (`reading_scores`): the two ρ's over the 120, 66, 54 and 91 pairs, the distance of paa from B-myr, the
   APL-effect loss with its mean, pattern, amplitude and shape parts, ρ and Pearson r on the per-pair APL effect, the
   mean distances (rhos.csv); E106's bound with its odor part and pair part on each population (`bound_scores`,
   bound.csv); the realism of each wiring under each reading's criterion (realism.csv);
2. the verdict "the reading moves ρ" (`moves_verdict`, verdict.json), with the per-wiring changes (changes.csv);
3. the yardstick of a reading (yardstick.csv): S's population against T's under R0, R1, R2 and R3 (under R4 once part
   B has run; in part A, the two halves of S beside), the data's split-half ceilings and ρ's, and the family against
   the flies' noise-free distances (E81) under each reading;
4. the saved gates read again (`reread_gate`): E109's held readings, E108's runs and E97's points under R1 to R3, their
   classes as written beside the re-read ones (gates.csv, gates.json) and the verdict "E109's gain holds like for like"
   (`e109_holds`);
5. the identities (identities.json; the run stops if one fails): R0 and R1 reproduce E106's and E72's population ρ's, R2
   reproduces E73's reading of C2 and E93's `zcell_distances`, the pipeline's deterministic distances are R0, R4 of the
   four reproduces E72's partE_noise.csv, the gates rebuilt under R0 and R1 reproduce each experiment's saved
   readings.csv and classes, the data matrices are E73's graded ones.

Part B (simulations): `--part B-prepare` writes results/E112_rho_like_for_like/fits.csv (one row per E109 input on each
held-out wiring: the wiring's E72 row with the input file and the held reading's absolute parameters) and links each
run's saved outputs under outputs/E112_rho_like_for_like/ so that E43's pipeline asserts its deterministic run against
them; then E43's part B runs per wiring (commands in `part_b_commands`, printed by `--part B-prepare`); `--part
B-score` scores R4 on all eight wirings and E109's inputs and rewrites the verdicts.

Readings chosen where the entry is open (fixed here before any reading was scored):
- "the population's ρ is higher than under R0 by 0.020 or more in at least one condition, on 6 of 8 wirings": the
  population's change in that condition is 0.020 or more and so is the change of at least 6 of the 8 wirings (each
  wiring read alone) in the same condition; "the other is not lower by more than 0.005" is the population's other
  condition. On fewer wirings (R4 in part A: 4) the same share, rounded up (3 of 4; E108's `wirings_needed`), and the
  verdict is marked provisional until part B.
- "E109's gain holds like for like": in each condition, E109's carrying gate (the six at 0.6 against the lab's input, on
  T against the unseen halves) re-read with the two ρ's under the reading and every other measure as written; its mean
  gain under the reading at least half of the same gain on raw counts (+0.046 / +0.038 as written), and each control
  (permuted, inverted, the control search) below half the re-read candidate's gain; "holds" when both conditions do.
- The gates re-read: the class "rho" re-reads only the ρ's (rho_*, rho15_*, the paa / B-myr distance) under the
  reading, every other measure as written (the APL-effect measures under the calibrated readout, the raw-count loss, the
  bands); the class "all" re-reads every measure under the reading (the entry's "constraints reported under each
  reading"), beside.
- R3 is not available for the archived jittered trials, so R4 stays on raw counts (the entry's definition) and R4z is
  given beside.
- Part B's E109 inputs: the two control searches (ctlA, ctlB) and the first draw of each half's search (A1, B1), at 0.6
  of their size, as their held readings (runs_R.csv's absolute parameters); each model's noise calibrated by the
  pipeline as E72's part E calibrated its four (E43's defaults otherwise, as part E: below-levels tuned, pseudo-recordings
  only at the two matched levels). E43 draws each model's noise by its position in the run, so on a wiring the inputs'
  draws are not the fit's (the changes under R4 are not paired draw for draw, as R3's are).

- The outcome letter: (a) / (b) when the reading moves ρ and E109's gain holds / shrinks under R3 (R4 of E109's inputs
  is part B's, given beside); (c) when no reading moves ρ and no reading's wirings split in sign; (d) when no reading
  moves ρ and R3's or R4's wirings split in sign (at least two of them up and two down in one condition). Provisional
  until part B (`outcome_provisional`).
- E97 saves E92's block in its 'cal' rows and part (R)'s pending block in its 'raw' rows under one column name; the rows
  rebuilt here carry the pending columns under every reading (as part (R) reads them), so the identity against E97's
  'cal' rows leaves those columns out.

Departures from the entry: none so far.

Concentrations: paa's input is the 10^-4 recording for the flies' 10^-3, euc and -aPine are at 10^-1.75 and 10^-1.5 for
the flies' 10^-2; every 120-pair and 54-pair number is affected (`H.APPROX_CONC_NOTE` in the provenance). The holds
(CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the family's hemibrain wiring and its re-draws, no other connectome
table, no raw-only recording; nothing is published.

Outputs: <shared>/results/E112_rho_like_for_like/ (rhos.csv, changes.csv, bound.csv, realism.csv, yardstick.csv,
gates.csv, gates.json, verdict.json, identities.json, provenance_<part>.json); the per-run distance cache under
outputs/E112_rho_like_for_like/.
Usage (from the repository root):
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E112_rho_like_for_like.py --part A --jobs 16
    ./venv/bin/python -u experiments/E112_rho_like_for_like.py --part B-prepare
    (E43's part B per held-out wiring: the commands B-prepare prints)
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E112_rho_like_for_like.py --part B-score --jobs 16
"""
from __future__ import annotations

import os

if __name__ == '__main__':
    for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(_v, '1')

import argparse
import hashlib
import json
import pickle
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E112_rho_like_for_like'
SCRIPT = 'experiments/E112_rho_like_for_like.py'
CONDITIONS = ('on', 'off')
READINGS = ('R0', 'R1', 'R2', 'R3', 'R4')
COUNT_READINGS = ('R0', 'R1', 'R2', 'R3')      # the readings of a run's saved spike counts
REREAD = ('R1', 'R2', 'R3')                    # the readings the saved gates are re-read under
READING_LABEL = {'R0': 'raw counts', 'R1': 'the calibrated saturating readout',
                 'R2': 'the per-KC session z-score (saturating map)',
                 'R3': 'the per-KC session z-score through measurement noise on pseudo-recordings',
                 'R4': 'the threshold-jitter trial noise (raw counts)',
                 'R4z': 'the per-KC session z-score of the threshold-jitter trials'}
N_PSEUDO = {'on': 11, 'off': 8}                # the flies' recordings: 11 control, 8 APL output blocked (TNT)
N_TRIALS = 3
N_R3_DRAWS = 20                                # draws of R3 for its own SD (the yardstick's), beside draw 0
SUBSETS = ('all120', 'megamat66', 'v2odors54', 'no_paa_bmyr91')
PAA_BMYR = ('paa', 'B-myr')
PAA_BMYR_2BUT = ('paa', 'B-myr', '2-but')
MOVE = 0.020                                   # "the reading moves ρ": a population ρ this much higher than under R0
MOVE_WIRINGS = 6                               # ... on 6 of the 8 wirings
OTHER_SLACK = 0.005                            # the other ρ not lower by more than this
E109_HALF = 0.5                                # E109's gain holds: at least half of its raw-count gain
CONTROL_RATIO = 0.5
PART_E = 'E72_wiring_population'               # the folder of part E's pipeline tables
E109_INPUTS_B = ('ctlA', 'ctlB', 'A1', 'B1')   # part B: E109's control searches and one searched input per half
E109_SIZE = 0.6
GATE_EXPERIMENTS = ('E109_input_search_rho', 'E108_paa_extrapolation', 'E97_glomerulus_claw_shares')

_CACHE: dict = {}


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


def _cs(s: pd.Series) -> pd.Series:
    from experiments import E83_cross_wiring_fit as E83
    return E83._cs(s)


def _sha(df: pd.DataFrame) -> str:
    h = hashlib.sha256(np.ascontiguousarray(df.to_numpy(dtype=float), dtype='<f8').tobytes())
    h.update(repr(list(df.index)).encode())
    h.update(repr(list(df.columns)).encode())
    return h.hexdigest()


def wiring_seed(wiring: str) -> int:
    """R3's seed for a wiring: its position among the eight (and the ceiling family's four after them), so that every run
    on one wiring gets the same draws and a change on a wiring is paired."""
    from experiments import E83_cross_wiring_fit as E83
    order = list(E83.ALL_WIRINGS) + ['fc_real', 'fc_P1', 'fc_P2', 'fc_M1']
    return order.index(wiring) if wiring in order else int(hashlib.sha256(wiring.encode()).hexdigest()[:8], 16)


# ---------------------------------------------------------------------------
# the readings

def _counts(o: Mapping[str, tuple], c: str) -> np.ndarray:
    return H.abbrev_columns(o[c][1].astype(float))[list(H.APL16_ODORS)].to_numpy()


def zcell_calibration(on: np.ndarray) -> dict:
    """E73's calibration of the per-KC session z-score (saturating map) on APL-on counts, as E74's `_calibrate_zcell`
    (the control flies' response rate under 1 z, the responders' mean amplitude, the silent KCs' trial sd; q the share of
    the session inside the response windows, the slow drift), for (units x odors) counts or (recordings, trials, units,
    odors) jittered trials."""
    from experiments import E73_per_kc_normalization as E73
    t = E73.refs()['targets']
    q = E73.q_values()[0]
    drift = E73.refs()['drift']
    cal = H.calibrate_session_zscore(on, q=q, amp='saturating', target_rate=t['rate'], target_amp=t['amp'],
                                     target_silent_sd=t['silent_sd'], drift=drift)
    return dict(K=float(cal['K']), noise=float(cal['noise']), noise_ratio=float(cal['noise_ratio']), q=float(q),
                drift=float(drift), converged=bool(cal['converged']), K_at_bound=bool(cal['K_at_bound']),
                noise_at_bound=bool(cal['noise_at_bound']))


def zcell_readings(on: np.ndarray, off: np.ndarray, cal: Optional[dict] = None) -> Tuple[dict, dict, dict]:
    """(z_on, z_off, calibration): both conditions read through `session_zscore` at the calibration set on APL on."""
    cal = zcell_calibration(on) if cal is None else cal
    kw = dict(q=cal['q'], noise=cal['noise'], amp='saturating', K=cal['K'], noise_ratio=cal['noise_ratio'],
              drift=cal['drift'])
    return H.session_zscore(on, **kw), H.session_zscore(off, **kw), cal


def pseudo_recording_distances(z: Mapping[str, object], n_recordings: int, n_trials: int = N_TRIALS,
                               seed: int = 0) -> pd.Series:
    """The mean over `n_recordings` pseudo-recordings of the cosine distances of their trial means (canonical pair
    Series): a `session_zscore` reading's noise-free values (its first recording, units x odors, the APL16 odors in
    order) plus independent Gaussian measurement noise of each KC's trial sd, per trial; seeded."""
    C = np.asarray(z['center'], dtype=float)[0]
    sd = np.asarray(z['trial_sd'], dtype=float)[0]
    rng = np.random.default_rng(seed)
    ds = []
    for _ in range(int(n_recordings)):
        tr = C[None] + rng.standard_normal((int(n_trials),) + C.shape) * sd[None, :, None]
        m = pd.DataFrame(tr.mean(axis=0), columns=list(H.APL16_ODORS))
        ds.append(_cs(H.pair_distances(m, 'cosine')))
    return _cs(pd.concat(ds, axis=1).mean(axis=1))


def r3_from_z(z_on: dict, z_off: dict, seed: int = 0) -> Dict[str, pd.Series]:
    """R3 of a run read as `zcell_readings`: APL on through `N_PSEUDO['on']` pseudo-recordings, APL off through
    `N_PSEUDO['off']`, independent draws per condition (seeds from `seed`)."""
    ss = np.random.SeedSequence(int(seed)).spawn(2)
    return {c: pseudo_recording_distances(z, N_PSEUDO[c], N_TRIALS, int(s.generate_state(1)[0]))
            for (c, z), s in zip((('on', z_on), ('off', z_off)), ss)}


def reading_distances(o: Mapping[str, tuple], reading: str, seed: int = 0) -> Dict[str, pd.Series]:
    """{'on', 'off'} pair distances (canonical order) of a saved run `o` ({'on' / 'off': (responses, spike counts)})
    under reading R0 (raw counts), R1 (E65's calibrated readout, K on the run's APL-on counts), R2 (E73's per-KC session
    z-score, saturating map, graded, in expectation over the measurement noise) or R3 (R2 through sampled measurement
    noise on pseudo-recordings: 11 with APL on and 8 with APL off, 3 trials each, the distance averaged over them)."""
    return run_readings(o, (reading,), seed)[reading]


def run_readings(o: Mapping[str, tuple], readings: Sequence[str] = COUNT_READINGS, seed: int = 0,
                 r3_draws: Sequence[int] = ()) -> Dict[str, Dict[str, pd.Series]]:
    """`reading_distances` of one run under several readings at once (R2 and R3 share one calibration); with
    `r3_draws`, also R3 at seeds seed + 1000 x k for each k (keys 'R3_d<k>')."""
    from experiments import E83_cross_wiring_fit as E83
    out = {}
    for r in readings:
        if r == 'R0':
            out[r] = E83.pair_distances(o, 'raw')
        elif r == 'R1':
            out[r] = E83.pair_distances(o, 'cal')
        elif r not in ('R2', 'R3'):
            raise ValueError(f'{r!r}: a run\'s counts are read under {COUNT_READINGS}')
    if {'R2', 'R3'} & set(readings) or r3_draws:
        z_on, z_off, _ = zcell_readings(_counts(o, 'on'), _counts(o, 'off'))
        if 'R2' in readings:
            out['R2'] = {'on': _cs(H.zcell_pair_distances(z_on, 'graded')), 'off': _cs(H.zcell_pair_distances(z_off, 'graded'))}
        if 'R3' in readings:
            out['R3'] = r3_from_z(z_on, z_off, seed)
        for k in r3_draws:
            out[f'R3_d{k}'] = r3_from_z(z_on, z_off, seed + 1000 * int(k))
    return out


def population(ds: Sequence[Mapping[str, pd.Series]]) -> Dict[str, pd.Series]:
    """The mean over runs of their pair distances per condition (E72's population; canonical order)."""
    return {c: _cs(pd.concat([_cs(d[c]) for d in ds], axis=1).mean(axis=1)) for c in CONDITIONS}


# --- R4: the pipeline's trial noise

def _model_wiring(model: str) -> Optional[str]:
    """The family's wiring label of an E43 pipeline model name (E<n>-<tag>), or None."""
    from experiments import E83_cross_wiring_fit as E83
    for w in E83.ALL_WIRINGS:
        if model.split('-', 1)[-1] == E83.saved_fit(w)['tag']:
            return w
    return None


def _levels(pairs_long: pd.DataFrame) -> pd.DataFrame:
    """The pipeline table with a level 'match_on+match_off' (a model whose two matched sigmas coincide) counted as both."""
    both = pairs_long[pairs_long['level'] == 'match_on+match_off']
    return pd.concat([pairs_long[pairs_long['level'] != 'match_on+match_off'], both.assign(level='match_on'),
                      both.assign(level='match_off')], ignore_index=True)


def trial_noise_distances(pairs_long: pd.DataFrame, column: str = 'raw') -> Dict[str, Dict[str, pd.Series]]:
    """R4 from E43's pipeline table (partB/pairs_noise_long.csv): {model: {'on', 'off'}} the mean over the
    pseudo-recordings of the distances of their three-trial means (`column` 'raw', graded spike counts; 'det' gives the
    deterministic run), APL on at the level matched to the control flies (match_on) and APL off at the level matched to
    APL output blocked (TNT) (match_off)."""
    d = _levels(pairs_long)
    d = d[d['transform'] == 'graded']
    out = {}
    for m, g in d.groupby('model', sort=False):
        r = {}
        for c, cond, level in (('on', 'on', 'match_on'), ('off', 'acute_off', 'match_off')):
            x = g[(g['condition'] == cond) & (g['level'] == level)]
            if len(x) != 120:
                break
            r[c] = _cs(pd.Series(x[column].to_numpy(dtype=float),
                                 index=pd.MultiIndex.from_arrays([x['odor_a'], x['odor_b']], names=['odor_a', 'odor_b'])))
        if len(r) == 2:
            out[m] = r
    return out


def trial_sigmas(res_dir: Path) -> Dict[str, Dict[str, float]]:
    """{model: {'on': match_on sigma, 'off': match_off sigma}} of a pipeline results folder (noise_levels.json via
    provenance.json's 'levels')."""
    prov = json.loads((res_dir / 'provenance.json').read_text())
    out = {}
    for m, lv in prov['levels'].items():
        s = lv.get('match_on+match_off')
        out[m] = {'on': lv.get('match_on', s), 'off': lv.get('match_off', s)}
    return out


def trial_zcell_distances(trials_dir: Path, model: str, sigma_on: float, sigma_off: float,
                          n_recordings: int = 10) -> Optional[Dict[str, pd.Series]]:
    """R4z: the per-KC session z-score (saturating map, graded) of a pipeline model's saved jittered trials, APL on from
    the trials at `sigma_on` and APL off from those at `sigma_off`, calibrated on the APL-on trials (E73's `_trials`
    reading); None when the trials are not saved."""
    def load(sig, key):
        fs = [trials_dir / f'{model}_sigma{sig:g}_n{N_TRIALS}_rep{r:02d}.npz' for r in range(n_recordings)]
        if not all(f.exists() for f in fs):
            return None
        return np.stack([np.load(f)[key].astype(float) for f in fs])
    t_on, t_off = load(sigma_on, 'on_spk'), load(sigma_off, 'acute_off_spk')
    if t_on is None or t_off is None:
        return None
    z_on, z_off, _ = zcell_readings(t_on, t_off)
    return {'on': _cs(H.zcell_pair_distances(z_on, 'graded')), 'off': _cs(H.zcell_pair_distances(z_off, 'graded'))}


# ---------------------------------------------------------------------------
# the scores

def data_matrices() -> Dict[str, pd.Series]:
    from experiments import E106_rho_gap_parts as E106
    return {c: _cs(v) for c, v in E106.data_matrices().items()}


def _pair_value(s: pd.Series, a: str, b: str) -> float:
    s = _cs(s)
    for k in ((a, b), (b, a)):
        if k in s.index:
            return float(s.loc[k])
    return np.nan


def reading_scores(d: Mapping[str, pd.Series], data: Optional[Mapping[str, pd.Series]] = None,
                   effect: Optional[pd.Series] = None) -> dict:
    """The scores of one unit's distances under a reading: the two ρ's over the 120, 66, 54 and 91 pairs
    ('rho_<on|off>_<subset>'), the distance of paa from B-myr ('paa_bmyr_<on|off>'), the APL-effect loss with its mean,
    pattern, amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair APL effect, the mean APL effect
    and the mean distances. `data` defaults to the flies' mean matrices, `effect` to the flies' per-pair APL effect."""
    from experiments import E106_rho_gap_parts as E106
    D = data_matrices() if data is None else data
    eff = H.yang_apl_effect() if effect is None else effect
    on, off = _cs(d['on']), _cs(d['off'])
    row = {}
    for c, m in (('on', on), ('off', off)):
        for s in SUBSETS:
            row[f'rho_{c}_{s}'] = float(E106.rho(m, D[c], s))
        row[f'paa_bmyr_{c}'] = _pair_value(m, *PAA_BMYR)
        row[f'mean_{c}'] = float(m.mean())
    st, _ = H.apl_effect_decomposition(on, off, eff)
    sp = H.apl_effect_spread_shape(on, off, eff)
    row.update(loss=float(st['apl_effect_rmse']), mean_effect=float(st['apl_effect_mean']),
               pattern=float(st['apl_effect_rmse_centered']), rho_effect=float(st['apl_effect_spearman']),
               pearson_effect=float(st['apl_effect_pearson']), amplitude=float(sp['apl_effect_amplitude']),
               shape=float(sp['apl_effect_shape']), sd_ratio=float(sp['apl_effect_sd_ratio']))
    return row


HYBRID_KEYS = ('model', 'flies_odor_model_pair', 'model_odor_flies_pair', 'flies_odor_model_pair_quantile',
               'model_odor_flies_pair_quantile', 'model_pair_optimal_odor', 'odor_shares_rho', 'pair_parts_rho')


def bound_scores(d: Mapping[str, pd.Series], data: Optional[Mapping[str, pd.Series]] = None) -> dict:
    """E106's bound under a reading: per condition the two ρ's with the pairs of paa and B-myr (29) and of paa, B-myr and
    2-but (42) at the flies' values ('at_<set>_<raw|quantile>_<on|off>', E106's `at_data_values`), and the odor part and
    pair part hybrids (E106's `odor_pair_hybrids`: '<key>_<on|off>' for HYBRID_KEYS)."""
    from experiments import E106_rho_gap_parts as E106
    D = data_matrices() if data is None else data
    row = {}
    for c in CONDITIONS:
        m = _cs(d[c])
        for name, ods in (('paa_bmyr', PAA_BMYR), ('paa_bmyr_2but', PAA_BMYR_2BUT)):
            for how in ('raw', 'quantile'):
                row[f'at_{name}_{how}_{c}'] = H.spearman(E106.at_data_values(m, D[c], ods, how), D[c])
        hy = E106.odor_pair_hybrids(m, D[c])
        row.update({f'{k}_{c}': float(hy[k]) for k in HYBRID_KEYS})
    return row


def wirings_needed(n: int) -> int:
    """MOVE_WIRINGS of 8, or the same share of n rounded up."""
    return int(np.ceil(MOVE_WIRINGS / 8 * n))


def moves_verdict(rows: pd.DataFrame, reading: str, base: str = 'R0', pop: Optional[str] = None) -> dict:
    """The entry's "the reading moves ρ" for one `reading` against `base` from rhos.csv's layout (columns reading, unit,
    kind ('wiring' / 'population'), rho_<on|off>_all120): the population `pop` (default the eight where present, else
    the largest population the reading has) and the wirings in it. Per condition: the population's change, the changes of
    its wirings, the count at MOVE or more, `moves_<c>` (the population's change at MOVE or more, as many wirings as
    `wirings_needed`, the other condition's population change not below -OTHER_SLACK); `moves` (either condition),
    `provisional` (fewer than 8 wirings), the sign split (wirings up and down) and the paired SD over the wirings."""
    from experiments import E83_cross_wiring_fit as E83
    pops = {'all8': E83.ALL_WIRINGS, 'S': E83.S_WIRINGS, 'T': E83.T_WIRINGS}
    have = set(rows.loc[(rows['reading'] == reading) & (rows['kind'] == 'wiring'), 'unit'])
    if pop is None:
        pop = next((p for p in ('all8', 'S', 'T') if set(pops[p]) <= have), None)
    out = dict(reading=reading, base=base, population=pop)
    if pop is None:
        out.update(moves=None, note='no population')
        return out
    ws = list(pops[pop])
    v = lambda r, u, c: float(rows.loc[(rows['reading'] == r) & (rows['unit'] == u), f'rho_{c}_all120'].iloc[0])
    other = {'on': 'off', 'off': 'on'}
    dpop = {c: v(reading, pop, c) - v(base, pop, c) for c in CONDITIONS}
    for c in CONDITIONS:
        dw = {w: v(reading, w, c) - v(base, w, c) for w in ws}
        n_up = int(sum(x >= MOVE for x in dw.values()))
        out.update({f'd_pop_{c}': dpop[c], f'd_wirings_{c}': dw, f'n_wirings_move_{c}': n_up,
                    f'n_up_{c}': int(sum(x > 0 for x in dw.values())), f'n_down_{c}': int(sum(x < 0 for x in dw.values())),
                    f'paired_sd_{c}': float(np.std(list(dw.values()), ddof=1)) if len(dw) > 1 else np.nan,
                    f'moves_{c}': bool(dpop[c] >= MOVE and n_up >= wirings_needed(len(ws)) and dpop[other[c]] >= -OTHER_SLACK)})
    out['n_wirings'] = len(ws)
    out['wirings_needed'] = wirings_needed(len(ws))
    out['moves'] = bool(out['moves_on'] or out['moves_off'])
    out['provisional'] = len(ws) < 8
    out['sign_split'] = bool(any(min(out[f'n_up_{c}'], out[f'n_down_{c}']) >= 2 for c in CONDITIONS))
    return out


# ---------------------------------------------------------------------------
# per-run distances, cached and in parallel

def _cache() -> dict:
    if 'store' not in _CACHE:
        p = local_dir() / 'distances_cache.pkl'
        _CACHE['store'] = pickle.loads(p.read_bytes()) if p.exists() else {}
    return _CACHE['store']


def _flush_cache() -> None:
    if _CACHE.get('dirty'):
        (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['store']))
        _CACHE['dirty'] = False


_LOADERS: Dict[str, Callable[[str, str], Mapping[str, tuple]]] = {}


def loader(experiment: str) -> Callable[[str, str], Mapping[str, tuple]]:
    """The saved run of (variant, wiring) of a gate experiment ('family' for the family's saved fits)."""
    if experiment not in _LOADERS:
        if experiment == 'family':
            from experiments import E83_cross_wiring_fit as E83
            _LOADERS[experiment] = lambda v, w: E83.saved_fit(w)['saved']
        elif experiment == 'E109_input_search_rho':
            from experiments import E109_input_search_rho as X
            _LOADERS[experiment] = X.run_outputs
        elif experiment == 'E108_paa_extrapolation':
            from experiments import E108_paa_extrapolation as X
            _LOADERS[experiment] = X.run_outputs
        elif experiment == 'E97_glomerulus_claw_shares':
            from experiments import E97_glomerulus_claw_shares as X
            _LOADERS[experiment] = X.run_outputs
        else:
            raise ValueError(experiment)
    return _LOADERS[experiment]


def _job(job: tuple) -> Tuple[tuple, dict, Optional[str]]:
    """One run's readings (worker): job (experiment, variant, wiring, readings, seed, r3_draws)."""
    experiment, variant, wiring, readings, seed, draws = job
    warnings.simplefilter('ignore', RuntimeWarning)
    try:
        o = loader(experiment)(variant, wiring)
        key = (_sha(o['on'][1]), _sha(o['off'][1]), int(seed))
        return key, run_readings(o, readings, seed, draws), None
    except Exception as e:
        import traceback
        return (experiment, variant, wiring), {}, f'{type(e).__name__}: {e}\n{traceback.format_exc()}'


def run_distances(runs: Sequence[Tuple[str, str, str]], readings: Sequence[str] = COUNT_READINGS, jobs: int = 1,
                  r3_draws: Sequence[int] = ()) -> Dict[Tuple[str, str, str], Dict[str, Dict[str, pd.Series]]]:
    """{(experiment, variant, wiring): {reading: {'on', 'off'}}} for every run, R3 seeded by `wiring_seed`; cached by
    the sha256 of the run's spike counts and the seed (outputs/.../distances_cache.pkl), computed in `jobs` processes."""
    store = _cache()
    want = [f'R3_d{k}' for k in r3_draws]
    out, todo = {}, []
    keys = {}
    for ex, v, w in runs:
        o = loader(ex)(v, w)
        k = (_sha(o['on'][1]), _sha(o['off'][1]), wiring_seed(w))
        keys[(ex, v, w)] = k
        have = store.get(k, {})
        need = [r for r in readings if r not in have]
        need_d = [int(x.split('_d')[1]) for x in want if x not in have]
        if need or need_d:
            todo.append((ex, v, w, tuple(need), wiring_seed(w), tuple(need_d)))
    if todo:
        print(f'E112: {len(todo)} runs to read ({jobs} processes)', flush=True)
        if jobs > 1:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=jobs, mp_context=mp.get_context('fork')) as ex_:
                res = list(ex_.map(_job, todo, chunksize=1))
        else:
            res = [_job(j) for j in todo]
        errs = [e for _, _, e in res if e]
        if errs:
            raise RuntimeError(f'{len(errs)} runs failed; the first:\n{errs[0]}')
        for k, d, _ in res:
            store.setdefault(k, {}).update(d)
        _CACHE['dirty'] = True
        _flush_cache()
    for r, k in keys.items():
        out[r] = {x: store[k][x] for x in list(readings) + want}
    return out


# ---------------------------------------------------------------------------
# part A: the family

def family_runs() -> List[Tuple[str, str, str]]:
    from experiments import E83_cross_wiring_fit as E83
    return [('family', 'fit', w) for w in E83.ALL_WIRINGS]


def part_e_r4() -> Tuple[Dict[str, Dict[str, pd.Series]], Dict[str, Dict[str, pd.Series]], Dict[str, Dict[str, pd.Series]]]:
    """(R4 by wiring, the pipeline's deterministic distances by wiring, R4z by wiring) of E72's part E (the four wirings
    of S), plus part B's wirings once their pipeline tables exist (results/E112_rho_like_for_like/partB_<wiring>/)."""
    dirs = [(H.results_dir(PART_E) / 'partB', H.shared_dir() / 'data' / 'model_outputs' / PART_E / 'partB' / 'trials')]
    dirs += [(p, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / p.name / 'trials')
             for p in sorted(out_dir().glob('partB_*')) if (p / 'pairs_noise_long.csv').exists()]
    r4, det, r4z = {}, {}, {}
    for res, trials in dirs:
        pl = pd.read_csv(res / 'pairs_noise_long.csv')
        R, Dt = trial_noise_distances(pl), trial_noise_distances(pl, 'det')
        sig = trial_sigmas(res)
        if not trials.exists():
            trials = H.outputs_dir(EXPERIMENT) / res.name / 'trials'
        for m, d in R.items():
            w = _model_wiring(m)
            if w is None:
                continue
            r4[w], det[w] = d, Dt[m]
            z = trial_zcell_distances(trials, m, sig[m]['on'], sig[m]['off'])
            if z is not None:
                r4z[w] = z
    return r4, det, r4z


def family_readings(jobs: int = 1) -> Dict[str, Dict[str, Dict[str, pd.Series]]]:
    """{reading: {wiring: {'on', 'off'}}}: R0 to R3 of the eight saved fits, R3's further draws ('R3_d<k>'), R4 and R4z
    where the pipeline ran."""
    rd = run_distances(family_runs(), COUNT_READINGS, jobs, r3_draws=range(1, N_R3_DRAWS))
    out: Dict[str, Dict[str, Dict[str, pd.Series]]] = {}
    for (_, _, w), d in rd.items():
        for r, x in d.items():
            out.setdefault(r, {})[w] = x
    r4, det, r4z = part_e_r4()
    out['R4'], out['R4det'], out['R4z'] = r4, det, r4z
    return out


def unit_distances(per: Mapping[str, Dict[str, pd.Series]]) -> Dict[str, Dict[str, pd.Series]]:
    """Each wiring alone and every population (the eight, S, T) whose wirings all have the reading."""
    from experiments import E83_cross_wiring_fit as E83
    out = dict(per)
    for p, ws in (('all8', E83.ALL_WIRINGS), ('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS)):
        if all(w in per for w in ws):
            out[p] = population([per[w] for w in ws])
    return out


def rhos_table(fam: Mapping[str, Mapping[str, Dict[str, pd.Series]]]) -> pd.DataFrame:
    """rhos.csv: per reading (R0 to R4, R4z) and unit (each wiring, the populations) `reading_scores`; R3 also as the
    mean and SD over its draws of the two ρ's ('r3_draws_*')."""
    rows = []
    for r in ('R0', 'R1', 'R2', 'R3', 'R4', 'R4z'):
        for unit, d in unit_distances(fam.get(r, {})).items():
            row = dict(reading=r, unit=unit, kind='population' if unit in ('all8', 'S', 'T') else 'wiring',
                       **reading_scores(d))
            if r == 'R3':
                draws = [k for k in fam if k.startswith('R3_d')]
                vals = {c: [] for c in CONDITIONS}
                for k in ['R3'] + draws:
                    du = unit_distances(fam[k]).get(unit)
                    if du is not None:
                        for c in CONDITIONS:
                            vals[c].append(H.spearman(du[c], data_matrices()[c]))
                for c in CONDITIONS:
                    row[f'r3_draws_rho_{c}_mean'] = float(np.mean(vals[c]))
                    row[f'r3_draws_rho_{c}_sd'] = float(np.std(vals[c], ddof=1)) if len(vals[c]) > 1 else np.nan
                    row['r3_n_draws'] = len(vals[c])
            rows.append(row)
    return pd.DataFrame(rows)


def bound_table(fam: Mapping[str, Mapping[str, Dict[str, pd.Series]]]) -> pd.DataFrame:
    """bound.csv: `bound_scores` of each reading's populations."""
    rows = []
    for r in ('R0', 'R1', 'R2', 'R3', 'R4', 'R4z'):
        for unit, d in unit_distances(fam.get(r, {})).items():
            if unit in ('all8', 'S', 'T'):
                rows.append(dict(reading=r, unit=unit, **bound_scores(d)))
    return pd.DataFrame(rows)


def changes_table(rt: pd.DataFrame) -> pd.DataFrame:
    """changes.csv: each wiring's and population's two ρ's (120 pairs) under each reading minus under R0."""
    base = rt[rt['reading'] == 'R0'].set_index('unit')
    rows = []
    for r in rt['reading'].unique():
        if r == 'R0':
            continue
        for x in rt[rt['reading'] == r].itertuples():
            rows.append(dict(reading=r, unit=x.unit, kind=x.kind, **{f'd_rho_{c}': getattr(x, f'rho_{c}_all120') -
                                                                  base.loc[x.unit, f'rho_{c}_all120'] for c in CONDITIONS}))
    return pd.DataFrame(rows)


def realism_table(jobs: int = 1) -> pd.DataFrame:
    """realism.csv: per wiring the response rates, their ratio, response rate ρ and the odors-per-KC distances from the
    flies' with APL on and off: on two spikes (the bands' rule, E83's `realism`) and under the per-KC session z-score's
    criterion (above 1 z, in expectation over the noise; E73's `criterion_stats`), beside the flies' values."""
    from experiments import E73_per_kc_normalization as E73
    from experiments import E83_cross_wiring_fit as E83
    rows = []
    for w in E83.ALL_WIRINGS:
        o = E83.saved_fit(w)['saved']
        rows.append(dict(wiring=w, reading='two spikes (R0, R1; the bands)', **E83.realism(o)))
        z_on, z_off, cal = zcell_readings(_counts(o, 'on'), _counts(o, 'off'))
        cs = E73.criterion_stats(H.zcell_stats(z_on), H.zcell_stats(z_off), breadth_checks=False)
        rows.append(dict(wiring=w, reading='above 1 z (R2, R3)', **{k: v for k, v in cs.items() if np.isscalar(v)},
                         **{f'cal_{k}': v for k, v in cal.items()}))
    dt = H.data_target
    rows.append(dict(wiring='flies', reading='KC data (baseline-noise units, 1.97 SD)', rate_on=dt('ctrl_frac'),
                     rate_off=dt('tnt_frac'), ratio=dt('tnt_frac') / dt('ctrl_frac')))
    return pd.DataFrame(rows)


def yardstick_table(fam: Mapping[str, Mapping[str, Dict[str, pd.Series]]]) -> pd.DataFrame:
    """yardstick.csv: per reading, S's population against T's (what a model that was the flies would score against itself
    through the reading: ρ between the two populations' matrices, per condition), and, where a reading has only S (R4 in
    part A), its two halves (real + P1 against P2 + M1) beside; the family's population against the flies' noise-free
    distances (E81); the data's split-half ρ's and ceilings (registry) and the mean matrices against the noise-free ones."""
    from experiments import E106_rho_gap_parts as E106
    nf = {c: _cs(v) for c, v in E106.noise_free_data().items()}
    D = data_matrices()
    rows = []
    for r in ('R0', 'R1', 'R2', 'R3', 'R4', 'R4z'):
        U = unit_distances(fam.get(r, {}))
        if 'S' in U and 'T' in U:
            rows.append(dict(reading=r, what='S population against T population', **{
                f'rho_{c}': H.spearman(U['S'][c], U['T'][c]) for c in CONDITIONS}))
        per = fam.get(r, {})
        if all(w in per for w in ('real', 'P1', 'P2', 'M1')):
            a, b = population([per['real'], per['P1']]), population([per['P2'], per['M1']])
            rows.append(dict(reading=r, what='real + P1 against P2 + M1 (two wirings each)',
                             **{f'rho_{c}': H.spearman(a[c], b[c]) for c in CONDITIONS}))
        pop = U.get('all8', U.get('S'))
        if pop is not None:
            rows.append(dict(reading=r, what=f"the family ({'eight' if 'all8' in U else 'S, four'}) against the flies' "
                             'noise-free distances (E81)', **{f'rho_{c}': H.spearman(pop[c], nf[c]) for c in CONDITIONS}))
            rows.append(dict(reading=r, what=f"the family ({'eight' if 'all8' in U else 'S, four'}) against the flies' "
                             'mean distances', **{f'rho_{c}': H.spearman(pop[c], D[c]) for c in CONDITIONS}))
    dt = H.data_target
    rows.append(dict(reading='data', what='split-half ρ of the recordings (registry ctrl_ / tnt_split_half_rho)',
                     rho_on=dt('ctrl_split_half_rho'), rho_off=dt('tnt_split_half_rho')))
    rows.append(dict(reading='data', what='measurement ceiling (registry ctrl_ / tnt_ceiling)',
                     rho_on=dt('ctrl_ceiling'), rho_off=dt('tnt_ceiling')))
    rows.append(dict(reading='data', what="the flies' mean distances against their noise-free distances (E81)",
                     **{f'rho_{c}': H.spearman(D[c], nf[c]) for c in CONDITIONS}))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the saved gates read again

def _rho_cols(cols: Sequence[str]) -> List[str]:
    return [c for c in cols if c.startswith(('rho_', 'rho15_', 'paa_bmyr_'))]


def gate_units(experiment: str) -> List[Tuple[str, str, List[Tuple[str, str]]]]:
    """[(variant, unit, [(variant, wiring), ...])] of every reading a gate experiment's readings.csv holds."""
    if experiment == 'E109_input_search_rho':
        from experiments import E109_input_search_rho as X
        runs = X.runs_table()
        out = []
        for ch, g in runs.groupby('change', sort=False):
            if set(g['wiring']) != set(X.T_WIRINGS):
                continue
            out += [(ch, w, [(ch, w)]) for w in X.T_WIRINGS] + [(ch, 'T', [(ch, w) for w in X.T_WIRINGS])]
        return out
    if experiment == 'E108_paa_extrapolation':
        from experiments import E108_paa_extrapolation as X
        runs = X.runs_table()
        out = []
        for name, g in runs.groupby('input', sort=False):
            have = list(g['unit'])
            out += [(name, u, [(name, u)]) for u in have]
            out += [(name, u, [(name, x) for x in us]) for u, us in X.UNITS if all(x in have for x in us)]
        return out
    if experiment == 'E97_glomerulus_claw_shares':
        from experiments import E97_glomerulus_claw_shares as X
        return [(v, u, list(mem)) for v, u, mem in X.populations(X.runs_table())]
    raise ValueError(experiment)


def gate_scores(experiment: str, d: Mapping[str, pd.Series], n_members: int) -> dict:
    """The experiment's own scores of one unit's distances (E108's `unit_scores` for E108 and E109, with E109's distance
    of paa from B-myr; E97's `score_row` and `rho_scores`, the pending columns on populations)."""
    if experiment in ('E109_input_search_rho', 'E108_paa_extrapolation'):
        from experiments import E108_paa_extrapolation as E108
        row = E108.unit_scores(d)
        if experiment == 'E109_input_search_rho':
            row.update(paa_bmyr_on=_pair_value(d['on'], *PAA_BMYR), paa_bmyr_off=_pair_value(d['off'], *PAA_BMYR))
        return row
    from experiments import E97_glomerulus_claw_shares as X
    return dict(X.score_row(d), **X.rho_scores(d, pending=n_members > 1))


def _key_cols(experiment: str) -> Tuple[str, str]:
    return {'E109_input_search_rho': ('change', 'unit'), 'E108_paa_extrapolation': ('input', 'unit'),
            'E97_glomerulus_claw_shares': ('variant', 'unit')}[experiment]


def gate_readings(experiment: str, readings: Sequence[str], jobs: int = 1) -> Dict[str, pd.DataFrame]:
    """{reading: rows} of the experiment's readings rebuilt from its saved runs under each reading (the variant / unit
    columns of its readings.csv, readout = the reading, and the scores of `gate_scores`)."""
    readings = tuple(dict.fromkeys(readings))
    units = gate_units(experiment)
    runs = sorted({(experiment, v, w) for _, _, mem in units for v, w in mem})
    _CACHE['gate_rd'] = run_distances(runs, readings, jobs)
    _CACHE['gate_units'] = units
    kv, ku = _key_cols(experiment)
    todo = [(experiment, r, i) for r in readings for i in range(len(units))]
    if jobs > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs, mp_context=mp.get_context('fork')) as ex_:
            scored = list(ex_.map(_score_job, todo, chunksize=4))
    else:
        scored = [_score_job(j) for j in todo]
    out = {}
    for (_, r, i), sc in zip(todo, scored):
        v, u, _ = units[i]
        out.setdefault(r, []).append({kv: v, ku: u, 'readout': r, **sc})
    return {r: pd.DataFrame(rows) for r, rows in out.items()}


def _score_job(job: tuple) -> dict:
    """`gate_scores` of one unit under one reading (worker; the distances are `gate_readings`' in this process)."""
    experiment, r, i = job
    warnings.simplefilter('ignore', RuntimeWarning)
    _, _, mem = _CACHE['gate_units'][i]
    ds = [_CACHE['gate_rd'][(experiment, mv, mw)][r] for mv, mw in mem]
    d = ds[0] if len(ds) == 1 else population(ds)
    return gate_scores(experiment, d, len(mem))


def saved_readings(experiment: str) -> pd.DataFrame:
    return pd.read_csv(H.results_dir(experiment) / 'readings.csv', low_memory=False)


def substitute(saved: pd.DataFrame, new: pd.DataFrame, experiment: str, which: str) -> pd.DataFrame:
    """The saved readings with the 'raw' rows' ρ columns (`which` 'rho') or the 'raw' and 'cal' rows' every score
    column (`which` 'all') taken from `new` (one reading's rows), matched on the variant and unit; the run means, the
    runs flagged and every column `new` lacks kept as saved."""
    kv, ku = _key_cols(experiment)
    out = saved.copy()
    nn = new.set_index([kv, ku])
    score_cols = [c for c in new.columns if c not in (kv, ku, 'readout')]
    cols = _rho_cols(score_cols) if which == 'rho' else score_cols
    targets = ('raw',) if which == 'rho' else ('raw', 'cal')
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
        m = out['readout'].isin(targets)
        keys = list(zip(out.loc[m, kv], out.loc[m, ku]))
        out.loc[m, c] = [nn[c].get(k, np.nan) for k in keys]
    return out


def compare_rows(saved: pd.DataFrame, new: pd.DataFrame, experiment: str, readout: str,
                 skip: Sequence[str] = ()) -> dict:
    """The largest absolute difference, over the score columns both have (but `skip`), between the saved readings of
    `readout` and rebuilt rows (an identity: 0 up to float noise), and the number of rows matched."""
    kv, ku = _key_cols(experiment)
    s = saved[saved['readout'] == readout].set_index([kv, ku])
    n = new.set_index([kv, ku])
    idx = s.index.intersection(n.index)
    cols = [c for c in n.columns if c in s.columns and c != 'readout' and c not in skip
            and pd.api.types.is_numeric_dtype(s[c])]
    diff = 0.0
    worst = None
    for c in cols:
        a, b = s.loc[idx, c].astype(float).to_numpy(), n.loc[idx, c].astype(float).to_numpy()
        if not np.isfinite(a).any():
            # (a column the experiment saved only for other readouts, e.g. E97's pending columns, raw rows only)
            continue
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.any():
            x = float(np.max(np.abs(a[ok] - b[ok])))
            if x > diff:
                diff, worst = x, c
        if (np.isfinite(a) != np.isfinite(b)).any():
            diff, worst = np.inf, f'{c} (NaN pattern)'
    return dict(max_abs_diff=diff, worst_column=worst, n_rows=len(idx), n_saved=len(s), n_new=len(n))


def gate_classes(experiment: str, rd: pd.DataFrame) -> Dict[str, dict]:
    """The experiment's own gates on readings `rd`: {key: gate dict with its class under 'cls'}. E109: the six at 0.6
    against the lab's input (carrying; 'six:I0:s0.6'), at 1, against E108's input, and the control searches; E108: each
    rule and control on the eight and on the ceiling family's four; E97: part (R)'s class of the mean proportions (and of
    the second draw where its runs are there)."""
    out = {}
    if experiment == 'E109_input_search_rho':
        from experiments import E109_input_search_rho as X
        runs = X.runs_table()
        refs = ['I0'] + (['E108ref'] if (rd['change'] == 'E108ref').any() else [])
        for ref in refs:
            for size in (X.GATE_SIZE, X.FULL_SIZE):
                out[f'six:{ref}:s{size:g}'] = X.gate(rd, runs, ref=ref, size=size)
            out[f'ctl:{ref}:s{X.GATE_SIZE:g}'] = X.gate(rd, runs, ref=ref, main=X.CONTROL_SEARCHES, second=(),
                                                         size=X.GATE_SIZE, with_controls=False)
        return out
    if experiment == 'E108_paa_extrapolation':
        from experiments import E108_paa_extrapolation as X
        runs = X.runs_table()
        for pop, units in (('all8', X.ALL_WIRINGS), ('fc4', X.FC_WIRINGS)):
            for name in X.RULES + X.CONTROLS:
                g = X.gate(rd, runs, name, units=units, pop=pop)
                if g.get('complete'):
                    g['cls'] = X.class_of(g)
                out[f'{pop}:{name}'] = g
        return out
    from experiments import E97_glomerulus_claw_shares as X
    out['B2'] = X.rho_class(rd, f'B2_{X.MEAN}')
    d2 = X.rho_class(rd, X.D2)
    if d2.get('cls') != 'incomplete':
        out['draw2'] = d2
    return out


def _flat_gate(g: Mapping) -> dict:
    keep = {}
    for k, v in g.items():
        if isinstance(v, (bool, int, float, str, np.floating, np.integer)) or v is None:
            keep[k] = v
        elif k in ('fails', 'pending_fails') and isinstance(v, dict):
            keep[k] = '; '.join(f'{a} {b:+.4f}' if isinstance(b, float) else f'{a} {b}' for a, b in v.items())
        elif k == 'controls' and isinstance(v, dict):
            for c, x in v.items():
                if isinstance(x, dict):
                    keep.update({f'ctl_{c}_{a}': b for a, b in x.items() if isinstance(b, (bool, int, float))})
                else:
                    keep[f'ctl_{c}'] = x
    return keep


def reread_gate(experiment: str, readings: Sequence[str] = REREAD, jobs: int = 1) -> Tuple[pd.DataFrame, dict]:
    """A saved gate (E109's held readings, E108's runs, E97's points) re-scored under `readings`: the identities first
    (the rows rebuilt under R0 and R1 against the saved 'raw' and 'cal' rows, and the classes recomputed from the saved
    readings against the saved ones), then per reading the classes with the ρ's re-read ('rho') and with every measure
    re-read ('all'). Returns (one row per gate, reading and kind; the identities and the full gate dicts)."""
    saved = saved_readings(experiment)
    rows_by = gate_readings(experiment, ('R0', 'R1') + tuple(readings), jobs)
    # (E97 saves E92's block in its 'cal' rows and part (R)'s pending block in its 'raw' rows, under one name; the rows
    # rebuilt here carry the pending columns under every reading, as the gate reads them, so they are not compared
    # with the 'cal' rows)
    skip_cal = ()
    if experiment == 'E97_glomerulus_claw_shares':
        from experiments import E97_glomerulus_claw_shares as X
        skip_cal = X.PENDING_COLS
    ident = {f'{r}_vs_saved_{ro}': compare_rows(saved, rows_by[r], experiment, ro, skip)
             for r, ro, skip in (('R0', 'raw', ()), ('R1', 'cal', skip_cal))}
    as_written = gate_classes(experiment, saved)
    rows, full = [], {'identities': ident, 'as_written': as_written}
    for key, g in as_written.items():
        rows.append(dict(experiment=experiment, gate=key, reading='as written', kind='as written', **_flat_gate(g)))
    for r in readings:
        for kind in ('rho', 'all'):
            g_all = gate_classes(experiment, substitute(saved, rows_by[r], experiment, kind))
            full[f'{r}:{kind}'] = g_all
            for key, g in g_all.items():
                rows.append(dict(experiment=experiment, gate=key, reading=r, kind=kind, **_flat_gate(g)))
    for r, df in rows_by.items():
        df.to_csv(out_dir() / f'gate_readings_{experiment}_{r}.csv', index=False)
    return pd.DataFrame(rows), full


def saved_classes(experiment: str) -> Dict[str, str]:
    """The classes as the experiment wrote them (gate.csv / verdict_rho.json), by `gate_classes`' keys."""
    res = H.results_dir(experiment)
    if experiment == 'E109_input_search_rho':
        g = pd.read_csv(res / 'gate.csv')
        return {f"{r.which}:{r.reference}:s{r.size:g}": r.cls for r in g.itertuples()}
    if experiment == 'E108_paa_extrapolation':
        g = pd.read_csv(res / 'gate.csv').rename(columns={'class': 'cls'})
        return {f'{r.population}:{r.input}': r.cls for r in g.itertuples()}
    v = json.loads((res / 'verdict_rho.json').read_text())
    out = {'B2': v.get('cls')}
    if v.get('draw2'):
        out['draw2'] = v['draw2'].get('cls')
    return out


def _cls(g: Mapping) -> Optional[str]:
    return g.get('cls') if g.get('complete', True) else None


def e109_holds(full: Mapping[str, Mapping[str, dict]], readings: Sequence[str] = REREAD) -> dict:
    """The entry's "E109's gain holds like for like" per reading (the module docstring): E109's carrying gate (the six
    at 0.6 against the lab's input) re-read with the ρ's under the reading; per condition the re-read mean gain against
    E109_HALF of the gain as written, each control below CONTROL_RATIO of the re-read gain."""
    from experiments import E109_input_search_rho as X
    key = 'six:I0:s0.6'
    g0 = full['as_written'][key]
    out = {'as_written': {f'mean_gain_{c}': g0.get(f'mean_gain_{c}') for c in CONDITIONS}}
    for r in readings:
        g = full[f'{r}:rho'][key]
        row = {}
        for c in CONDITIONS:
            gain, gain0 = g.get(f'mean_gain_{c}', np.nan), g0.get(f'mean_gain_{c}', np.nan)
            # (the gate reads its controls only on a ρ that gains; read here on both, as its gate computes them)
            vals = _e109_controls(full[f'{r}:rho:rd'], X, c)
            ctl_ok = bool(np.isfinite(gain) and gain > 0 and all(v is not None and np.isfinite(v) and v < CONTROL_RATIO * gain
                                                                for v in vals.values()))
            row.update({f'mean_gain_{c}': gain, f'mean_gain_{c}_as_written': gain0, f'half_{c}': E109_HALF * gain0,
                        f'gain_holds_{c}': bool(np.isfinite(gain) and gain >= E109_HALF * gain0),
                        f'controls_{c}': vals, f'controls_ok_{c}': ctl_ok,
                        f'holds_{c}': bool(np.isfinite(gain) and gain >= E109_HALF * gain0 and ctl_ok)})
        row['holds'] = bool(row['holds_on'] and row['holds_off'])
        row['cls_rho'] = _cls(g)
        out[r] = row
    return out


def _e109_controls(rd: pd.DataFrame, X, c: str) -> dict:
    """E109's controls' mean gains on ρ `c` (permuted, inverted, control search), as its gate computes them."""
    runs = X.runs_table()
    vals = {}
    for label, kind, seeds, names in (('permuted', 'perm', X.PERM_SEEDS, X.FIRST_SPLIT), ('inverted', 'inverted', (None,), X.FIRST_SPLIT),
                                      ('control search', 'search', (None,), X.CONTROL_SEARCHES)):
        per = []
        for half in sorted({X.SEARCHES[n]['half'] for n in names}):
            dv = X.direction_values(rd, runs, [n for n in names if X.SEARCHES[n]['half'] == half], 'I0', kind, X.GATE_SIZE, seeds)
            per.append(np.nan if dv is None else dv[f'd_rho_{c}'])
        vals[label] = float(np.mean(per)) if per and np.isfinite(per).all() else np.nan
    return vals


def all_gates(jobs: int = 1, experiments: Sequence[str] = GATE_EXPERIMENTS) -> Tuple[pd.DataFrame, dict]:
    """gates.csv and gates.json: every experiment's gates as written and re-read, the identities, E109's verdict."""
    frames, full = [], {}
    for ex in experiments:
        df, f = reread_gate(ex, REREAD, jobs)
        sv = saved_classes(ex)
        df['class_saved'] = df['gate'].map(sv)
        frames.append(df)
        full[ex] = f
        if ex == 'E109_input_search_rho':
            for r in REREAD:
                f[f'{r}:rho:rd'] = substitute(saved_readings(ex), pd.read_csv(out_dir() / f'gate_readings_{ex}_{r}.csv'), ex, 'rho')
    gt = pd.concat(frames, ignore_index=True)
    return gt, full


# ---------------------------------------------------------------------------
# identities

def identities(fam: Mapping[str, Mapping[str, Dict[str, pd.Series]]], gates_full: Optional[Mapping] = None,
               gates_table: Optional[pd.DataFrame] = None, tol: float = 1e-9) -> dict:
    """The checks of part A (the module docstring); 'all_exact' False if any fails."""
    from experiments import E106_rho_gap_parts as E106
    from experiments import E72_wiring_population as E72
    from experiments import E73_per_kc_normalization as E73
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E93_missing_glomeruli_input as E93
    out = {}
    D = data_matrices()
    r = E73.refs()['graded']
    out['data_vs_e73_graded'] = max(float((D[c] - _cs(r[c]).reindex(D[c].index)).abs().max()) for c, k in (('on', 'on'), ('off', 'off')))
    # R0: E106's family population ρ's (conditions.csv, 120 pairs) and each wiring's
    cond = pd.read_csv(H.results_dir('E106_rho_gap_parts') / 'conditions.csv')
    cond = cond[(cond['family'] == 'family') & (cond['subset'] == 'all120')].set_index('reading')
    U0 = unit_distances(fam['R0'])
    lab = {'population': 'all8'}
    out['r0_vs_e106'] = max(abs(H.spearman(U0[lab.get(k, k)][c], D[c]) - float(cond.loc[k, f'rho_{c}']))
                            for k in cond.index for c in CONDITIONS)
    # R1: E72's 'cal' distances of each wiring
    out['r1_vs_e72_cal'] = max(float((fam['R1'][w][c] - _cs(E72.distances(E72.load_fit(E83.saved_fit(w)['folder'],
                                E83.saved_fit(w)['tag']), 'cal')[c]).reindex(fam['R1'][w][c].index)).abs().max())
                               for w in ('real', 'P3') for c in CONDITIONS)
    # R2: E93's zcell_distances, and E73's reading of C2 (by_variant.csv, zcell_saturating, graded)
    o = E83.saved_fit('real')['saved']
    z93 = E93.zcell_distances(o)
    out['r2_vs_e93'] = max(float((fam['R2']['real'][c] - z93[c].reindex(fam['R2']['real'][c].index)).abs().max()) for c in CONDITIONS)
    bv = pd.read_csv(H.results_dir('E73_per_kc_normalization') / 'by_variant.csv', low_memory=False)
    m = bv[(bv['folder'] == 'E67_apl_effect_fits') & (bv['tag'] == E83.saved_fit('real')['tag']) &
           (bv['reading_def'] == 'zcell_saturating') & (bv['reading'] == 'graded')]
    if len(m):
        out['r2_vs_e73_rho'] = max(abs(H.spearman(fam['R2']['real'][c], D[c]) - float(m[f'rho_{c}'].iloc[0])) for c in CONDITIONS)
        out['r2_vs_e73_loss'] = abs(reading_scores(fam['R2']['real'])['loss'] - float(m['rmse'].iloc[0]))
    # R4: the pipeline's deterministic distances are R0; R4 of the four reproduces E72's partE_noise.csv (raw, graded)
    out['r4_det_vs_r0'] = max(float((fam['R4det'][w][c] - fam['R0'][w][c].reindex(fam['R4det'][w][c].index)).abs().max())
                              for w in fam['R4det'] for c in CONDITIONS)
    pe = pd.read_csv(H.results_dir(PART_E) / 'partE_noise.csv')
    pe = pe[(pe['group'] == 'population') & (pe['transform'] == 'graded') & (pe['distances'] == 'raw')]
    pl = pd.read_csv(H.results_dir(PART_E) / 'partB' / 'pairs_noise_long.csv')
    pl = _levels(pl)
    dd = []
    for level in ('match_on', 'match_off'):
        g = pl[(pl['level'] == level) & (pl['transform'] == 'graded')]
        on = g[g.condition == 'on'].groupby(['odor_a', 'odor_b']).mean(numeric_only=True)
        off = g[g.condition == 'acute_off'].groupby(['odor_a', 'odor_b']).mean(numeric_only=True)
        st = H.apl_effect_stats(on['raw'], off['raw'], H.apl_effect(on['data_raw'], off['data_raw']))
        dd.append(abs(st['apl_effect_rmse'] - float(pe.loc[pe['level'] == level, 'apl_effect_rmse'].iloc[0])))
    out['r4_partE_rebuild'] = max(dd)
    S = population([fam['R4'][w] for w in E83.S_WIRINGS])
    out['r4_data_raw_vs_data'] = max(float((_cs(pd.Series(pl[(pl['condition'] == cc) & (pl['level'] == 'match_on')
                                        & (pl['transform'] == 'graded')].groupby(['odor_a', 'odor_b'])['data_raw'].mean()))
                                        - D[c]).abs().max()) for c, cc in (('on', 'on'), ('off', 'acute_off')))
    out['r4_S_rows'] = len(S['on'])
    checks = {k: v for k, v in out.items() if k not in ('r4_S_rows', 'r4_data_raw_vs_data')}
    if gates_full is not None:
        for ex, f in gates_full.items():
            for k, v in f['identities'].items():
                out[f'{ex}:{k}'] = v
                checks[f'{ex}:{k}'] = v['max_abs_diff']
    if gates_table is not None:
        aw = gates_table[gates_table['reading'] == 'as written']
        mism = aw[aw['class_saved'].notna() & (aw['cls'].astype(str) != aw['class_saved'].astype(str))]
        out['classes_as_written_mismatches'] = mism[['experiment', 'gate', 'cls', 'class_saved']].to_dict('records')
        checks['classes_as_written_mismatches'] = float(len(mism))
    # (E73's table rounds nothing, but its calibration converged on another process; allow a small tolerance there)
    tols = {'r2_vs_e73_rho': 1e-6, 'r2_vs_e73_loss': 1e-6}
    fails = {k: v for k, v in checks.items() if not (np.isfinite(v) and v <= tols.get(k, tol))}
    out['fails'] = fails
    out['all_exact'] = not fails
    return out


# ---------------------------------------------------------------------------
# the parts

def provenance(part: str, **extra) -> None:
    prov = dict(generated=datetime.now(timezone.utc).isoformat(), script=SCRIPT, part=part, host=os.uname().nodename,
                approx_conc_note=H.APPROX_CONC_NOTE, readings=READING_LABEL, n_pseudo=N_PSEUDO, n_trials=N_TRIALS,
                n_r3_draws=N_R3_DRAWS, **H.code_commits(), model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()),
                **extra)
    _write_json(out_dir() / f'provenance_{part}.json', prov)


def verdicts(rt: pd.DataFrame, holds: Optional[dict] = None) -> dict:
    """verdict.json: "the reading moves ρ" under every reading against R0 (R3 and R4 decide; R1, R2 and R4z beside),
    "E109's gain holds like for like" per reading, and the outcome letter."""
    mv = {r: moves_verdict(rt, r) for r in ('R1', 'R2', 'R3', 'R4', 'R4z') if (rt['reading'] == r).any()}
    moves = bool(any(mv[r].get('moves') for r in ('R3', 'R4') if r in mv))
    out = dict(moves=mv, reading_moves_rho=moves, deciding=('R3', 'R4'),
               provisional=bool(mv.get('R4', {}).get('provisional', True)))
    if holds is not None:
        out['e109_holds'] = holds
        dec = [r for r in ('R3',) if r in holds]
        e109 = bool(all(holds[r]['holds'] for r in dec)) if dec else None
        out['e109_gain_holds'] = e109
    sign = any(mv[r].get('sign_split') for r in ('R3', 'R4') if r in mv)
    if not moves:
        out['outcome'] = 'd' if sign else 'c'
    else:
        out['outcome'] = 'a' if out.get('e109_gain_holds') else 'b'
    out['sign_split'] = bool(sign)
    # (the outcome is final only once R4 covers the eight wirings, part B)
    out['outcome_provisional'] = out['provisional']
    return out


def part_a(jobs: int = 1, gates: bool = True) -> dict:
    """Part A (no simulations): the family's readings and tables, the gates re-read, the identities, the verdicts."""
    warnings.simplefilter('ignore', RuntimeWarning)
    out = out_dir()
    fam = family_readings(jobs)
    rt = rhos_table(fam)
    rt.to_csv(out / 'rhos.csv', index=False)
    changes_table(rt).to_csv(out / 'changes.csv', index=False)
    bound_table(fam).to_csv(out / 'bound.csv', index=False)
    yardstick_table(fam).to_csv(out / 'yardstick.csv', index=False)
    realism_table().to_csv(out / 'realism.csv', index=False)
    gt, full, holds = None, None, None
    if gates:
        gt, full = all_gates(jobs)
        holds = e109_holds(full['E109_input_search_rho'])
        gt.to_csv(out / 'gates.csv', index=False)
        _write_json(out / 'gates.json', {ex: {k: v for k, v in f.items() if not k.endswith(':rd')} for ex, f in full.items()})
    ident = identities(fam, full, gt)
    _write_json(out / 'identities.json', ident)
    v = verdicts(rt, holds)
    _write_json(out / 'verdict.json', dict(v, **H.code_commits()))
    provenance('A', all_exact=ident['all_exact'])
    show = rt[rt['unit'].isin(['all8', 'S', 'T'])][['reading', 'unit', 'rho_on_all120', 'rho_off_all120', 'paa_bmyr_on',
                                                     'paa_bmyr_off', 'loss', 'pattern', 'mean_effect']]
    print(show.round(4).to_string(index=False), flush=True)
    print(json.dumps({k: v[k] for k in ('reading_moves_rho', 'outcome', 'provisional', 'sign_split')}, default=str), flush=True)
    if not ident['all_exact']:
        raise RuntimeError(f"E112 part A: identities failed: {ident['fails']}")
    return v


# --- part B

def _e109_row(change: str, wiring: str) -> pd.Series:
    r = pd.read_csv(H.results_dir('E109_input_search_rho') / 'runs_R.csv')
    m = r[(r['change'] == change) & (r['wiring'] == wiring)]
    if len(m) != 1:
        raise ValueError(f'E109 runs_R.csv: {len(m)} rows for {change} on {wiring}')
    return m.iloc[0]


def part_b_tag(change: str, wiring: str) -> str:
    return f'E109-{change}_{wiring}'


def part_b_prepare() -> pd.DataFrame:
    """Part B's model table: results/E112_rho_like_for_like/fits.csv with one row per E109 input (E109_INPUTS_B at
    E109_SIZE) on each held-out wiring, the wiring's E72 fits.csv row with the input file (orn_input) and the held
    reading's absolute parameters (runs_R.csv's params_on) and the run's own response rates; each run's saved outputs
    linked under outputs/E112_rho_like_for_like/ by the row's tag, so that E43's pipeline asserts its deterministic run
    against them."""
    from experiments import E109_input_search_rho as X
    from experiments import E83_cross_wiring_fit as E83
    f72 = pd.read_csv(H.results_dir('E72_wiring_population') / 'fits.csv')
    rows = []
    loc = local_dir()
    for w in E83.T_WIRINGS:
        base = f72[f72['tag'] == E83.saved_fit(w)['tag']].iloc[0]
        for name in E109_INPUTS_B:
            ch = f'{name}_s{E109_SIZE:g}'
            r = _e109_row(ch, w)
            p = json.loads(r['params_on'])
            o = X.load_run('R', f'{ch}_{w}')
            if o is None:
                raise RuntimeError(f'no saved E109 run {ch}_{w}')
            row = base.copy()
            tag = part_b_tag(ch, w)
            row['tag'] = tag
            row['variant'] = tag
            row['orn_input'] = str(r['orn_input'])
            for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale'):
                row[k] = float(p[k])
            row['kc_kc_inh'] = float(p['wKCKC_scale'])
            row['rate_on'] = float(o['on'][0].mean().mean())
            row['rate_off'] = float(o['off'][0].mean().mean())
            row['init_from'] = f'E109_input_search_rho runs_R.csv {ch} on {w} (held reading)'
            rows.append(row)
            for c, cc in (('on', 'on'), ('off', 'acute_off')):
                for kind in ('responses', 'spike_counts'):
                    src = next(d / f'{ch}_{w}_{cc}_{kind}.parquet' for d in (X.local_dir('R'),
                               H.shared_dir() / 'data' / 'model_outputs' / X.EXPERIMENT / 'R')
                               if (d / f'{ch}_{w}_{cc}_{kind}.parquet').exists())
                    dst = loc / f'{tag}_{cc}_{kind}.parquet'
                    if dst.is_symlink() or dst.exists():
                        dst.unlink()
                    dst.symlink_to(src)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'fits.csv', index=False)
    print('\n'.join(part_b_commands()), flush=True)
    return df


def part_b_check(tags: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Builds part B's models of fits.csv (`tags`, default all) as E43's pipeline builds them and compares each
    deterministic run with the saved held reading entry for entry (part_b_check.csv; E43 asserts the same when its lanes
    start, this finds a mismatch before a lane is launched). One model build and one deterministic run per tag."""
    from experiments import E43_model_trial_noise as E43
    f = pd.read_csv(out_dir() / 'fits.csv')
    tags = list(f['tag']) if tags is None else list(tags)
    orn = H.apl16_orn_deltas()
    wd = local_dir() / 'partB_check'
    wd.mkdir(exist_ok=True)
    rows = []
    for t in tags:
        name = E43.model_from(f'{EXPERIMENT}:{t}')
        m = E43.NoisyModel(name, orn, wd)
        row = dict(tag=t, model=name)
        for cond in ('on', 'acute_off'):
            r, s = (pd.read_parquet(local_dir() / f'{t}_{cond}_{k}.parquet') for k in ('responses', 'spike_counts'))
            pr, ps = E43.to_arrays(r, s)
            row[f'{cond}_responses_mismatch'] = int((pr != m.det[cond][0]).sum())
            row[f'{cond}_spike_counts_mismatch'] = int((ps != m.det[cond][1]).sum())
        rows.append(row)
        print(row, flush=True)
        m.close()
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'part_b_check.csv', index=False)
    return df


def part_b_commands() -> List[str]:
    """E43's part B per held-out wiring (one lane each, 16 threads): the wiring's fit and E109's four inputs on it, at
    E72's part E settings, into results/E112_rho_like_for_like/partB_<wiring>/."""
    from experiments import E83_cross_wiring_fit as E83
    cmds = []
    for w in E83.T_WIRINGS:
        specs = [f"--model-from E72_wiring_population:{E83.saved_fit(w)['tag']}"]
        specs += [f'--model-from {EXPERIMENT}:{part_b_tag(f"{n}_s{E109_SIZE:g}", w)}' for n in E109_INPUTS_B]
        cmds.append('OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E43_model_trial_noise.py '
                    f'--part B --tune-trials 3 --seed 43 --n-recordings 10 --n-trials 3 --levels match_on match_off '
                    f'{" ".join(specs)} --experiment {EXPERIMENT} '
                    f'--results-dir {H.results_dir(EXPERIMENT)}/partB_{w} --outputs-dir {H.outputs_dir(EXPERIMENT)}/partB_{w}')
    return cmds


def part_b_inputs_r4() -> pd.DataFrame:
    """R4 of E109's inputs on T against the wiring's fit through the same pipeline: per input and unit (each wiring of T,
    T's population) the two ρ's against all recordings and the input's unseen half (E108's `dist_targets`), and the change
    from the fit's (the gain under R4)."""
    from experiments import E108_paa_extrapolation as E108
    from experiments import E109_input_search_rho as X
    from experiments import E83_cross_wiring_fit as E83
    dt = E108.dist_targets()
    per: Dict[str, Dict[str, Dict[str, pd.Series]]] = {}
    for w in E83.T_WIRINGS:
        p = out_dir() / f'partB_{w}' / 'pairs_noise_long.csv'
        if not p.exists():
            continue
        R = trial_noise_distances(pd.read_csv(p))
        for m, d in R.items():
            if m.startswith('E112-E109-'):
                ch = m[len('E112-E109-'):].rsplit('_', 1)[0]
            elif _model_wiring(m) == w:
                ch = 'I0'
            else:
                continue
            per.setdefault(ch, {})[w] = d
    rows = []
    from experiments import E106_rho_gap_parts as E106
    for ch, ws in per.items():
        name = ch.split('_s')[0]
        unseen = X.UNSEEN.get(X.SEARCHES.get(name, {}).get('half'), None)
        units = dict(ws)
        if all(w in ws for w in E83.T_WIRINGS):
            units['T'] = population([ws[w] for w in E83.T_WIRINGS])
        for u, d in units.items():
            row = dict(change=ch, unit=u, unseen=unseen)
            for c in CONDITIONS:
                row[f'rho_{c}_all'] = float(E106.rho_pairs(d[c], dt['all'][c]))
                for h in ('A', 'B', 'C', 'D'):
                    row[f'rho_{c}_{h}'] = float(E106.rho_pairs(d[c], dt[h][c]))
            rows.append(row)
    df = pd.DataFrame(rows)
    if len(df) and (df['change'] == 'I0').any():
        base = df[df['change'] == 'I0'].set_index('unit')
        for c in CONDITIONS:
            df[f'd_rho_{c}_unseen'] = [row[f'rho_{c}_{row["unseen"]}'] - base.loc[row['unit'], f'rho_{c}_{row["unseen"]}']
                                       if isinstance(row['unseen'], str) and row['unit'] in base.index else np.nan
                                       for _, row in df.iterrows()]
    return df


def part_b_score(jobs: int = 1) -> dict:
    """Part B's scores: part A again (R4 now on every wiring the pipeline ran, the verdicts final when it covers the
    eight) and R4 of E109's inputs on T (inputs_r4.csv)."""
    v = part_a(jobs)
    df = part_b_inputs_r4()
    df.to_csv(out_dir() / 'inputs_r4.csv', index=False)
    provenance('B_score')
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B-prepare', 'B-check', 'B-score'), required=True)
    ap.add_argument('--tags', nargs='*', default=None, help='B-check: the fits.csv tags to build (default all)')
    ap.add_argument('--jobs', type=int, default=1, help='processes for the per-run readings (one OpenMP thread each)')
    ap.add_argument('--no-gates', action='store_true', help='part A without the gates re-read (a quick look)')
    a = ap.parse_args()
    if a.part == 'A':
        part_a(a.jobs, gates=not a.no_gates)
    elif a.part == 'B-prepare':
        part_b_prepare()
    elif a.part == 'B-check':
        part_b_check(a.tags)
    else:
        part_b_score(a.jobs)


if __name__ == '__main__':
    main()
