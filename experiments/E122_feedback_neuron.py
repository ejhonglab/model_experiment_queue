#!/usr/bin/env python3
"""E122: one feedback neuron that survives the block in place of KC-KC inhibition: the family's fitted KC-KC matrix
replaced by a rank-1 (or per-type block rank-1) matrix of the same total weight, through olfsysm's existing KC-KC path.

The entry is item E122 of experiment_priorities.md (proposed by viki 2026-11-20 22:09; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-20 21:59; moved into the list 2026-11-21 00:12 under the idle-host rule).
"The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's
`saved_fit`); S is real, P1, P2 and M1 and T is P3, P4, M2 and M3 (E83's S_WIRINGS / T_WIRINGS); "the population" of a
set of wirings is the mean over them of the pair distances per condition (E72's, `E83.population`). The family's bouton
base pins multibouton_claws='first' (E3's BOUTON, through E10_joint_fit.base_kws), so its fits are reproduced although
olfsysm's default is now 'argmax'.

The matrices (KC x KC, post x pre, the model's KCs; `arm_matrix`). The family's KC-KC input is the fitted scale times the
hemibrain KC->KC synapse counts normalized to a mean row sum of 1 (`hemibrain_kckc`, fit_mb_model's own). Every arm's
matrix has the same total weight (the sum of all entries, N, the number of KCs), so at the fit's KC-KC scale it delivers
the same total inhibition for uniform activity, and the scale stays comparable:
- U (DPM-like): one neuron, all KCs in and all KCs out, uniform: every entry 1/N (rank 1, the diagonal included: a pooled
  neuron feeds back onto the KCs that drive it); every KC's row sum 1.
- TY: one neuron per KC type (γ 'g', α′β′ "a'b'", αβ 'ab'), each pooling its own type's KCs uniformly and inhibiting them
  uniformly: within type t every entry R_t / n_t, R_t the hemibrain matrix's mean row sum over type t's KCs (so each
  type's total incoming weight is the hemibrain's), 0 between types (block-diagonal, rank 1 per block). ("TY", not "T",
  which names the held-out wirings.)
- G (MB-C1-like): all KCs in, the γ KCs out: rows of γ KCs every entry 1 / n_γ (row sum N / n_γ), other rows 0.
- SH (control): the hemibrain matrix shuffled within each (post type, pre type) block: the block's entries permuted
  among its positions (the diagonal excluded, kept 0), seed SH_SEED. Each block keeps its weights and total; who
  inhibits whom within the block, and each KC's own amount, are drawn again.
- GAB (control): G with its output moved to the αβ KCs (rows of αβ KCs 1 / n_αβ).
The wirings re-draw the PN->claw wiring of the same hemibrain KCs (E72), so every wiring has the same KCs and one matrix
per arm serves all (checked in step 0).

Parts (CLI `--part`):
- 0 (no simulation): (i) the hemibrain matrix's weight by (post type, pre type) block: share of the total weight, share
  of the KC pairs, mean row sum, and the within-type share against the share of KC pairs within type; (ii) how well each
  KC's received KC-KC input in the family's saved runs (the input the fitted matrix delivers, W @ spike counts, per KC
  and odor, APL on and APL off; the counts are the integral of the activity the matrix filters, so this is the total
  input over the run) is predicted by pooled activity: the arms' matrices applied to the same counts ('U', 'TY', 'G',
  'GAB'), the block-mean matrix ('blocks', SH's expectation) and each KC's own per-pre-type amounts times each type's
  mean count ('own_amounts'). Per wiring and condition: R² (1 - SSE / SST, no free scale), Pearson r over KC x odor,
  and the mean over odors of Pearson r across KCs. step0_blocks.csv, step0_prediction.csv, step0.json.
- check: (i) the hemibrain matrix passed explicitly as `_wKCKC` at the fit's parameters reproduces the hemibrain fit's
  saved APL-on and APL-off spike counts exactly; (ii) a 5-evaluation smoke fit of E10_joint_fit.py with `--kckc-matrix`
  (arm U, hemibrain wiring). check.json.
- 1 (simulations): each arm on S's four wirings at the wiring's fitted scales (KC-KC at the fit's scale, the arm's
  matrix), both response rates held at the fit's (`E77_input_screen.hold_rates`: the somatic threshold bisected with
  APL off, then the tied APL scales with APL on; APL->bouton as fitted). Rows in runs_1.csv (absolute parameters, the
  factors, realism, the types' response rates, commits), runs under outputs/E122_feedback_neuron/step1/, each held point
  in starts.json as S1_<arm>_<w>.
- score1 (no simulation): step-1 readings (each wiring and S's population; raw counts and the calibrated readout) beside
  the family and E110's part A (KC-KC 0 at fixed scales) on the same wirings; the gate: an arm goes on if the pattern
  part under the calibrated readout of its S population is GATE1 (0.055) or lower (the entry's number; the family has
  0.050 there and E110's part A 0.061, read from E110's readings.csv and recorded beside). gate1.json, step1_readings.csv.
- 2 (fits): each arm that passed, on S's four wirings: E10_joint_fit.py with E110's / E72's part-B command, the arm's
  matrix (`--kckc-matrix`), the threshold, APL->KC (with KC->APL and bouton->APL tied), APL->bouton (cap 3) and the
  neuron's scale (the KC-KC scale) free, started at the step-1 held point, 150 evaluations, each in its own folder
  (results and outputs of E122_feedback_neuron/fits/<arm>_<w>). Main arms first (U, TY, G, ordered by their step-1
  pattern part), then the controls.
- T (simulations, after an arm's four S fits): the read on T. The refits' carried parameters: the neuron's scale (the
  mean of the four S refits' KC-KC scales, absolute) and the APL->bouton factor (the mean over S of the refit's
  APL->bouton over the family fit's); on each T wiring the family's fit with these, the arm's matrix, and both response
  rates re-held (hold_rates: threshold, then the tied APL scales). Runs under outputs/.../T/, rows in runs_T.csv.
- score (no simulation): readings.csv (every kind of run: the family, E110's A and B, step 1, the refits and the T
  reads; each wiring and the S and T populations; raw and calibrated: E108's `unit_scores` against all recordings and
  halves A to D, the two ρ's, the runs' response rates, the γ / α′β′ / αβ response rates and off / on ratios, bands),
  verdict.json, comparison_table.md, provenance_score.json. The class of each refitted arm, on T against the family's
  T fits: "replaces KC-KC" if its pattern part under the readout is not higher than the family's by more than WITHIN
  (0.003) on T's population and on 3 of the 4 T wirings, against all recordings and each half (A to D), with every T
  run inside the four bands; "below the family" (E121's pass rule) if it is lower by BELOW (0.003) or more on the
  population and on 3 of 4 wirings against all recordings and each half, the mean APL effect under the readout not
  lower by 0.01 and the raw-count loss not higher by 0.005 (population), every T run in the bands; "neither" otherwise.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.
- chain: check, 1, score1, 2 and T for the passing arms (gated), score, variants, in one process (the detached launch).

Choices fixed before any run: the matched total (the sum of the matrix, not each KC's row sum: U, TY and G give KCs the
same input whatever their hemibrain amount, which is what separates them from SH); the read on T by carried parameters
with the rates re-held (the T wirings have their own family fits, which the carried parameters adjust); the gate's
0.055 as the entry states it; the halves A to D of E108's `effect_targets`.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws
of it, and the model loader's hemibrain KC->KC counts and KC types; no other connectome table; no raw-only recording.
Nothing is published to claude.ai.

Usage (from the repository root; 16 threads, under nice, one model job at a time):
    ./venv/bin/python -u experiments/E122_feedback_neuron.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E122_feedback_neuron.py --part chain
then
    ./venv/bin/python archive_tools/archive_model_outputs.py E122_feedback_neuron --exclude 'matrices,distances_cache.pkl'
    ./venv/bin/python archive_tools/migrate_model_output_params.py E122_feedback_neuron --apply
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import socket
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108

EXPERIMENT = 'E122_feedback_neuron'
ARMS = ('U', 'TY', 'G', 'SH', 'GAB')
MAIN_ARMS = ('U', 'TY', 'G')
ARM_LABELS = {'U': 'one neuron, all KCs in and out, uniform (DPM-like)',
              'TY': 'one neuron per KC type, each pooling and inhibiting its own type',
              'G': 'all KCs in, the gamma KCs out (MB-C1-like)',
              'SH': 'the hemibrain matrix shuffled within each type pair (control)',
              'GAB': "G's output moved to the alpha-beta KCs (control)"}
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
S_WIRINGS = E83.S_WIRINGS
T_WIRINGS = E83.T_WIRINGS
WIRINGS = E83.ALL_WIRINGS
K = 2.0
MAXFEV = 150
SH_SEED = 0
GATE1 = 0.055
WITHIN, BELOW = 0.003, 0.003
MEAN_SLACK, RAW_SLACK = 0.01, 0.005
N_T_NEEDED = 3
TARGETS = ('all',) + tuple(E108.HALVES)
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS
FIT_PREFIX = f'{EXPERIMENT}/fits'
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E122_feedback_neuron/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=str) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _rows(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame(columns=['tag'])


# ---------------------------------------------------------------------------
# the matrices

def model_kc_index(w: str = 'real') -> pd.MultiIndex:
    """The model's KC index (kc_id, kc_type) of a wiring's saved fit."""
    return E83.saved_fit(w)['saved']['on'][1].index


def kc_types_of(kc_index: pd.Index) -> np.ndarray:
    """The KC type ('g', "a'b'", 'ab') of each model KC (the index's kc_type level)."""
    return np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)


def hemibrain_kckc(kc_ids: Sequence[int]) -> pd.DataFrame:
    """The normalized hemibrain KC x KC matrix (post x pre, mean row sum 1), as fit_mb_model builds it."""
    from pyolfsysm.mb_model import connectome_KCKC_weights
    return connectome_KCKC_weights(pd.Index(kc_ids))


def shuffle_within_blocks(ref: np.ndarray, types: np.ndarray, seed: int = SH_SEED) -> np.ndarray:
    """`ref` with each (post type, pre type) block's entries permuted among the block's positions (the diagonal of a
    same-type block excluded and kept as it is)."""
    rng = np.random.default_rng(seed)
    out = ref.copy()
    for tp in TYPES:
        rows = np.flatnonzero(types == tp)
        for tq in TYPES:
            cols = np.flatnonzero(types == tq)
            ii, jj = np.meshgrid(rows, cols, indexing='ij')
            ii, jj = ii.ravel(), jj.ravel()
            if tp == tq:
                keep = ii != jj
                ii, jj = ii[keep], jj[keep]
            vals = ref[ii, jj]
            out[ii, jj] = vals[rng.permutation(len(vals))]
    return out


def arm_array(arm: str, ref: np.ndarray, types: np.ndarray, seed: int = SH_SEED) -> np.ndarray:
    """The arm's KC x KC matrix (post x pre) as an array with the same total weight as `ref` (the module docstring)."""
    n = len(types)
    total = float(ref.sum())
    if arm == 'HB':
        return ref.copy()
    if arm == 'U':
        return np.full((n, n), total / n / n)
    if arm == 'TY':
        out = np.zeros((n, n))
        for t in TYPES:
            m = types == t
            out[np.ix_(m, m)] = ref[m].sum() / m.sum() / m.sum()
        return out
    if arm in ('G', 'GAB'):
        m = types == ('g' if arm == 'G' else 'ab')
        out = np.zeros((n, n))
        out[m, :] = total / m.sum() / n
        return out
    if arm == 'SH':
        return shuffle_within_blocks(ref, types, seed)
    raise ValueError(f'unknown arm {arm!r}')


def arm_matrix(arm: str, ref: pd.DataFrame, types: np.ndarray, seed: int = SH_SEED) -> pd.DataFrame:
    """The arm's KC x KC matrix (post x pre, `ref`'s index and columns) with the same total weight as `ref`."""
    return pd.DataFrame(arm_array(arm, ref.to_numpy(dtype=float), types, seed), index=ref.index, columns=ref.columns)


def matrix_path(arm: str) -> Path:
    return local_dir('matrices') / f'kckc_{arm}.parquet'


def write_matrix(arm: str) -> Tuple[Path, pd.DataFrame, str]:
    """The arm's matrix for the model's KCs, written to outputs/.../matrices/kckc_<arm>.parquet (columns as strings, as
    parquet needs; E10_joint_fit.py's --kckc-matrix reads them back as integers) unless there with the same values;
    returns (path, matrix, sha256 of its values)."""
    from experiments.E77_input_screen import frame_sha256
    idx = model_kc_index()
    ids = idx.get_level_values(H.KC_ID_LEVEL)
    ref = hemibrain_kckc(ids)
    M = arm_matrix(arm, ref, kc_types_of(idx))
    sha = frame_sha256(M)
    p = matrix_path(arm)
    if p.exists():
        old = pd.read_parquet(p)
        old.columns = pd.Index([int(c) for c in old.columns], name=old.columns.name)
        if frame_sha256(old) == sha:
            return p, M, sha
    W = M.copy()
    W.columns = W.columns.astype(str)
    W.to_parquet(p)
    return p, M, sha


def load_matrix(arm: str) -> pd.DataFrame:
    return write_matrix(arm)[1]


# ---------------------------------------------------------------------------
# step 0

def block_table(ref: np.ndarray, types: np.ndarray) -> pd.DataFrame:
    """Per (post type, pre type): the share of the matrix's total weight, the share of the KC pairs (diagonal
    excluded), the mean row sum of the post type's KCs from the pre type, and the number of nonzero entries."""
    n = len(types)
    tot = float(ref.sum())
    rows = []
    for tp in TYPES:
        mp = types == tp
        for tq in TYPES:
            mq = types == tq
            blk = ref[np.ix_(mp, mq)]
            pairs = mp.sum() * mq.sum() - (mp.sum() if tp == tq else 0)
            rows.append(dict(post_type=tp, pre_type=tq, n_post=int(mp.sum()), n_pre=int(mq.sum()),
                             weight_share=float(blk.sum() / tot), pair_share=float(pairs / (n * (n - 1))),
                             mean_row_sum=float(blk.sum(axis=1).mean()), nonzero=int((blk > 0).sum())))
    return pd.DataFrame(rows)


def predictors(ref: np.ndarray, types: np.ndarray) -> Dict[str, np.ndarray]:
    """The matrices whose input predicts the hemibrain matrix's: the arms U, TY, G, GAB, the block-mean matrix
    ('blocks', each block's entries replaced by its mean over its positions, the diagonal of same-type blocks 0) and
    'own_amounts' (each KC's own row sum from each pre type, spread evenly over that type's KCs)."""
    n = len(types)
    out = {a: arm_array(a, ref, types) for a in ('U', 'TY', 'G', 'GAB')}
    blocks = np.zeros((n, n))
    own = np.zeros((n, n))
    for tp in TYPES:
        mp = np.flatnonzero(types == tp)
        for tq in TYPES:
            mq = np.flatnonzero(types == tq)
            blk = ref[np.ix_(mp, mq)]
            npos = blk.size - (len(mp) if tp == tq else 0)
            blocks[np.ix_(mp, mq)] = blk.sum() / npos if npos else 0.0
            own[np.ix_(mp, mq)] = blk.sum(axis=1, keepdims=True) / len(mq)
    np.fill_diagonal(blocks, 0.0)
    out['blocks'] = blocks
    out['own_amounts'] = own
    return out


def prediction_scores(actual: np.ndarray, pred: np.ndarray) -> dict:
    """R² (1 - SSE / SST about the actual's mean; no free scale), Pearson r over all entries (KC x odor), and the mean
    over odors (columns) of Pearson r across KCs (odors whose actual or predicted input is constant are skipped)."""
    a, p = actual.ravel(), pred.ravel()
    sst = float(((a - a.mean()) ** 2).sum())
    r2 = float(1 - ((a - p) ** 2).sum() / sst) if sst > 0 else np.nan
    r = float(np.corrcoef(a, p)[0, 1]) if a.std() > 0 and p.std() > 0 else np.nan
    per = [np.corrcoef(actual[:, j], pred[:, j])[0, 1] for j in range(actual.shape[1])
           if actual[:, j].std() > 0 and pred[:, j].std() > 0]
    return dict(r2=r2, pearson_r=r, mean_r_within_odor=float(np.mean(per)) if per else np.nan, n_odors_within=len(per))


def part_0() -> Dict[str, pd.DataFrame]:
    """Step 0 (the module docstring): step0_blocks.csv, step0_prediction.csv, step0.json."""
    out = out_dir()
    idx = model_kc_index()
    for w in WIRINGS:
        if not model_kc_index(w).equals(idx):
            raise RuntimeError(f'{w}: the KCs differ from the hemibrain wiring\'s (one matrix per arm would not serve)')
    types = kc_types_of(idx)
    ref = hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL)).to_numpy(dtype=float)
    bt = block_table(ref, types)
    bt.to_csv(out / 'step0_blocks.csv', index=False)
    within = float(bt.loc[bt.post_type == bt.pre_type, 'weight_share'].sum())
    within_pairs = float(bt.loc[bt.post_type == bt.pre_type, 'pair_share'].sum())
    rowsum = ref.sum(axis=1)
    preds = predictors(ref, types)
    rows = []
    for w in WIRINGS:
        F = E83.saved_fit(w)
        scale = float(F['p']['wKCKC_scale'])
        for c in ('on', 'off'):
            S = F['saved'][c][1].to_numpy(dtype=float)
            actual = scale * (ref @ S)
            for name, P in preds.items():
                rows.append(dict(wiring=w, condition=c, predictor=name, kckc_scale=scale,
                                 **prediction_scores(actual, scale * (P @ S))))
    pr = pd.DataFrame(rows)
    pr.to_csv(out / 'step0_prediction.csv', index=False)
    mean = pr.groupby(['condition', 'predictor'])[['r2', 'pearson_r', 'mean_r_within_odor']].mean()
    summary = dict(
        generated=_now(), host=_host(), n_kcs=int(len(types)), n_by_type={t: int((types == t).sum()) for t in TYPES},
        total_weight=float(ref.sum()), nonzero=int((ref > 0).sum()), within_type_weight_share=within,
        within_type_pair_share=within_pairs, within_over_pairs=within / within_pairs,
        mean_row_sum_by_type={t: float(rowsum[types == t].mean()) for t in TYPES},
        row_sum_cv=float(rowsum.std() / rowsum.mean()), kcs_without_input=int((rowsum == 0).sum()),
        prediction_mean_over_wirings={f'{c}/{p}': {k: float(v) for k, v in r.items()} for (c, p), r in mean.iterrows()},
        note=("received input = the wiring's fitted KC-KC scale x the normalized hemibrain matrix @ the saved spike counts "
              "(per KC and odor; the integral of the filtered activity the matrix acts on); predictors as "
              "`predictors`; R² without a free scale"),
        provenance=H.code_commits())
    _write_json(out / 'step0.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('provenance', 'prediction_mean_over_wirings')},
                     indent=1), flush=True)
    print(mean.round(3).to_string(), flush=True)
    return dict(blocks=bt, prediction=pr)


# ---------------------------------------------------------------------------
# runs: saving, loading, per-type rates

def type_rates(o: Dict[str, tuple]) -> dict:
    """Per KC type: the response rate (two or more spikes) with APL on and off and their ratio."""
    row = {}
    for c in ('on', 'off'):
        r = o[c][0]
        t = kc_types_of(r.index)
        R = r.to_numpy(dtype=float) > 0.5
        for tp in TYPES:
            row[f'rate_{c}_{TYPE_COL[tp]}'] = float(R[t == tp].mean())
    for tp in TYPES:
        on = row[f'rate_on_{TYPE_COL[tp]}']
        row[f'ratio_{TYPE_COL[tp]}'] = row[f'rate_off_{TYPE_COL[tp]}'] / on if on > 0 else np.nan
    return row


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, archive_dir() / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def hold_point(part: str, tag: str, arm: str, w: str, params: Mapping[str, float], table: Path, prov_path: Path,
               cc: dict, note: str) -> dict:
    """One fixed-scale run of an arm on a wiring with both response rates held at the wiring's fit's
    (E77_input_screen.hold_rates from `params`), saved and recorded (a row in `table`, provenance in `prov_path`)."""
    from experiments.E77_input_screen import hold_rates
    E108._still_clean(cc)
    P = E83.wiring_parts(w)
    path, M, sha = write_matrix(arm)
    res = hold_rates(lambda: E83.build_simulator(w, params=dict(params), extra_kws=dict(_wKCKC=M)),
                     target_off=P['target_off'], target_on=P['target_on'], params=dict(params), label=f'E122 {tag}')
    o = {'on': res['on'], 'off': res['off']}
    save_run(part, tag, o)
    row = dict(part=part, tag=tag, arm=arm, wiring=w, fit_folder=P['folder'], fit_tag=P['tag'], matrix=str(path.name),
               matrix_sha256=sha, start_params=json.dumps({k: float(v) for k, v in params.items()}),
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
               target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), host=_host())
    row.update(E83.realism(o))
    row.update(type_rates(o))
    row.update(cc)
    E83._append_row(table, row)
    prov = _read_json(prov_path) or dict(part=part, script='experiments/E122_feedback_neuron.py', points={})
    prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
        change=note, matrix=dict(arm=arm, label=ARM_LABELS[arm], path=str(path.relative_to(H.HERE)), sha256=sha),
        fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
    _write_json(prov_path, prov)
    print(f"E122 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, "
          f"{res['n_off']}/{res['target_off']} off); gamma {row['rate_on_g']:.3f} / {row['rate_off_g']:.3f}; flag "
          f"'{row.get('realism_flag') or ''}'; {res.get('seconds', 0):.0f}s", flush=True)
    del res, o
    gc.collect()
    return row


def starts_path() -> Path:
    return out_dir() / 'starts.json'


def add_start(key: str, params: Mapping[str, float], source: str) -> None:
    cur = _read_json(starts_path())
    e = {k: float(params[k]) for k in PARAM_KEYS}
    if key in cur and any(float(cur[key][k]) != e[k] for k in PARAM_KEYS):
        raise RuntimeError(f'starts.json: {key} would change (a start is fixed once written)')
    cur[key] = dict(e, source=source)
    _write_json(starts_path(), cur)


# ---------------------------------------------------------------------------
# check

def part_check(smoke: bool = True) -> dict:
    """(i) the explicit hemibrain matrix reproduces the hemibrain fit's saved runs; (ii) a smoke fit with --kckc-matrix."""
    out = out_dir()
    res = _read_json(out / 'check.json')
    if not res.get('explicit_matrix', {}).get('exact'):
        F = E83.saved_fit('real')
        _, M, sha = write_matrix('HB')
        sim, _ = E83.build_simulator('real', extra_kws=dict(_wKCKC=M))
        p = dict(F['p'])
        diff = {}
        for c, q in (('on', p), ('off', dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))):
            r, s = sim(q)
            diff[c] = int((s.to_numpy(dtype=float) != F['saved'][c][1].to_numpy(dtype=float)).sum())
        res['explicit_matrix'] = dict(differing_counts=diff, exact=all(v == 0 for v in diff.values()), matrix_sha256=sha,
                                      checked=_now())
        del sim
        gc.collect()
        _write_json(out / 'check.json', res)
        print(f"E122 check (i): {res['explicit_matrix']}", flush=True)
        if not res['explicit_matrix']['exact']:
            raise RuntimeError('the explicit hemibrain matrix does not reproduce the fit')
    if smoke and not res.get('smoke', {}).get('passes'):
        path, M, sha = write_matrix('U')
        add_start('fit_real', E83.saved_fit('real')['p'], "the hemibrain wiring's family fit (E83.saved_fit('real'))")
        spec = FitSpec('smokeU_real', 'U', 'real', ('params', 'fit_real'), 5)
        cc = fit_clean()
        rec = launch(fit_command(spec, smoke=True), H.results_dir(f'{FIT_PREFIX}/{spec.name}'),
                     local_dir('logs') / f'{spec.name}.log', cc, spec.name)
        row = fit_row(spec.name, fit_tag(spec))
        h = None
        hp = H.results_dir(f'{FIT_PREFIX}/{spec.name}')
        hs = sorted(hp.glob('history_*.csv'))
        if hs:
            h = pd.read_csv(hs[-1])
        res['smoke'] = dict(returncode=rec['returncode'], seconds=rec['seconds'], row=row is not None,
                            sha_in_row=(None if row is None else row.get('kckc_matrix_sha256')), matrix_sha256=sha,
                            kckc_values=(None if h is None or 'wKCKC_scale' not in h else
                                         sorted({float(x) for x in h['wKCKC_scale']})[:5]),
                            passes=bool(rec['returncode'] == 0 and row is not None
                                        and str(row.get('kckc_matrix_sha256')) == sha))
        _write_json(out / 'check.json', res)
        print(f"E122 check (ii): {res['smoke']}", flush=True)
        if not res['smoke']['passes']:
            raise RuntimeError('the smoke fit with --kckc-matrix failed')
    return res


# ---------------------------------------------------------------------------
# step 1

def part_1(arms: Sequence[str] = ARMS, wirings: Sequence[str] = S_WIRINGS) -> pd.DataFrame:
    """Step 1 (the module docstring), resumable: a point is done when its row is in runs_1.csv."""
    table = out_dir() / 'runs_1.csv'
    cc = E108._commits_clean()
    for arm in arms:
        for w in wirings:
            tag = f'S1_{arm}_{w}'
            done = _rows(table)
            if tag in set(done['tag']):
                print(f'E122 {tag}: done', flush=True)
                continue
            P = E83.wiring_parts(w)
            row = hold_point('step1', tag, arm, w, P['p'], table, out_dir() / 'provenance_1.json', cc,
                             note=(f"the KC-KC matrix replaced by arm {arm}'s (same total weight) at the fit's KC-KC scale; "
                                   "fixed_thr re-tuned with APL off to the fit's APL-off response count, then the tied APL "
                                   "scales to its APL-on count (E77_input_screen.hold_rates)"))
            add_start(tag, json.loads(row['params_on']), f'step 1 held point {tag} (results/{EXPERIMENT}/runs_1.csv)')
    return _rows(table)


# ---------------------------------------------------------------------------
# scoring helpers

def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83.pair_distances), cached by the sha256 of its spike counts and the readout."""
    import pickle
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
        st[key] = E83.pair_distances(o, readout)
        try:
            path.write_bytes(pickle.dumps(st))
        except OSError:
            pass
    return st[key]


def scores_of(d: Dict[str, pd.Series]) -> dict:
    from experiments import E110_no_kckc as E110
    return E110.scores_of(d)


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale') + tuple(
    f'{k}_{TYPE_COL[t]}' for t in TYPES for k in ('rate_on', 'rate_off', 'ratio'))


def bands_ok(r: Mapping) -> bool:
    from experiments import E110_no_kckc as E110
    return all(E110.bands_of(float(r['rate_on']), float(r['rate_off']), float(r['rate_rho_on']), float(r['opk_on_tv']),
                             float(r['opk_off_tv'])).values())


def collect(kinds: Optional[Sequence[str]] = None) -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run there is: ({(kind, wiring): outputs}, runs table). Kinds: 'family', 'E110_A', 'E110_B' (E110's part A
    and its counted part-B fit), 'S1_<arm>' (step 1), 'R_<arm>' (the refits), 'T_<arm>' (the T reads)."""
    from experiments import E110_no_kckc as E110
    outs: Dict[Tuple[str, str], Dict[str, tuple]] = {}
    rows = []

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **type_rates(o))
        r['in_bands'] = bands_ok(r)
        rows.append(r)

    for w in WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], **F['p'])
        a = E110.load_a(w)
        if a is not None:
            add('E110_A', w, a, tag=f'A_{w}', folder='E110_no_kckc/A')
        b = E110.counted_b(w)
        if b is not None:
            o = E110.load_fit_outputs(b['name'], b['row']['tag'])
            if o is not None:
                add('E110_B', w, o, tag=b['row']['tag'], folder=f"E110_no_kckc/fits/{b['name']}",
                    **E110.params_of_row(b['row']))
    r1 = _rows(out_dir() / 'runs_1.csv')
    for r in r1.to_dict('records'):
        o = load_run('step1', r['tag'])
        if o is not None:
            add(f"S1_{r['arm']}", r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/step1',
                **{k: float(r[k]) for k in PARAM_KEYS})
    for spec in fit_specs(ARMS):
        row = fit_row(spec.name, fit_tag(spec))
        if row is None:
            continue
        o = load_fit_outputs(spec.name, fit_tag(spec))
        if o is not None:
            add(f'R_{spec.arm}', spec.wiring, o, tag=fit_tag(spec), folder=f'{FIT_PREFIX}/{spec.name}', name=spec.name,
                **E110.params_of_row(row), loss=float(row['loss']), n_evals=int(row['n_evals']))
    rt = _rows(out_dir() / 'runs_T.csv')
    for r in rt.to_dict('records'):
        o = load_run('T', r['tag'])
        if o is not None:
            add(f"T_{r['arm']}", r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/T', **{k: float(r[k]) for k in PARAM_KEYS})
    runs = pd.DataFrame(rows)
    if kinds is not None:
        keep = runs['kind'].isin(kinds)
        runs = runs[keep].reset_index(drop=True)
        outs = {k: v for k, v in outs.items() if k[0] in set(kinds)}
    return outs, runs


def readings(outs, runs) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4' / 'T4' when all four of the set have one) and readout: `scores_of`,
    the runs' means of RUN_COLS ('run_<col>'), the number of runs outside the bands."""
    rows = []
    for kind in dict.fromkeys(runs['kind']):
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in (('S4', S_WIRINGS), ('T4', T_WIRINGS))
                                             if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                d = dd[unit] if len(ws) == 1 and unit in dd else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, n_wirings=len(ws), readout=ro,
                                 n_outside_bands=int((~s['in_bands'].astype(bool)).sum()),
                                 **{f'run_{c}': float(s[c].astype(float).mean()) for c in RUN_COLS if c in s.columns},
                                 **scores_of(d)))
    return pd.DataFrame(rows)


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


# ---------------------------------------------------------------------------
# score1: the step-1 gate

def e110_reference() -> dict:
    """E110's readings.csv values on S's population under the readout: the family's and part A's pattern parts."""
    rd = pd.read_csv(H.shared_dir() / 'results' / 'E110_no_kckc' / 'readings.csv')

    def g(kind):
        r = rd[(rd.kind == kind) & (rd.unit == 'S4') & (rd.readout == 'cal')]
        return float(r['all_pattern'].iloc[0]) if len(r) else np.nan
    fam, a = g('family'), g('A')
    return dict(family_S4_pattern_cal=fam, E110_A_S4_pattern_cal=a, half_way=(fam + a) / 2, gate=GATE1)


def score1() -> dict:
    out = out_dir()
    outs, runs = collect(['family', 'E110_A'] + [f'S1_{a}' for a in ARMS])
    runs = runs[runs['wiring'].isin(S_WIRINGS)].reset_index(drop=True)
    outs = {k: v for k, v in outs.items() if k[1] in S_WIRINGS}
    rd = readings(outs, runs)
    rd.to_csv(out / 'step1_readings.csv', index=False)
    runs.to_csv(out / 'step1_runs.csv', index=False)
    ref = e110_reference()
    arms = {}
    for a in ARMS:
        k = f'S1_{a}'
        pat = _v(rd, k, 'S4', 'all_pattern')
        arms[a] = dict(pattern_cal_S4=pat, complete=bool(np.isfinite(pat)), passes=bool(np.isfinite(pat) and pat <= GATE1),
                       d_family=pat - _v(rd, 'family', 'S4', 'all_pattern'),
                       d_E110_A=pat - _v(rd, 'E110_A', 'S4', 'all_pattern'),
                       per_wiring={w: _v(rd, k, w, 'all_pattern') for w in S_WIRINGS},
                       mean_effect_cal=_v(rd, k, 'S4', 'mean_effect'), loss_raw=_v(rd, k, 'S4', 'all_loss', 'raw'),
                       rho_on=_v(rd, k, 'S4', 'rho_on', 'raw'), rho_off=_v(rd, k, 'S4', 'rho_off', 'raw'),
                       n_outside_bands=int(_v(rd, k, 'S4', 'n_outside_bands', 'raw')) if np.isfinite(
                           _v(rd, k, 'S4', 'n_outside_bands', 'raw')) else None)
    g = dict(generated=_now(), host=_host(), gate=GATE1, reference_from_E110=ref,
             recomputed=dict(family=_v(rd, 'family', 'S4', 'all_pattern'), E110_A=_v(rd, 'E110_A', 'S4', 'all_pattern')),
             arms=arms, passing=[a for a in ARMS if arms[a]['passes']], provenance=H.code_commits())
    _write_json(out / 'gate1.json', g)
    print(json.dumps({a: (round(v['pattern_cal_S4'], 4), v['passes']) for a, v in arms.items()}), flush=True)
    return g


# ---------------------------------------------------------------------------
# step 2: the refits

class FitSpec:
    """One refit: name (<arm>_<w>), arm, wiring, start (('params', key) in starts.json; default the step-1 held
    point), maxfev."""

    def __init__(self, name: str, arm: str, wiring: str, start: tuple = None, maxfev: int = MAXFEV):
        self.name, self.arm, self.wiring, self.maxfev = name, arm, wiring, maxfev
        self.start = start if start is not None else ('params', f'S1_{arm}_{wiring}')


def fit_specs(arms: Sequence[str]) -> List[FitSpec]:
    return [FitSpec(f'{a}_{w}', a, w) for a in arms for w in S_WIRINGS]


def fit_tag(spec: FitSpec) -> str:
    """The tag E10_joint_fit.py gives (E110's e10_tag with the kcm-<TAG>_ prefix after kckc<S>_)."""
    from experiments import E110_no_kckc as E110
    r = float(H.data_target('ctrl_frac'))
    suffix = '_cap3_smoke' if spec.name.startswith('smoke') else '_cap3_R'
    return f'{E110.wiring_prefix(spec.wiring)}rate{r:.3f}_aplboth_kcm-E122{spec.arm}_k{K:g}{suffix}'


def fit_command(spec: FitSpec, smoke: bool = False, python: Optional[str] = None) -> List[str]:
    from experiments import E110_no_kckc as E110
    py = str(H.HERE / 'venv' / 'bin' / 'python') if python is None else python
    suffix = '_cap3_smoke' if smoke else '_cap3_R'
    cmd = [py, '-u', 'experiments/E10_joint_fit.py', '-k', f'{K:g}', '--aplpn-max', '3', '--rate-target', 'ctrl_frac',
           '--loss', 'apl', '--apl-readout', 'both', '--tag', suffix, '--maxfev', str(int(spec.maxfev)),
           '--experiment', f'{FIT_PREFIX}/{spec.name}', '--kckc-matrix', str(matrix_path(spec.arm)),
           '--kckc-matrix-tag', f'E122{spec.arm}'] + E110.wiring_args(spec.wiring)
    if spec.start[0] == 'params':
        cmd += ['--init-params', f'results/{EXPERIMENT}/starts.json:{spec.start[1]}']
    if smoke:
        cmd += ['--smoke']
    return cmd


def fit_row(name: str, tag: str) -> Optional[pd.Series]:
    p = H.shared_dir() / 'results' / FIT_PREFIX / name / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    r = f[f['tag'] == tag]
    return r.iloc[-1] if len(r) else None


def load_fit_outputs(name: str, tag: str) -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / FIT_PREFIX / name, archive_dir() / 'fits' / name):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def fit_clean(wait_s: float = 6 * 3600, every: float = 300) -> dict:
    """`E108._commits_clean()` with no tracked file under experiments/ or kcq/ modified (an E10_joint_fit.py process
    loads modules lazily); waits up to `wait_s` for another session's edits to be committed. Untracked files (another
    experiment's new script) are allowed and recorded."""
    t0 = time.time()
    while True:
        r = subprocess.run(['git', '-C', str(H.HERE), 'diff', '--name-only', 'HEAD'], capture_output=True, text=True)
        bad = [p for p in r.stdout.split() if p.startswith(('experiments/', 'kcq/'))]
        if not bad:
            return E108._commits_clean()
        if time.time() - t0 > wait_s:
            raise RuntimeError(f'commit first: {bad} (an E10_joint_fit.py process may load them)')
        print(f'E122: waiting for {bad} to be committed', flush=True)
        time.sleep(every)


def _child_env() -> dict:
    env = dict(os.environ)
    if 'OMP_NUM_THREADS' in env and 'OMP_THREAD_LIMIT' not in env:
        env['OMP_THREAD_LIMIT'] = env['OMP_NUM_THREADS']
    env.setdefault('PYTHONUNBUFFERED', '1')
    return env


def launch(cmd: List[str], res: Path, log: Path, cc: dict, name: str) -> dict:
    """One E10_joint_fit.py process (log `log`) with its launch record res/launch.json."""
    rec = dict(name=name, command=cmd, cwd=str(H.HERE), host=_host(), started=_now(), commits=cc,
               log=str(log.relative_to(H.HERE)), launcher='experiments/E122_feedback_neuron.py',
               threads=dict(OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS'),
                            OMP_THREAD_LIMIT=_child_env().get('OMP_THREAD_LIMIT')))
    res.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    _write_json(res / 'launch.json', rec)
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT, env=_child_env())
    rec.update(finished=_now(), returncode=r.returncode, seconds=time.time() - t0)
    _write_json(res / 'launch.json', rec)
    return rec


def refit_order(g: dict) -> List[str]:
    """The passing arms, main arms first, each group by its step-1 pattern part."""
    p = g['passing']
    key = lambda a: g['arms'][a]['pattern_cal_S4']
    return sorted([a for a in p if a in MAIN_ARMS], key=key) + sorted([a for a in p if a not in MAIN_ARMS], key=key)


def part_2(arms: Sequence[str]) -> None:
    """The refits of `arms` on S, one process at a time (skipping those done), each arm's T read right after its four."""
    for a in arms:
        write_matrix(a)
        for spec in fit_specs([a]):
            if fit_row(spec.name, fit_tag(spec)) is not None:
                print(f'E122 fit {spec.name}: done', flush=True)
                continue
            cc = fit_clean()
            print(f'E122 fit {spec.name}: start {_now()}', flush=True)
            rec = launch(fit_command(spec), H.results_dir(f'{FIT_PREFIX}/{spec.name}'),
                         local_dir('logs') / f'fit_{spec.name}.log', cc, spec.name)
            print(f"E122 fit {spec.name}: return code {rec['returncode']}, {rec['seconds'] / 60:.0f} min", flush=True)
            if rec['returncode'] != 0:
                raise RuntimeError(f'fit {spec.name} failed (log {rec["log"]})')
        part_t(a)


def carried(arm: str) -> dict:
    """The parameters the S refits of an arm carry to T: the mean neuron scale (absolute) and the mean APL->bouton
    factor against the family's fits."""
    from experiments import E110_no_kckc as E110
    sc, ap = [], []
    for spec in fit_specs([arm]):
        r = fit_row(spec.name, fit_tag(spec))
        if r is None:
            raise RuntimeError(f'{spec.name}: no refit yet')
        p = E110.params_of_row(r)
        sc.append(p['wKCKC_scale'])
        ap.append(p['wAPLPN_scale'] / float(E83.saved_fit(spec.wiring)['p']['wAPLPN_scale']))
    return dict(wKCKC_scale=float(np.mean(sc)), aplpn_factor=float(np.mean(ap)), scales=sc, aplpn_factors=ap)


def part_t(arm: str) -> pd.DataFrame:
    """The read on T of an arm (the module docstring), resumable."""
    table = out_dir() / 'runs_T.csv'
    car = carried(arm)
    cc = E108._commits_clean()
    for w in T_WIRINGS:
        tag = f'T_{arm}_{w}'
        if tag in set(_rows(table)['tag']):
            print(f'E122 {tag}: done', flush=True)
            continue
        p = dict(E83.saved_fit(w)['p'])
        p['wKCKC_scale'] = car['wKCKC_scale']
        p['wAPLPN_scale'] = p['wAPLPN_scale'] * car['aplpn_factor']
        hold_point('T', tag, arm, w, p, table, out_dir() / 'provenance_T.json', cc,
                   note=(f"arm {arm}'s matrix at the S refits' mean neuron scale ({car['wKCKC_scale']:.4g}) and APL->bouton "
                         f"factor ({car['aplpn_factor']:.4g} x the family fit's), from the family's fit of {w}; both response "
                         "rates re-held (E77_input_screen.hold_rates)"))
    return _rows(table)


# ---------------------------------------------------------------------------
# score

def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """The class of a refitted arm on T (the module docstring)."""
    k = f'T_{arm}'
    if not all(((runs['kind'] == k) & (runs['wiring'] == w)).any() for w in T_WIRINGS):
        return dict(arm=arm, complete=False)
    d = {}
    for t in TARGETS:
        col = f'{t}_pattern'
        d[t] = dict(pop=_v(rd, k, 'T4', col) - _v(rd, 'family', 'T4', col),
                    per_wiring={w: _v(rd, k, w, col) - _v(rd, 'family', w, col) for w in T_WIRINGS})
    within = all(x['pop'] <= WITHIN and sum(v <= WITHIN for v in x['per_wiring'].values()) >= N_T_NEEDED for x in d.values())
    below_pat = all(x['pop'] <= -BELOW and sum(v <= -BELOW for v in x['per_wiring'].values()) >= N_T_NEEDED
                    for x in d.values())
    d_mean = _v(rd, k, 'T4', 'mean_effect') - _v(rd, 'family', 'T4', 'mean_effect')
    d_raw = _v(rd, k, 'T4', 'all_loss', 'raw') - _v(rd, 'family', 'T4', 'all_loss', 'raw')
    bands = bool(runs[(runs['kind'] == k)]['in_bands'].all())
    below = below_pat and d_mean >= -MEAN_SLACK and d_raw <= RAW_SLACK and bands
    cls = 'below the family' if below else ('replaces KC-KC' if (within and bands) else 'neither')
    return dict(arm=arm, complete=True, cls=cls, pattern_changes=d, within=within, below_pattern=below_pat,
                d_mean_effect_cal=d_mean, d_loss_raw=d_raw, all_T_in_bands=bands, carried=carried(arm))


def comparison_table(rd: pd.DataFrame, kinds: Sequence[str], unit: str) -> str:
    stats = (('pattern part, readout', 'all_pattern', 'cal'), ('amplitude part, readout', 'all_amplitude', 'cal'),
             ('shape part, readout', 'all_shape', 'cal'), ('spread ratio, readout', 'all_sd_ratio', 'cal'),
             ('Pearson r on the per-pair APL effect, readout', 'all_pearson_effect', 'cal'),
             ('mean APL effect, readout', 'mean_effect', 'cal'), ('APL-effect loss, raw', 'all_loss', 'raw'),
             ('mean APL effect, raw', 'mean_effect', 'raw'), ('ρ on odor-odor distances, APL on', 'rho_on', 'raw'),
             ('ρ on odor-odor distances, APL off', 'rho_off', 'raw'), ('response rate, APL on', 'run_rate_on', 'raw'),
             ('APL-off / APL-on response rate', 'run_ratio', 'raw'), ('γ response rate, APL on', 'run_rate_on_g', 'raw'),
             ('γ off / on ratio', 'run_ratio_g', 'raw'), ('KC-KC (neuron) scale', 'run_wKCKC_scale', 'raw'),
             ('somatic threshold', 'run_fixed_thr', 'raw'), ('runs outside the bands', 'n_outside_bands', 'raw'))
    ks = [k for k in kinds if ((rd['kind'] == k) & (rd['unit'] == unit)).any()]
    lines = [f"| statistic ({unit}) | " + ' | '.join(ks) + ' |', '|---|' + '---|' * len(ks)]
    for lab, col, ro in stats:
        lines.append(f'| {lab} | ' + ' | '.join(f'{_v(rd, k, unit, col, ro):.3f}' for k in ks) + ' |')
    return '\n'.join(lines)


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    g = _read_json(out / 'gate1.json')
    refitted = [a for a in ARMS if f'T_{a}' in set(runs['kind'])]
    v = dict(generated=_now(), host=_host(), gate1=g.get('arms'), passing=g.get('passing'),
             classes={a: classify(rd, runs, a) for a in refitted})
    _write_json(out / 'verdict.json', v)
    md = ['# E122 comparison tables (generated)', '',
          comparison_table(rd, ['family', 'E110_A', 'E110_B'] + [f'S1_{a}' for a in ARMS] + [f'R_{a}' for a in ARMS],
                           'S4'), '',
          comparison_table(rd, ['family', 'E110_B'] + [f'T_{a}' for a in ARMS], 'T4'), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script='experiments/E122_feedback_neuron.py',
                                                    commits=H.code_commits()))
    print(json.dumps({a: c.get('cls') for a, c in v['classes'].items()}), flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every run of this experiment, absolute parameters) and model_pair_distances_long.csv (raw spike
    counts) for the model-output archive and kcq/apl_effect.py."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    for r in runs[runs['kind'].str.match(r'^(S1|R|T)_')].itertuples():
        o = outs[(r.kind, r.wiring)]
        arm = r.kind.split('_', 1)[1]
        what = {'S1': "at the fit's scales, both response rates held", 'R': 'refitted on S (threshold, APL scales, neuron scale)',
                'T': "the read on T: the S refits' neuron scale and APL->bouton factor, both response rates held"}[r.kind.split('_')[0]]
        subdir = r.folder.split('/', 1)[1] if '/' in r.folder else r.folder
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.kind.split('_')[0], arm=arm, arm_label=ARM_LABELS[arm], base='bouton',
                         k=K, wiring=r.wiring, kind=f"C2's family's {r.wiring}: KC-KC matrix of arm {arm} ({ARM_LABELS[arm]}), {what}",
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         kckc_matrix=f'matrices/kckc_{arm}.parquet', fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale, rate_on=r.rate_on, rate_off=r.rate_off,
                         max_count_on=r.max_count_on, max_count_off=r.max_count_off, realism_flag=r.realism_flag,
                         subdir=subdir))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def chain() -> None:
    """check, step 1, score1, the gated refits with their T reads, score, variants."""
    part_check()
    part_1()
    g = score1()
    order = refit_order(g)
    print(f'E122 step 1 gate: passing {order}', flush=True)
    if order:
        part_2(order)
    score()
    write_variants()
    print(f'E122 chain done {_now()}', flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'check', '1', 'score1', '2', 'T', 'score', 'variants', 'chain',
                                                      'commands'))
    ap.add_argument('--arms', nargs='+', default=None, help='arms for parts 1, 2 and T (default: all / the passing ones)')
    a = ap.parse_args()
    if a.part == '0':
        part_0()
    elif a.part == 'check':
        part_check()
    elif a.part == '1':
        part_1(a.arms or ARMS)
    elif a.part == 'score1':
        score1()
    elif a.part == '2':
        part_2(a.arms or refit_order(_read_json(out_dir() / 'gate1.json')))
    elif a.part == 'T':
        for arm in a.arms or []:
            part_t(arm)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'chain':
        chain()
    elif a.part == 'commands':
        for s in fit_specs(a.arms or ARMS):
            print(' '.join(fit_command(s)))


if __name__ == '__main__':
    main()
