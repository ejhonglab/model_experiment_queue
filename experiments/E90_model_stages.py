#!/usr/bin/env python3
"""E90: the model's stages against the flies' stages. The model's PN, bouton and claw stages read as the PN-dendrite,
bouton and claw recordings were made; whether the ester and the alcohol pairs leave the flies before the claws, at the
step from claws to somata with APL on, or at that step without APL; and the step in wild-type flies as a second reading
of the main target's pattern that uses no fly with APL output blocked (TNT).

The entry is item E90 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from atlas
2026-11-13 11:00). C2's family throughout (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72
part B's eight wirings at its own fit (E83's `saved_fit` / `build_simulator`); a population is the mean over wirings of
the pair distances. The 66 pairs of the panel's 12 megamat odors, in groups (GROUPS): the ester group (15 pairs), the
alcohol group (6), ester with alcohol (24), the other pairs (21), and E78's six near pairs.

Parts (CLI `--part`):
- A (no simulations): A1 the flies' stage table (stages_flies.csv; per pair and group, 90 % intervals from resampling
  each stage's flies, N_BOOT resamples): ORNs (E85's per-fly reader, olfdb dataset 1, the model's 36 glomeruli) and PN
  dendrites (GH146, al_analysis's sent_to_remy/2023-10-29 tables, beside the pebbled ORNs of the same table) across
  glomeruli, per fly then averaged; Remy's PN boutons and KC claws from the per-frame matrices over STAGE_WINDOW, as
  measured (the z-scored `stim` matrices, as data/remy_rdm_timeseries/windows/) and with the noise removed through the
  block averages (each fly's window-mean block-average correlations divided by the geometric mean of the two odors'
  reliabilities, the diagonal); the single-trial correlations of rdm_timeseries_explainers.html beside where they decode;
  the somata (Remy's wild type, Yang's control flies, Yang's flies with APL output blocked (TNT); per fly / recording,
  1 - Pearson and 1 - cosine). A2 the two steps per pair (steps_flies.csv): somata with APL on minus claws (Remy's
  somata; Yang's control flies) and somata with APL output blocked (TNT) minus claws, against Yang's per-pair APL effect
  and the error common to all wirings (E72's common_error_pairs.csv), the step's floor (flies resampled at both stages)
  and the agreement of the claw flies' halves (the three splits of two and two); step_target.json. A3 the antennal lobe
  as recorded (al_transfer.csv): PN dendrites against ORNs on the shared glomeruli (the transfer, each odor's summed
  response, the pair distances, per fly), also under a saturating map from dF/F to rate (the bound for the indicator).
- B0 (simulations): identity: each wiring's run with the rates saved reproduces its saved fit's counts exactly, APL on
  and off; on the real wiring the claw inputs formed as wiring x bouton rate minus APL's inhibition, clipped at 0,
  reproduce a run's `claw_sims` exactly (identity.json).
- B (simulations, resumable per wiring): each wiring's PN rates (pn_sims), bouton rates after APL's inhibition
  (bouton_sims) and APL's activity (inh_sims), APL on and off, read as window means over MODEL_WINDOW minus the
  pre-odor means; the stage matrices (glomeruli / boutons / claws x odors) saved under outputs/E90_model_stages/B/.
  B5: the real wiring at the flies' 2 s pulse at E82's re-held B1 point (its saved counts reproduced first), the stages
  read over PULSE_WINDOW.
- Bscore (no simulations): B1 the PN stage against A3 (entries, per-odor totals, pair distances; pn_stage.csv); B2 the
  bouton stage and B3 the claw stage (the drive alone and with CLAW_SPIKE_SHARES of the KC's spike signal) in the three
  READINGS ('rate', 'session_z', 'noisy': session_z through measurement noise at the flies' block-to-block
  reliability, the noise then removed as in A1); B4 the soma stage (raw counts and the calibrated readout, Pearson and
  cosine); the check of the noise removal (noise_check.csv); stages_model.csv, steps_model.csv.
- C (no simulations): the verdict per group, fixed before any run (VERDICT_* below; verdict.json) and the step target's
  usability (step_target.json).
- D1 / D2 (only as the verdict says): D1 the firing-ceiling family (E69's fit 1, E72 part G's fits on P1, P2, M1) read
  the same way beside C2's; D2 the PN stage of E77's screen's normalization (S1) and exponent (S3) settings beside their
  saved main-target scores (d2_pn_settings.csv).
- table: comparison_table.md and the figures.

Outputs: results/E90_model_stages/ (stages_flies.csv, steps_flies.csv, step_target.json, al_transfer.csv,
stages_model.csv, steps_model.csv, noise_check.csv, verdict.json, comparison_table.md, figures, provenance_<part>.json);
stage matrices and runs under outputs/E90_model_stages/ (archived to data/model_outputs/E90_model_stages/).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E90_model_stages.py --part A
    OMP_NUM_THREADS=12 OMP_THREAD_LIMIT=12 nice -n 10 ./venv/bin/python -u experiments/E90_model_stages.py --part B0
    OMP_NUM_THREADS=12 OMP_THREAD_LIMIT=12 nice -n 10 ./venv/bin/python -u experiments/E90_model_stages.py --part B
    ./venv/bin/python -u experiments/E90_model_stages.py --part Bscore
    ./venv/bin/python -u experiments/E90_model_stages.py --part C
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E90_model_stages'

# ---------------------------------------------------------------------------
# the panel's megamat odors, the pairs and the groups

ODORS = tuple(H.MEGAMAT12)
ESTER = ('2h', 'IaA', 'pa', 'eb', 'ep', 'ms')
ALCOHOL = ('1-8ol', '1-5ol', '1-6ol', '6al')
NEAR = (('eb', 'ep'), ('6al', 't2h'), ('IaA', 'pa'), ('2h', 'pa'), ('1-5ol', '1-6ol'), ('2h', 'IaA'))   # E78's SIX
GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near', 'all66')
GROUP_NAMES = {'ester': 'ester group (15 pairs)', 'alcohol': 'alcohol group (6)', 'ester_alcohol': 'ester with alcohol (24)',
               'other': 'the other pairs (21)', 'near': "E78's near pairs (6)", 'all66': 'all 66 megamat pairs'}

# ---------------------------------------------------------------------------
# part A: the flies

STAGE_WINDOW = (0.5, 3.0)               # Remy's boutons and claws: per-frame matrices averaged over this window (s)
N_BOOT = 1000
CI = (5.0, 95.0)                        # 90 % intervals
BOOT_SEED = 0
REMY_VARIANT = {'bouton': 'pn_{kind}_F_zscore', 'claw': 'kc_{kind}_Fc_zscore'}   # the matrices Remy's figures use
SOMA_SOURCES = ('remy', 'ctrl', 'tnt')  # Remy's wild type; Yang's control flies; Yang's APL output blocked (TNT)
# A3's saturating map from dF/F to rate (the bound for the indicator): r = x / (1 - x / (SAT_HEADROOM * x_max)), x_max the
# stage's largest fly-mean entry
SAT_HEADROOM = 1.2

# ---------------------------------------------------------------------------
# part B: the model

MODEL_WINDOW = (0.0, 0.5)               # the 0.5 s pulse: the stages' window (s)
PRE_WINDOW = (-0.5, 0.0)                # the pre-odor window whose mean is subtracted (s)
PULSE = 2.0                             # B5: the flies' pulse (E82)
PULSE_WINDOW = (0.5, 2.0)               # B5: the stages read over this window (s)
READINGS = ('rate', 'session_z', 'noisy')
CLAW_SPIKE_SHARES = (0.0, 0.25, 0.5)    # B3's versions of what a claw's calcium reports
N_NOISE_REPS = 20                       # replicate noise draws per reading (each of the stage's number of flies)
NOISE_SEED = 0
CLAW_IDENTITY_ODORS = ('eb', '1-8ol')   # B0: the odors of the run with claw_sims saved (memory)

# ---------------------------------------------------------------------------
# parts C and D: the verdict, fixed before any run

VERDICT_DIFF = 0.05                     # the population's group mean this far or more from the flies'
VERDICT_SAME_SIDE = 6                   # of the 8 wirings on the same side
VERDICT_EFFECT_ERR = 0.03               # "not resolved" when nothing differs and the group's error in the effect is above
NOISE_CHECK_TOL = 0.03                  # the noise removal within this of the truth for every group (E88's A0)
STEP_HALVES_MIN = 0.6                   # the step target: the claw flies' halves agree at this or more on average
STEP_FLOOR_RATIO = 0.5                  # and its floor at most this times its SD over the 66 pairs
VERDICTS = ('upstream', 'the step with APL on', 'the step without APL', 'both steps', 'not resolved', 'none')


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# pairs and groups

def _canon(a: str, b: str) -> Tuple[str, str]:
    return tuple(sorted((str(a), str(b))))


def _cs(s: pd.Series) -> pd.Series:
    """A pair Series in canonical order (each pair's odors sorted, then the pairs; E83's `_cs`)."""
    s = s.copy()
    s.index = pd.MultiIndex.from_tuples([_canon(*x) for x in s.index], names=['odor_a', 'odor_b'])
    return s.sort_index()


def megamat_pairs() -> pd.MultiIndex:
    """The 66 megamat pairs in canonical (sorted) order, as E83's `_cs`."""
    return pd.MultiIndex.from_tuples(sorted({_canon(a, b) for a, b in combinations(ODORS, 2)}), names=['odor_a', 'odor_b'])


def group_masks(pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """Boolean masks over `pairs` for every name in GROUPS ('all66': both odors megamat)."""
    a = pd.Index(pairs.get_level_values(0))
    b = pd.Index(pairs.get_level_values(1))
    mm = np.asarray(a.isin(ODORS) & b.isin(ODORS))
    ester = np.asarray(a.isin(ESTER) & b.isin(ESTER))
    alcohol = np.asarray(a.isin(ALCOHOL) & b.isin(ALCOHOL))
    ea = np.asarray((a.isin(ESTER) & b.isin(ALCOHOL)) | (a.isin(ALCOHOL) & b.isin(ESTER)))
    near = {_canon(*q) for q in NEAR}
    return dict(ester=ester, alcohol=alcohol, ester_alcohol=ea, other=mm & ~(ester | alcohol | ea),
                near=np.array([_canon(x, y) in near for x, y in pairs]), all66=mm)


def group_means(pairs: pd.Series) -> Dict[str, float]:
    """A pair Series' mean over each group's pairs (NaN pairs left out)."""
    s = _cs(pairs)
    m = group_masks(s.index)
    return {g: float(np.nanmean(s.to_numpy()[m[g]])) if m[g].any() else np.nan for g in GROUPS}


def _pairs_from_square(mat: np.ndarray, odors: Sequence[str]) -> pd.Series:
    """The 66 megamat pairs of an odor x odor matrix (labels `odors`), canonical order."""
    ix = {o: i for i, o in enumerate(odors)}
    idx = megamat_pairs()
    return pd.Series([float(mat[ix[a], ix[b]]) for a, b in idx], index=idx)


def _rows_to_pairs(rows: Dict[str, pd.Series]) -> pd.DataFrame:
    """{unit: pair Series} -> units x the 66 pairs (canonical order)."""
    idx = megamat_pairs()
    return pd.DataFrame({k: _cs(v).reindex(idx) for k, v in rows.items()}).T


def _pearson_dist_matrix(X: np.ndarray) -> np.ndarray:
    """1 - Pearson r between the rows (odors) of X (odors x units)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        return 1.0 - np.corrcoef(X)


def _cosine_dist_matrix(X: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(X, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return 1.0 - (X @ X.T) / np.outer(n, n)


def odor_distances(X: pd.DataFrame, measure: str = 'pearson') -> pd.Series:
    """The 66 megamat pair distances of a units x odors table (columns odor abbreviations), across its units."""
    M = X[list(ODORS)].to_numpy(dtype=float).T
    D = _pearson_dist_matrix(M) if measure == 'pearson' else _cosine_dist_matrix(M)
    return _pairs_from_square(D, ODORS)


# ---------------------------------------------------------------------------
# A: readers

def remy_dir() -> Path:
    return H.shared_dir() / 'data' / 'remy_rdm_timeseries'


def remy_matrices(stage: str, kind: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Remy's per-fly per-frame matrices of `stage` ('bouton' / 'claw'), `kind` 'stim' or 'blockavg' (REMY_VARIANT):
    (flies x frames x 17 x 17 distances 1 - r, frame times, odors)."""
    d = remy_dir()
    a = np.load(d / (REMY_VARIANT[stage].format(kind=kind) + '.npy'))
    t = np.load(d / ('times_pn.npy' if stage == 'bouton' else 'times_kc.npy'))
    odors = json.loads((d / 'odors.json').read_text())
    if a.shape[1] != len(t) or a.shape[2] != len(odors):
        raise RuntimeError(f'{stage} {kind}: shape {a.shape}, {len(t)} frames, {len(odors)} odors')
    return a, t, odors


def _frames(t: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    """The frames of a closed window (as window_csvs.py selects its closed windows)."""
    return (t >= window[0] - 1e-6) & (t <= window[1] + 1e-6)


def remy_window_means(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> np.ndarray:
    """Per fly, the mean over the window's frames of the `stim` matrices (flies x 17 x 17, 1 - r), as
    decode_scripts/window_csvs.py makes data/remy_rdm_timeseries/windows/ (closed window)."""
    a, t, _ = remy_matrices(stage, 'stim')
    return a[:, _frames(t, window)].mean(axis=1)


def remy_reliability(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> pd.DataFrame:
    """Per fly and odor, the window-mean block-average reliability (1 - the diagonal of the 1 - r matrices); flies x 17."""
    a, t, odors = remy_matrices(stage, 'blockavg')
    r = 1.0 - a[:, _frames(t, window)].mean(axis=1)
    return pd.DataFrame(np.diagonal(r, axis1=1, axis2=2), columns=odors,
                        index=pd.Index([f'fly{i + 1}' for i in range(a.shape[0])], name='fly'))


def remy_noise_removed(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> np.ndarray:
    """Per fly, 1 - (window-mean block-average correlation of two odors) / sqrt(product of the two odors' window-mean
    reliabilities, the diagonal); flies x 17 x 17 (NaN where a reliability is not positive)."""
    a, t, _ = remy_matrices(stage, 'blockavg')
    r = 1.0 - a[:, _frames(t, window)].mean(axis=1)
    rel = np.diagonal(r, axis1=1, axis2=2)
    den = rel[:, :, None] * rel[:, None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, 1.0 - r / np.sqrt(np.where(den > 0, den, 1.0)), np.nan)


_SINGLE: Dict[str, object] = {}


def remy_single_trial(stage: str, window: Tuple[float, float] = STAGE_WINDOW) -> Optional[dict]:
    """The single-trial correlations of rdm_timeseries_explainers.html (its data script's `trial`: per fly and frame
    the correlation of every pair of the 51 trials), over the window's frames: 'between' (flies x 17 x 17, the mean over
    the 3 x 3 pairs of trials of two odors; the diagonal the mean over the 3 pairs of different trials of the odor, the
    reliability), 'odors', 'n_frames'; None if the HTML or its data do not decode."""
    import re
    if 'D' not in _SINGLE:
        try:
            t = (H.shared_dir() / 'from_remy' / 'rdm_timeseries_explainers.html').read_text(encoding='utf8')
            m = re.search(r'<script[^>]*id="data"[^>]*>', t)
            _SINGLE['D'] = json.loads(t[m.end():t.index('</script>', m.end())])
        except Exception as e:                  # (the HTML is optional: its numbers are given beside only)
            print(f'single-trial correlations not decoded: {e!r}', flush=True)
            _SINGLE['D'] = None
    D = _SINGLE['D']
    if D is None:
        return None
    tr = D['trial']['pn' if stage == 'bouton' else 'kc']
    r = np.asarray(tr['r'], dtype=float)                       # flies x frames x trial pairs
    tt = np.asarray(tr['t'], dtype=float)
    sel = _frames(tt, window)
    rm = r[:, sel].mean(axis=1)                                # flies x trial pairs
    oo = np.asarray(tr['odor_of_trial'], dtype=int)
    n = len(D['odors'])
    sums = np.zeros((r.shape[0], n, n))
    cnt = np.zeros((n, n))
    for k, (i, j) in enumerate(tr['pairs']):
        a, b = oo[i], oo[j]
        sums[:, a, b] += rm[:, k]
        cnt[a, b] += 1
        if a != b:
            sums[:, b, a] += rm[:, k]
            cnt[b, a] += 1
    return dict(between=sums / cnt[None], odors=list(D['odors']), n_frames=int(sel.sum()))


def gh146_dir() -> Path:
    return H.al_analysis_data_dir() / 'sent_to_remy' / '2023-10-29'


def al_table(kind: str) -> pd.DataFrame:
    """The megamat trial-mean responses (dF/F) of al_analysis's sent_to_remy/2023-10-29 table `kind` ('GH146' /
    'pebbled'): rows (fly 'date/fly_num', odor abbreviation), columns glomeruli (NaN where a fly did not have one)."""
    df = pd.read_pickle(gh146_dir() / f'{kind}_ij_certain-roi_stats.p')
    m = df.loc['megamat']
    m = m[m.index.get_level_values('odor2') == 'solvent']
    m = m.groupby(level='odor1').mean()                        # mean over the three repeats
    m.index = [H.odor_abbrev(o) for o in m.index]
    rows = {}
    for (date, fly), sub in m.T.groupby(level=['date', 'fly_num']):
        x = sub.droplevel(['date', 'fly_num']).T               # odors x glomeruli
        x.columns = list(x.columns)
        rows[f'{pd.Timestamp(date).date()}/{int(fly)}'] = x
    out = pd.concat(rows, names=['fly', 'odor'])
    missing = [o for o in ODORS if o not in out.index.get_level_values('odor')]
    if missing:
        raise RuntimeError(f'{kind}: megamat odors missing {missing}')
    return out


def shared_glomeruli() -> List[str]:
    """The 28 glomeruli the GH146 recordings and the model's input share (GH146, pebbled and the model's 36 input
    glomeruli; the note of 2026-11-07)."""
    gh = set(al_table('GH146').columns)
    pb = set(al_table('pebbled').columns)
    model = set(H.apl16_orn_deltas().index)
    return sorted(gh & pb & model)


def per_fly_distances(table: pd.DataFrame, glomeruli: Optional[Sequence[str]] = None,
                      measure: str = 'pearson') -> pd.DataFrame:
    """Per fly the 66 pair distances across the fly's glomeruli (those in `glomeruli` it has for every megamat odor);
    flies x pairs."""
    rows = {}
    for fly, x in table.groupby(level='fly', sort=False):
        x = x.droplevel('fly').reindex(list(ODORS))
        cols = [g for g in x.columns if (glomeruli is None or g in glomeruli) and x[g].notna().all()]
        rows[fly] = odor_distances(x[cols].T, measure)
    return _rows_to_pairs(rows)


def orn_fly_distances() -> pd.DataFrame:
    """The ORN flies' (E85's reader, olfdb dataset 1) per-fly 66 pair distances, 1 - Pearson of the trial means across
    the fly's measured glomeruli of the model's 36; flies x pairs."""
    from experiments import E85_orn_fly_population as E85
    F = E85.flies()
    rows = {}
    with np.errstate(invalid='ignore'):
        m = np.nanmean(F['x'], axis=3)                         # flies x odors x glomeruli
    for f, lab in enumerate(F['labels']):
        x = pd.DataFrame(m[f], index=F['odors'], columns=F['gloms'])
        x = x.loc[:, x.notna().all(axis=0)]
        rows[lab] = odor_distances(x.T, 'pearson')
    return _rows_to_pairs(rows)


def remy_soma_fly_distances() -> pd.DataFrame:
    """Remy's wild-type megamat KC somata, per fly (18; the per-fly correlation matrices she gave Anoop,
    pyolfsysm's data/from_remy/2024-11-12/xrda_stim_rdm_concat__iscell_good_xid0__correlation.nc), the 66 pair
    correlation distances 1 - Pearson r as stored (NaN where a fly lacks an odor); flies x pairs. Remy's somata exist as
    correlation distances only, so there is no cosine reading of them."""
    import xarray as xr
    from pyolfsysm.mb_model import data_root
    da = xr.load_dataarray(Path(data_root) / 'from_remy' / '2024-11-12' /
                           'xrda_stim_rdm_concat__iscell_good_xid0__correlation.nc')
    odors = [str(o) for o in da.coords['abbrev_row'].values]
    if odors != [str(o) for o in da.coords['abbrev_col'].values]:
        raise RuntimeError('row and column odors differ')
    names = [str(x) for x in da.coords['title_str'].values]
    rows = {n: _pairs_from_square(da.values[i], odors) for i, n in enumerate(names)}
    return _rows_to_pairs(rows)


def yang_soma_distances(genotype: str, measure: str = 'pearson') -> pd.DataFrame:
    """Yang's per-recording 66 megamat pair distances (`H.yang_recording_rdms`); recordings x pairs."""
    r = H.yang_recording_rdms(genotype, measure)
    rows = {'/'.join(map(str, k)): _cs(r.loc[k]) for k in r.index}
    return _rows_to_pairs(rows)


def _flies_pairs(mats: np.ndarray, odors: Sequence[str], prefix: str = 'fly') -> pd.DataFrame:
    return _rows_to_pairs({f'{prefix}{i + 1}': _pairs_from_square(mats[i], odors) for i in range(mats.shape[0])})


_SAMPLES: Dict[str, pd.DataFrame] = {}

# every fly stage: key -> (stage, reading, measure, description)
STAGE_KEYS = {
    'orn': ('orn', 'olfdb', 'pearson', "ORNs across the model's 36 glomeruli (olfdb dataset 1, 9 flies)"),
    'orn_pebbled28': ('orn', 'pebbled_shared28', 'pearson', 'ORNs across the 28 shared glomeruli (pebbled table)'),
    'pn_dend': ('pn_dend', 'gh146', 'pearson', 'PN dendrites across glomeruli (GH146, 7 flies, 29 glomeruli)'),
    'pn_dend28': ('pn_dend', 'gh146_shared28', 'pearson', 'PN dendrites across the 28 shared glomeruli'),
    'bouton_measured': ('bouton', 'measured', 'pearson', "PN boutons as measured (Remy's z-scored stim matrices)"),
    'bouton_noise_removed': ('bouton', 'noise_removed', 'pearson', 'PN boutons, noise removed through the block averages'),
    'claw_measured': ('claw', 'measured', 'pearson', "KC claws as measured (Remy's z-scored stim matrices)"),
    'claw_noise_removed': ('claw', 'noise_removed', 'pearson', 'KC claws, noise removed through the block averages'),
    'soma_remy': ('soma', 'remy', 'pearson', "KC somata, Remy's wild type (18 flies)"),
    'soma_ctrl': ('soma', 'ctrl', 'pearson', "KC somata, control flies (Yang, 11 recordings)"),
    'soma_ctrl_cos': ('soma', 'ctrl', 'cosine', "KC somata, control flies (Yang), cosine"),
    'soma_tnt': ('soma', 'tnt', 'pearson', 'KC somata, APL output blocked (TNT) (Yang, 8 recordings)'),
    'soma_tnt_cos': ('soma', 'tnt', 'cosine', 'KC somata, APL output blocked (TNT) (Yang), cosine'),
}


def stage_samples() -> Dict[str, pd.DataFrame]:
    """Every fly stage of STAGE_KEYS as units of resampling (flies; Yang's recordings) x the 66 pairs; cached."""
    if _SAMPLES:
        return _SAMPLES
    gh, pb = al_table('GH146'), al_table('pebbled')
    sh = shared_glomeruli()
    s = dict(orn=orn_fly_distances(), orn_pebbled28=per_fly_distances(pb, sh),
             pn_dend=per_fly_distances(gh), pn_dend28=per_fly_distances(gh, sh))
    for st in ('bouton', 'claw'):
        _, _, odors = remy_matrices(st, 'stim')
        s[f'{st}_measured'] = _flies_pairs(remy_window_means(st), odors)
        s[f'{st}_noise_removed'] = _flies_pairs(remy_noise_removed(st), odors)
    s['soma_remy'] = remy_soma_fly_distances()
    for g, k in (('TNT3c', 'ctrl'), ('TNT', 'tnt')):
        s[f'soma_{k}'] = yang_soma_distances(g, 'pearson')
        s[f'soma_{k}_cos'] = yang_soma_distances(g, 'cosine')
    _SAMPLES.update(s)
    return _SAMPLES


def bootstrap_mean(samples: pd.DataFrame, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> np.ndarray:
    """n_boot x pairs means (NaN left out) over resampled rows (flies or recordings, with replacement)."""
    X = samples.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, X.shape[0], size=(n_boot, X.shape[0]))
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.stack([np.nanmean(X[i], axis=0) for i in idx])


def _group_boot(boot: np.ndarray, pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    m = group_masks(pairs)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return {g: np.nanmean(boot[:, m[g]], axis=1) for g in GROUPS}


def _summary_rows(key_cols: dict, mean: pd.Series, boot: np.ndarray, n: int) -> List[dict]:
    """Per pair and per group: the mean, the bootstrap SD and the 90 % interval."""
    rows = []
    lo, hi = np.nanpercentile(boot, CI, axis=0)
    sd = np.nanstd(boot, axis=0)
    for i, (a, b) in enumerate(mean.index):
        rows.append(dict(key_cols, n=n, kind='pair', odor_a=a, odor_b=b, group='', mean=float(mean.iloc[i]),
                         lo=float(lo[i]), hi=float(hi[i]), sd_boot=float(sd[i])))
    gm = group_means(mean)
    gb = _group_boot(boot, mean.index)
    for g in GROUPS:
        glo, ghi = np.nanpercentile(gb[g], CI)
        rows.append(dict(key_cols, n=n, kind='group', odor_a='', odor_b='', group=g, mean=gm[g], lo=float(glo),
                         hi=float(ghi), sd_boot=float(np.nanstd(gb[g]))))
    return rows


# ---------------------------------------------------------------------------
# A: tables

def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def part_a1() -> pd.DataFrame:
    """stages_flies.csv: per stage key (STAGE_KEYS), pair and group, the mean over the stage's flies (Yang's
    recordings) and its 90 % interval and SD from N_BOOT resamples of them."""
    S = stage_samples()
    rows = []
    for k, X in S.items():
        stage, reading, measure, desc = STAGE_KEYS[k]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            mean = pd.Series(np.nanmean(X.to_numpy(dtype=float), axis=0), index=X.columns)
        rows += _summary_rows(dict(key=k, stage=stage, reading=reading, measure=measure, description=desc), mean,
                              bootstrap_mean(X), X.shape[0])
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'stages_flies.csv', index=False)
    return df


def flies_stage_mean(key: str) -> pd.Series:
    """A fly stage's mean pair distances (the 66 pairs, canonical order)."""
    X = stage_samples()[key]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return pd.Series(np.nanmean(X.to_numpy(dtype=float), axis=0), index=X.columns)


# the flies' steps: name -> (soma key, claw key's stage); the claw readings STEP_CLAW_READINGS
FLY_STEPS = {'on_remy': 'soma_remy', 'on_ctrl': 'soma_ctrl', 'off_tnt': 'soma_tnt',
             'on_ctrl_cos': 'soma_ctrl_cos', 'off_tnt_cos': 'soma_tnt_cos'}
STEP_CLAW_READINGS = ('measured', 'noise_removed')
# the step target's rule is applied to Remy's wild-type somata, in both of these claw readings
STEP_TARGET = 'on_remy'
STEP_TARGET_READINGS = ('measured', 'noise_removed')
SOMA_HALF_SEED = 0


def _step_boot(soma: pd.DataFrame, claw: pd.DataFrame, seed: int = BOOT_SEED) -> np.ndarray:
    """N_BOOT x pairs: the step (soma minus claw) with the flies of both stages resampled independently."""
    return bootstrap_mean(soma, seed=seed) - bootstrap_mean(claw, seed=seed + 1)


def _corr(a, b, method: str = 'pearson') -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    if method == 'spearman':
        from scipy.stats import spearmanr
        return float(spearmanr(a[ok], b[ok])[0])
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def claw_half_splits(n: int = 4) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """The three splits of four claw flies into two and two (each split once)."""
    out = []
    for h in combinations(range(n), n // 2):
        rest = tuple(i for i in range(n) if i not in h)
        if 0 in h:
            out.append((h, rest))
    return out


def step_halves(soma: pd.DataFrame, claw: pd.DataFrame) -> dict:
    """The agreement (Pearson r over the 66 pairs) of the step computed from each half of the claw flies, for the three
    splits: 'common' with all the somata's flies in both halves' steps (the claws' halves only), 'split' with the somata's
    flies also split in two at random (seed SOMA_HALF_SEED), so that the halves share no fly."""
    s_all = soma.mean(axis=0, skipna=True)
    rng = np.random.default_rng(SOMA_HALF_SEED)
    perm = rng.permutation(soma.shape[0])
    s_h = (soma.iloc[perm[: soma.shape[0] // 2]].mean(axis=0), soma.iloc[perm[soma.shape[0] // 2:]].mean(axis=0))
    out = dict(common=[], split=[], common_spearman=[])
    for h1, h2 in claw_half_splits(claw.shape[0]):
        c1, c2 = claw.iloc[list(h1)].mean(axis=0), claw.iloc[list(h2)].mean(axis=0)
        out['common'].append(_corr(s_all - c1, s_all - c2))
        out['common_spearman'].append(_corr(s_all - c1, s_all - c2, 'spearman'))
        out['split'].append(_corr(s_h[0] - c1, s_h[1] - c2))
    return {k: dict(per_split=v, mean=float(np.mean(v))) for k, v in out.items()}


def _target_effect() -> pd.Series:
    """Yang's per-pair APL effect (the main target: cosine, all recordings), the 66 megamat pairs."""
    return _cs(H.yang_apl_effect()).reindex(megamat_pairs())


def _common_error() -> pd.Series:
    """E72's error common to all wirings (calibrated readout), the 66 megamat pairs."""
    d = pd.read_csv(H.results_dir('E72_wiring_population') / 'common_error_pairs.csv')
    return _cs(d.set_index(['odor_a', 'odor_b'])['common']).reindex(megamat_pairs())


def part_a2() -> pd.DataFrame:
    """steps_flies.csv (per step, claw reading, pair and group: the step's mean, interval and SD; and per step and claw
    reading its correlations with Yang's per-pair APL effect and with the common error, its floor, its SD over the 66
    pairs and the agreement of the claw flies' halves) and step_target.json."""
    S = stage_samples()
    eff, ce = _target_effect(), _common_error()
    rows, summ = [], []
    for step, sk in FLY_STEPS.items():
        for cr in STEP_CLAW_READINGS:
            ck = f'claw_{cr}'
            if ck not in S:
                continue
            soma, claw = S[sk], S[ck]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                st = pd.Series(np.nanmean(soma.to_numpy(dtype=float), axis=0) - np.nanmean(claw.to_numpy(dtype=float), axis=0),
                               index=soma.columns)
            boot = _step_boot(soma, claw)
            rows += _summary_rows(dict(step=step, soma=sk, claw=ck), st, boot, int(claw.shape[0]))
            dev = boot - st.to_numpy()[None]
            dev_c = dev - np.nanmean(dev, axis=1, keepdims=True)
            halves = step_halves(soma, claw)
            gb = _group_boot(boot, st.index)
            d_ea = gb['ester'] - gb['alcohol']
            gm = group_means(st)
            summ.append(dict(step=step, soma=sk, claw=ck, n_soma=int(soma.shape[0]), n_claw=int(claw.shape[0]),
                             sd_over_pairs=float(np.std(st.to_numpy())),
                             floor=float(np.mean(np.sqrt(np.nanmean(dev_c ** 2, axis=1)))),
                             floor_uncentred=float(np.mean(np.sqrt(np.nanmean(dev ** 2, axis=1)))),
                             median_pair_sd=float(np.nanmedian(np.nanstd(boot, axis=0))),
                             halves_common=halves['common']['mean'], halves_common_spearman=halves['common_spearman']['mean'],
                             halves_split=halves['split']['mean'],
                             halves_common_per_split=json.dumps([round(x, 4) for x in halves['common']['per_split']]),
                             halves_split_per_split=json.dumps([round(x, 4) for x in halves['split']['per_split']]),
                             r_effect=_corr(st, eff), rho_effect=_corr(st, eff, 'spearman'),
                             r_common_error=_corr(st, ce), rho_common_error=_corr(st, ce, 'spearman'),
                             ester_minus_alcohol=gm['ester'] - gm['alcohol'],
                             ester_minus_alcohol_lo=float(np.nanpercentile(d_ea, CI[0])),
                             ester_minus_alcohol_hi=float(np.nanpercentile(d_ea, CI[1])),
                             **{f'mean_{g}': v for g, v in group_means(st).items()}))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'steps_flies.csv', index=False)
    sm = pd.DataFrame(summ)
    sm.to_csv(out_dir() / 'steps_flies_summary.csv', index=False)
    target = {}
    for cr in STEP_TARGET_READINGS:
        r = sm[(sm.step == STEP_TARGET) & (sm.claw == f'claw_{cr}')].iloc[0]
        target[cr] = dict(halves_agreement=float(r.halves_common), halves_agreement_split_somata=float(r.halves_split),
                          floor=float(r.floor), sd_over_pairs=float(r.sd_over_pairs),
                          halves_ok=bool(r.halves_common >= STEP_HALVES_MIN),
                          floor_ok=bool(r.floor <= STEP_FLOOR_RATIO * r.sd_over_pairs))
        target[cr]['usable'] = target[cr]['halves_ok'] and target[cr]['floor_ok']
    beside = {}
    for _, r in sm[sm.step != STEP_TARGET].iterrows():
        beside[f'{r.step}__{r.claw}'] = dict(halves_agreement=float(r.halves_common), floor=float(r.floor),
                                             sd_over_pairs=float(r.sd_over_pairs),
                                             usable=bool(r.halves_common >= STEP_HALVES_MIN and r.floor <= STEP_FLOOR_RATIO * r.sd_over_pairs))
    _write_json(out_dir() / 'step_target.json', dict(
        rule=f'usable if the claw flies\' halves agree on the step at {STEP_HALVES_MIN} or more on average (Pearson r over '
             f'the 66 pairs, the three splits of two and two, all of the somata\'s flies in both halves) and its floor (the '
             f'root mean squared difference over the 66 pairs, means removed, between a resampled step and the measured '
             f'one, flies of both stages resampled; mean over {N_BOOT} resamples) is at most {STEP_FLOOR_RATIO} times its '
             f'SD over the 66 pairs; applied to the step from the claws to Remy\'s wild-type somata in both claw readings '
             f'{list(STEP_TARGET_READINGS)} (fixed before the flies were read)',
        step=STEP_TARGET, readings=target, usable=bool(all(v['usable'] for v in target.values())), beside=beside,
        part_d1_due=bool(all(v['usable'] for v in target.values()))))
    return df


def _sat_map(x: np.ndarray, x_max: float) -> np.ndarray:
    """A3's saturating map from dF/F to rate (the bound for the indicator): x / (1 - x / (SAT_HEADROOM x_max)), x clipped
    below 0.99 SAT_HEADROOM x_max."""
    c = SAT_HEADROOM * x_max
    x = np.minimum(np.asarray(x, dtype=float), 0.99 * c)
    return x / (1.0 - x / c)


def _fly_mean(table: pd.DataFrame, gl: Sequence[str]) -> pd.DataFrame:
    """The mean over flies (glomeruli x odors) of a table's entries on `gl` (each entry over the flies that have it)."""
    return table[list(gl)].groupby(level='odor').mean().reindex(list(ODORS)).T


def _cv(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.std(x) / np.mean(x))


def al_tables(mapped: bool = False) -> Dict[str, pd.DataFrame]:
    """GH146 and pebbled tables on the 28 shared glomeruli, linear or through the saturating map (x_max of each stage's
    largest fly-mean entry)."""
    gl = shared_glomeruli()
    out = {}
    for k, kind in (('pn', 'GH146'), ('orn', 'pebbled')):
        t = al_table(kind)[gl]
        if mapped:
            xm = float(np.nanmax(_fly_mean(t, gl).to_numpy()))
            t = pd.DataFrame(_sat_map(t.to_numpy(), xm), index=t.index, columns=t.columns)
        out[k] = t
    return out


def al_statistics(pn: pd.DataFrame, orn: pd.DataFrame) -> Dict[str, object]:
    """A3's statistics of a pair of tables (fly x odor rows, glomeruli columns): per glomerulus the transfer (slope
    through the origin of the fly-mean PN entries on the fly-mean ORN entries over the 12 odors, and their r); per odor
    the summed response over the glomeruli (fly means), its CV over odors; the correlation of the entries; per pair the
    distances (mean over flies of per-fly 1 - r)."""
    gl = list(pn.columns)
    P, O = _fly_mean(pn, gl), _fly_mean(orn, gl)
    slope = (P * O).sum(axis=1) / (O * O).sum(axis=1)
    r_g = pd.Series({g: _corr(P.loc[g], O.loc[g]) for g in gl})
    tp, to = P.sum(axis=0), O.sum(axis=0)
    return dict(slope=slope, r_glom=r_g, total_pn=tp, total_orn=to, cv_total_pn=_cv(tp), cv_total_orn=_cv(to),
                r_totals=_corr(tp, to), r_entries=_corr(P.to_numpy().ravel(), O.to_numpy().ravel()),
                dist_pn=per_fly_distances(pn, gl).mean(axis=0), dist_orn=per_fly_distances(orn, gl).mean(axis=0), P=P, O=O)


def _resample_flies(t: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    flies = list(dict.fromkeys(t.index.get_level_values('fly')))
    pick = rng.integers(0, len(flies), len(flies))
    parts = []
    for j, i in enumerate(pick):
        x = t.xs(flies[i], level='fly', drop_level=True)
        x.index = pd.MultiIndex.from_product([[f'r{j}'], x.index], names=['fly', 'odor'])
        parts.append(x)
    return pd.concat(parts)


def part_a3(n_boot: int = 200) -> pd.DataFrame:
    """al_transfer.csv: A3's statistics (`al_statistics`) of the GH146 PN dendrites against the pebbled ORNs on the 28
    shared glomeruli, linear and through the saturating map; 90 % intervals from n_boot resamples of each stage's flies
    (fewer than A1's: each resample recomputes the per-fly distances). Long format: section, map, key, stat, value, lo,
    hi."""
    rows = []
    for mapped in (False, True):
        T = al_tables(mapped)
        st = al_statistics(T['pn'], T['orn'])
        rng = np.random.default_rng(BOOT_SEED)
        boots = [al_statistics(_resample_flies(T['pn'], rng), _resample_flies(T['orn'], rng)) for _ in range(n_boot)]
        mp = 'saturating' if mapped else 'linear'

        def add(section, key, stat, value, bvals):
            bv = np.asarray(bvals, dtype=float)
            lo, hi = (np.nanpercentile(bv, CI) if np.isfinite(bv).any() else (np.nan, np.nan))
            rows.append(dict(section=section, map=mp, key=key, stat=stat, value=float(value), lo=float(lo), hi=float(hi)))
        for k in ('cv_total_pn', 'cv_total_orn', 'r_totals', 'r_entries'):
            add('summary', 'all', k, st[k], [b[k] for b in boots])
        for g in st['slope'].index:
            add('transfer', g, 'slope', st['slope'][g], [b['slope'][g] for b in boots])
            add('transfer', g, 'r', st['r_glom'][g], [b['r_glom'][g] for b in boots])
        for o in ODORS:
            add('totals', o, 'total_pn', st['total_pn'][o], [b['total_pn'][o] for b in boots])
            add('totals', o, 'total_orn', st['total_orn'][o], [b['total_orn'][o] for b in boots])
        for which in ('pn', 'orn'):
            d = st[f'dist_{which}']
            for (a, b), v in d.items():
                add('pairs', f'{a}|{b}', f'dist_{which}', v, [bb[f'dist_{which}'][(a, b)] for bb in boots])
            gm = group_means(d)
            gb = [group_means(bb[f'dist_{which}']) for bb in boots]
            for g in GROUPS:
                add('groups', g, f'dist_{which}', gm[g], [x[g] for x in gb])
        # the fly-mean entries themselves (B1 compares the model's PN stage with them)
        for g in st['P'].index:
            for o in ODORS:
                rows.append(dict(section='entries', map=mp, key=f'{g}|{o}', stat='pn', value=float(st['P'].loc[g, o]),
                                 lo=np.nan, hi=np.nan))
                rows.append(dict(section='entries', map=mp, key=f'{g}|{o}', stat='orn', value=float(st['O'].loc[g, o]),
                                 lo=np.nan, hi=np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'al_transfer.csv', index=False)
    return df


def part_a() -> dict:
    """Part A: A1, A2, A3, the identity checks of the readers (part_a.json) and the provenance."""
    t0 = time.time()
    checks = {}
    # the per-fly window means reproduce data/remy_rdm_timeseries/windows/ (closed 0.5-3 s windows)
    for st, name in (('bouton', 'pn_boutons'), ('claw', 'kc_claws')):
        wm = remy_window_means(st, (0.5, 3.0))
        _, _, odors = remy_matrices(st, 'stim')
        dmax = 0.0
        for i in range(wm.shape[0]):
            c = pd.read_csv(remy_dir() / 'windows' / f'{name}_fly{i + 1}_corrdist_odor_0.5-3.csv', index_col=0)
            dmax = max(dmax, float(np.max(np.abs(c.loc[odors, odors].to_numpy() - wm[i]))))
        checks[f'{st}_windows_max_abs_diff'] = dmax
    for st in ('bouton', 'claw'):
        one = remy_single_trial(st)
        if one is not None:
            a, t, odors = remy_matrices(st, 'blockavg')
            if one['odors'] != odors:
                raise RuntimeError('the explainer\'s odors differ from the matrices\'')
            blk = 1.0 - a[:, _frames(t, STAGE_WINDOW)].mean(axis=1)
            checks[f'{st}_single_trial_vs_blockavg_max_abs_diff'] = float(np.max(np.abs(one['between'] - blk)))
            checks[f'{st}_single_trial_frames'] = one['n_frames']
    S = stage_samples()
    pb36 = per_fly_distances(al_table('pebbled'), list(H.apl16_orn_deltas().index))
    checks['orn_pebbled_vs_olfdb_max_abs_diff'] = float(np.max(np.abs(pb36.mean(axis=0) - flies_stage_mean('orn'))))
    orn_c = _cs(H.orn16_consensus_pairs()).reindex(megamat_pairs())
    checks['orn_vs_orn16_consensus_max_abs_diff'] = float(np.max(np.abs(flies_stage_mean('orn') - orn_c)))
    ref = H.reference_pairs()
    ref.index = pd.MultiIndex.from_tuples([_canon(*p) for p in ref.index], names=['odor_a', 'odor_b'])
    ref = ref.sort_index().reindex(megamat_pairs())
    for k, col in (('soma_remy', 'kc_soma'), ('pn_dend', 'pn_dend'), ('bouton_measured', 'pn_bout'),
                   ('claw_measured', 'kc_claw'), ('orn', 'orn')):
        m = flies_stage_mean(k)
        checks[f'{k}_vs_reference_{col}'] = dict(max_abs_diff=float(np.max(np.abs(m - ref[col]))), r=_corr(m, ref[col]))
    checks['n'] = {k: int(v.shape[0]) for k, v in S.items()}
    checks['shared_glomeruli'] = shared_glomeruli()
    rel = {st: remy_reliability(st) for st in ('bouton', 'claw')}
    checks['reliability_mean'] = {st: float(r[list(ODORS)].to_numpy().mean()) for st, r in rel.items()}
    a1 = part_a1()
    a2 = part_a2()
    a3 = part_a3()
    _write_json(out_dir() / 'part_a.json', checks)
    _write_json(out_dir() / 'provenance_A.json', dict(part='A', script=H.script_path(__file__), **H.code_commits(),
        settings=dict(STAGE_WINDOW=STAGE_WINDOW, N_BOOT=N_BOOT, CI=CI, BOOT_SEED=BOOT_SEED, REMY_VARIANT=REMY_VARIANT,
                      SAT_HEADROOM=SAT_HEADROOM, STEP_TARGET=STEP_TARGET, STEP_TARGET_READINGS=STEP_TARGET_READINGS,
                      SOMA_HALF_SEED=SOMA_HALF_SEED),
        inputs=dict(remy=str(remy_dir()), gh146=str(gh146_dir()), yang=str(H.yang_data_dir()),
                    orn='olfdb dataset 1 (E85_orn_fly_population.flies)'),
        seconds=round(time.time() - t0, 1)))
    print(json.dumps({k: v for k, v in checks.items() if k != 'shared_glomeruli'}, indent=1, default=str), flush=True)
    return dict(checks=checks, a1=a1, a2=a2, a3=a3)


# ---------------------------------------------------------------------------
# B: the model's stages

ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's ALL_WIRINGS (S then T)
CONDS = (('on', 'on'), ('off', 'acute_off'))
STAGES = ('pn', 'bouton', 'claw')


def claw_inputs(W: np.ndarray, boutons: np.ndarray, inh: np.ndarray, wAPLKC: np.ndarray) -> np.ndarray:
    """The claws' net drive as olfsysm's sim_KC_layer forms it (no claw dynamics, subtractive APL inhibition):
    wiring x bouton rate at step t minus the claw's scaled APL weight times APL's activity at step t - 1, clipped at 0
    (claws x steps). `boutons` boutons x steps (olfsysm's bouton_sims, after APL's inhibition at the bouton), `inh` APL's
    activity at the steps before (steps; olfsysm's inh_sims shifted by one), `wAPLKC` the claws' scaled APL weights."""
    x = W @ boutons
    x -= np.outer(np.asarray(wAPLKC, dtype=float).ravel(), np.asarray(inh, dtype=float).ravel())
    np.maximum(x, 0.0, out=x)
    return x


def _windows(t: np.ndarray, window: Tuple[float, float]) -> Dict[str, np.ndarray]:
    """The stage window's steps (closed, as E78's who_selects reads 0-0.5 s) and the pre-odor window's (PRE_WINDOW,
    open on the left: olfsysm's first step of the run has no bouton value)."""
    return dict(win=np.flatnonzero((t >= window[0] - 1e-9) & (t <= window[1] + 1e-9)),
                pre=np.flatnonzero((t > PRE_WINDOW[0] + 1e-9) & (t < PRE_WINDOW[1] - 1e-9)))


def step_times(mp, rv) -> np.ndarray:
    """The time of every simulation step, pre_start + k dt (or start + k dt once the pretime is removed), from olfsysm's
    parameters. (`get_time_index` spreads the steps evenly from the first time to time_end inclusive, which puts every
    step a fraction of a step later than olfsysm's, enough to move a window's first step.)"""
    n = int(np.asarray(rv.orn.sims[0]).shape[-1])
    dt = float(mp.time_dt)
    for t0 in (float(mp.time_pre_start), float(mp.time_start)):
        if n == int(round((float(mp.time_end) - t0) / dt)):
            return t0 + np.arange(n) * dt
    raise RuntimeError(f'{n} steps match neither the pretime nor the start')


def glomerulus_of_claws(W: pd.DataFrame) -> np.ndarray:
    """Each claw's glomerulus (that of its largest PN->claw weight)."""
    gl = np.asarray(W.columns.get_level_values('glomerulus'))
    return gl[np.asarray(W.to_numpy(dtype=float).argmax(axis=1))]


def read_stages(mp, rv, W: pd.DataFrame, order: np.ndarray, odors: Sequence[str], window: Tuple[float, float]) -> dict:
    """The stage matrices of the last simulation (olfsysm's variables `mp`, `rv`; the simulator's output column j is
    olfsysm's odor order[j]): per stage ('pn' glomeruli, 'bouton' boutons, 'claw' claws; x odors) the window means
    ('<stage>_win') and the pre-odor means ('<stage>_pre'), and 'apl' (APL's activity, window and pre-odor means, per
    odor)."""
    from pyolfsysm.mb_model import get_sim_var
    t = step_times(mp, rv)
    w = _windows(t, window)
    pn = get_sim_var(rv.pn, 'pn_sims')
    gl = list(dict.fromkeys(W.columns.get_level_values('glomerulus')))
    if np.asarray(pn[0]).shape[0] != len(gl):
        raise RuntimeError(f'{np.asarray(pn[0]).shape[0]} PN rows for {len(gl)} glomeruli in the wiring')
    btn = get_sim_var(rv.pn, 'bouton_sims')
    inh = get_sim_var(rv.kc, 'inh_sims')
    if btn.shape[1] != W.shape[1] or inh.shape[0] != btn.shape[0]:
        raise RuntimeError(f'bouton_sims {btn.shape}, inh_sims {inh.shape}, wiring {W.shape}')
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    if len(wk) != W.shape[0]:
        raise RuntimeError(f'{len(wk)} APL->claw weights for {W.shape[0]} claws')
    Wm = W.to_numpy(dtype=float)
    out = {k: {} for k in ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')}
    apl = {}
    for j, o in enumerate(odors):
        i = int(order[j])
        for k in ('win', 'pre'):
            idx = w[k]
            out[f'pn_{k}'][o] = np.asarray(pn[i])[:, idx].mean(axis=1)
            b = np.asarray(btn[i])[:, idx]
            if np.isnan(b).any():
                raise RuntimeError(f'NaN bouton rates in the {k} window')
            out[f'bouton_{k}'][o] = b.mean(axis=1)
            a_prev = np.asarray(inh[i]).ravel()[idx - 1]
            out[f'claw_{k}'][o] = claw_inputs(Wm, b, a_prev, wk).mean(axis=1)
        a = np.asarray(inh[i]).ravel()
        apl[o] = dict(win=float(a[w['win']].mean()), pre=float(a[w['pre']].mean()))
    res = {}
    for k, v in out.items():
        idx = (pd.Index(gl, name='glomerulus') if k.startswith('pn') else W.columns if k.startswith('bouton') else W.index)
        res[k] = pd.DataFrame(v, index=idx)
    res['apl'] = pd.DataFrame(apl).T
    res['n_steps'] = dict(win=int(len(w['win'])), pre=int(len(w['pre'])), dt=float(t[1] - t[0]))
    return res


def _save_stage(d: Path, tag: str, st: dict, o: Tuple[pd.DataFrame, pd.DataFrame]) -> None:
    for k in ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre'):
        x = st[k].copy()
        x.columns = [str(c) for c in x.columns]
        if isinstance(x.index, pd.MultiIndex):
            x.index = x.index.set_names([str(n) for n in x.index.names])
        x.to_parquet(d / f'{tag}_{k}.parquet')
    st['apl'].to_csv(d / f'{tag}_apl.csv')
    o[0].to_parquet(d / f'{tag}_responses.parquet')
    o[1].to_parquet(d / f'{tag}_spike_counts.parquet')


def load_stage(part: str, tag: str) -> Optional[dict]:
    """A saved run's stage matrices and outputs (`_save_stage`), from the local outputs or the archive; None if absent."""
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        f = d / f'{tag}_claw_win.parquet'
        if f.exists():
            st = {k: pd.read_parquet(d / f'{tag}_{k}.parquet') for k in
                  ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')}
            st['apl'] = pd.read_csv(d / f'{tag}_apl.csv', index_col=0)
            st['responses'] = pd.read_parquet(d / f'{tag}_responses.parquet')
            st['spike_counts'] = pd.read_parquet(d / f'{tag}_spike_counts.parquet')
            return st
    return None


class StageRun:
    """One wiring at its fit (C2's family; E83's `build_simulator`, whose closure holds olfsysm's variables), APL's
    activity saved (`mp.kc.save_inh_sims`); `run(cond)` simulates the condition at the fit's parameters (APL off: the
    APL->KC and APL->bouton scales 0, as E83's `_check_saved`), checks the counts and responses against the saved fit's
    exactly, and returns its outputs and stage matrices (`read_stages`)."""

    def __init__(self, label: str):
        from experiments import E83_cross_wiring_fit as E83
        P = E83.wiring_parts(label)
        self.label, self.p, self.W, self.saved = label, dict(P['p']), P['W'], P['saved']
        self.sim, self.info = E83.build_simulator(label)
        cells = self.info['cells']
        self.mp, self.rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        self.order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        self.mp.kc.save_inh_sims = True
        wp = np.asarray(self.rv.kc.wPNKC, dtype=float)
        self.wiring_equal = bool(wp.shape == self.W.shape and np.array_equal(wp, self.W.to_numpy(dtype=float)))
        if not self.wiring_equal:
            raise RuntimeError(f"{label}: olfsysm's wPNKC differs from the wiring's claw x bouton matrix")

    def params(self, cond: str) -> Dict[str, float]:
        return dict(self.p) if cond == 'on' else dict(self.p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)

    def run(self, cond: str) -> dict:
        r, s = self.sim(self.params(cond))
        rs, ss = self.saved[cond]
        nd = int((s.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum())
        nr = int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum())
        odors = [H.odor_abbrev(c) for c in r.columns]
        st = read_stages(self.mp, self.rv, self.W, self.order, odors, MODEL_WINDOW)
        st.update(responses=r, spike_counts=s, identity=dict(counts_differ=nd, responses_differ=nr,
                  index_equal=bool(s.index.equals(ss.index) and s.columns.equals(ss.columns)),
                  exact=bool(nd == 0 and nr == 0 and s.index.equals(ss.index))))
        return st


def claw_identity(odors: Sequence[str] = CLAW_IDENTITY_ODORS) -> dict:
    """B0's claw identity: the real wiring at its fit (C2) with the lab's input restricted to `odors` (the breadth held
    at the full input's), built with olfsysm's dynamics recorded (`return_dynamics=True`, pretime kept so the simulator
    can re-run), APL on: its spike counts must equal C2's saved counts of those odors, and the claw inputs formed by
    `claw_inputs` from its bouton_sims, inh_sims and APL->claw weights must equal its claw_sims at every step after the
    first."""
    from experiments import E10_joint_fit as J
    from experiments import E83_cross_wiring_fit as E83
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator, get_sim_var
    P = E83.wiring_parts('real')
    orn = P['orn']
    cols = [c for c in orn.columns if H.odor_abbrev(c) in odors]
    wnorm, W = P['wnorm'], P['W']
    p = dict(P['p'])
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], P['breadth'], E83.BREADTH_K)
    sim, _ = make_fixed_param_kc_simulator(orn[cols], initial_params=p, **J.base_kws(None), _wPNKC=W, _wAPLKC=wAPLKC,
                                           _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'], _wPNAPL=wnorm['wPNAPL'],
                                           kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU, return_dynamics=True,
                                           delete_pretime=False)
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    r, s = sim(p)
    saved = P['saved']['on'][1]
    sv = H.abbrev_columns(saved)
    sc = H.abbrev_columns(s)
    counts_differ = int((sc.to_numpy(dtype=float) != sv[list(sc.columns)].to_numpy(dtype=float)).sum())
    claws = get_sim_var(rv.kc, 'claw_sims')
    btn = get_sim_var(rv.pn, 'bouton_sims')
    inh = get_sim_var(rv.kc, 'inh_sims')
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    Wm = W.to_numpy(dtype=float)
    res = dict(odors=list(sc.columns), counts_differ=counts_differ, n_claws=int(Wm.shape[0]), n_boutons=int(Wm.shape[1]),
               steps=int(claws.shape[-1]))
    nd, mx, nsteps = 0, 0.0, 0
    for j in range(len(order)):
        i = int(order[j])
        b = np.asarray(btn[i])
        ok = np.flatnonzero(~np.isnan(b).any(axis=0))
        ok = ok[ok >= 1]
        x = claw_inputs(Wm, b[:, ok], np.asarray(inh[i]).ravel()[ok - 1], wk)
        c = np.asarray(claws[i])[:, ok]
        nd += int((x != c).sum())
        mx = max(mx, float(np.max(np.abs(x - c))))
        nsteps += int(len(ok))
    res.update(claw_values_differ=nd, claw_max_abs_diff=mx, steps_compared=nsteps,
               one_bouton_per_claw=bool(((Wm != 0).sum(axis=1) == 1).all()),
               exact=bool(counts_differ == 0 and nd == 0))
    del sim, cells, mp, rv, claws, btn
    return res


def part_b0() -> dict:
    """B0: every wiring's run with the rates saved reproduces its saved fit's counts exactly, APL on and off (the runs
    are part B's, saved as they are made), and the claw identity on the real wiring; identity.json."""
    import gc
    t0 = time.time()
    path = out_dir() / 'identity.json'
    res = json.loads(path.read_text()) if path.exists() else {}
    if 'claw_identity' not in res:
        res['claw_identity'] = claw_identity()
        gc.collect()
        _write_json(path, res)
        print(f"claw identity: {res['claw_identity']}", flush=True)
    part_b()
    res = json.loads(path.read_text())
    res['all_exact'] = bool(res['claw_identity']['exact'] and all(v.get('exact') for k, v in res.items()
                                                                   if k.startswith('run_')))
    _write_json(path, res)
    print(f"B0 done in {time.time() - t0:.0f} s; all exact: {res['all_exact']}", flush=True)
    return res


def part_b(wirings: Optional[Sequence[str]] = None) -> None:
    """Part B: each wiring's runs (APL on and off) with the stage matrices, saved under outputs/E90_model_stages/B/ as
    <wiring>_<on|acute_off>_<stage>_<win|pre>.parquet with the responses and counts (resumable: wirings whose four
    files exist are skipped); each run's identity against its saved fit goes to identity.json (run_<wiring>_<cond>)."""
    import gc
    d = local_dir('B')
    path = out_dir() / 'identity.json'
    for w in (wirings or ALL_WIRINGS):
        if all((d / f'{w}_{cc}_claw_win.parquet').exists() for _, cc in CONDS):
            continue
        t0 = time.time()
        sr = StageRun(w)
        unit_meta(w, sr.W)
        for c, cc in CONDS:
            st = sr.run(c)
            _save_stage(d, f'{w}_{cc}', st, (st['responses'], st['spike_counts']))
            res = json.loads(path.read_text()) if path.exists() else {}
            res[f'run_{w}_{c}'] = dict(st['identity'], n_steps=st['n_steps'], wiring_equal=sr.wiring_equal)
            _write_json(path, res)
            if not st['identity']['exact']:
                raise RuntimeError(f'{w} {c}: the run differs from the saved fit: {st["identity"]}')
            print(f'{w} {c}: exact, {st["n_steps"]}', flush=True)
            del st
        del sr
        gc.collect()
        print(f'{w}: {time.time() - t0:.0f} s', flush=True)
    _write_json(out_dir() / 'provenance_B.json', dict(part='B', script=H.script_path(__file__),
        model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), wirings=list(wirings or ALL_WIRINGS),
        settings=dict(MODEL_WINDOW=MODEL_WINDOW, PRE_WINDOW=PRE_WINDOW)))


def b5_spec() -> dict:
    """E82's B1_C2 point (the real wiring at 2 s, re-held) and its parameters (runs_params.json)."""
    from experiments import E82_protocol_runs as R
    spec = next(x for x in R.point_specs() if x['point_id'] == 'B1_C2')
    rec = json.loads((H.results_dir(R.EXPERIMENT) / 'runs_params.json').read_text())['B1_C2']
    return dict(spec=spec, params=rec['params'], params_off=rec['params_off'], rec=rec)


def b5_saved() -> Dict[str, tuple]:
    """E82's saved B1_C2 outputs (local outputs, else the archive)."""
    from experiments import E82_protocol_runs as R
    for d in (H.outputs_dir(R.EXPERIMENT) / 'runs', H.shared_dir() / 'data' / 'model_outputs' / R.EXPERIMENT / 'runs'):
        if (d / 'B1_C2_on_spike_counts.parquet').exists():
            return {c: tuple(pd.read_parquet(d / f'B1_C2_{cc}_{k}.parquet') for k in ('responses', 'spike_counts'))
                    for c, cc in CONDS}
    raise RuntimeError("E82's B1_C2 outputs not found")


def part_b5() -> dict:
    """B5: the real wiring at the flies' 2 s pulse at E82's re-held B1 point: built as E82_protocol_runs builds it (its
    breadth pre-pass at 2 s), run at the point's parameters APL on and off; the counts must reproduce E82's saved B1_C2
    counts exactly before the stages are read (over PULSE_WINDOW); saved under outputs/E90_model_stages/B5/."""
    import gc
    from experiments import E82_protocol_runs as R
    t0 = time.time()
    b = b5_spec()
    spec = b['spec']
    m = R.model_parts('C2')
    breadth = R.point_breadth(spec, m)
    from experiments import E80_pulse as E80P
    if E80P._sha(breadth.to_numpy()) != b['rec']['breadth_sha256']:
        raise RuntimeError("B5: the breadth pre-pass at 2 s differs from E82's B1_C2 point's")
    sim, info = R.build_simulator(R.model_kws(spec, m), wPNKC=None, breadth=breadth, params=b['params'],
                                  kc_kc_tau=spec['kc_kc_tau'])
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    mp.kc.save_inh_sims = True
    saved = b5_saved()
    W = info['wPNKC']
    d = local_dir('B5')
    res = dict(point='E82 B1_C2', params=b['params'], params_off=b['params_off'])
    for c, cc in CONDS:
        r, s = sim(b['params'] if c == 'on' else b['params_off'])
        rs, ss = saved[c]
        nd = int((s.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum())
        nr = int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum())
        res[f'{c}_counts_differ'], res[f'{c}_responses_differ'] = nd, nr
        if nd or nr:
            raise RuntimeError(f'B5 {c}: {nd} counts and {nr} responses differ from E82\'s B1_C2')
        odors = [H.odor_abbrev(x) for x in r.columns]
        st = read_stages(mp, rv, W, order, odors, PULSE_WINDOW)
        _save_stage(d, f'real_pulse2_{cc}', st, (r, s))
        res[f'{c}_n_steps'] = st['n_steps']
        print(f'B5 {c}: exact; {st["n_steps"]}', flush=True)
    res['exact'] = True
    _write_json(out_dir() / 'b5_point.json', res)
    _write_json(out_dir() / 'provenance_B5.json', dict(part='B5', script=H.script_path(__file__),
        model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), spec={k: v for k, v in spec.items()},
        settings=dict(PULSE=PULSE, PULSE_WINDOW=PULSE_WINDOW, PRE_WINDOW=PRE_WINDOW), seconds=round(time.time() - t0, 1)))
    del sim, info, cells, mp, rv
    gc.collect()
    return res


# ---------------------------------------------------------------------------
# B: readings

# the session z-score's share of the session's frames inside the response windows: Remy's 3 s pulse over her trials'
# 25 s of frames (-5 to 20 s), the model's response taken as confined to its window
Q_SESSION = 3.0 / 25.0
SUBSETS = ('shared', 'all')             # the shared glomeruli's boutons / claws, and all of them
MODEL_READINGS = ('rate', 'session_z', 'noisy_nr', 'noisy_measured')
NOISE_FREE = ('rate', 'session_z')


def disattenuate(r: np.ndarray) -> np.ndarray:
    """1 - r_ij / sqrt(r_ii r_jj) of a (..., n, n) correlation matrix whose diagonal holds the reliabilities (NaN where
    their product is not positive)."""
    rel = np.diagonal(r, axis1=-2, axis2=-1)
    den = rel[..., :, None] * rel[..., None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, 1.0 - r / np.sqrt(np.where(den > 0, den, 1.0)), np.nan)


def block_average_distances(trials: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """A pseudo-fly's two readings from its trials (trials x odors x units), as Remy's matrices are made: 'as measured',
    1 - r of the trial means (her `stim` matrices); 'noise removed', `disattenuate` of the block-average correlations
    (between two odors the mean over all pairs of their trials of the single-trial correlations, on the diagonal the mean
    over the pairs of different trials of the odor; her `blockavg` matrices). Returns two odors x odors arrays."""
    x = np.asarray(trials, dtype=float)
    T = x.shape[0]
    z = x - x.mean(axis=2, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = z / np.linalg.norm(z, axis=2, keepdims=True)
    C = np.einsum('aiu,bju->abij', z, z)                      # trials x trials x odors x odors
    r = C.mean(axis=(0, 1))
    off = ~np.eye(T, dtype=bool)
    d = np.einsum('abii->abi', C)[off].mean(axis=0)
    r[np.diag_indices_from(r)] = d
    return _pearson_dist_matrix(x.mean(axis=0)), disattenuate(r)


def session_z(values: pd.DataFrame) -> pd.DataFrame:
    """Each unit divided by its SD over a simulated session of the panel's odors and one blank, the response confined to
    a share Q_SESSION of the frames (`pyolfsysm.readout.session_zscore`, linear amplitude, no noise; units silent through
    the session stay 0)."""
    from pyolfsysm.readout import session_zscore
    z = session_zscore(values.to_numpy(dtype=float), q=Q_SESSION, noise=0.0, n_trials=1, n_blank=1)
    return pd.DataFrame(z['center'][0], index=values.index, columns=values.columns)


def _drop_silent(x: pd.DataFrame) -> pd.DataFrame:
    """Units that are 0 for every odor are left out (Remy's units are responsive ROIs; they would only add a constant)."""
    return x.loc[(x.abs() > 1e-12).any(axis=1)]


def noisy_reading(values: pd.DataFrame, reliability: pd.Series, n_flies: int, rng: np.random.Generator,
                  n_reps: int = N_NOISE_REPS) -> dict:
    """`values` (units x the 12 megamat odors, session z) through measurement noise at the flies' per-odor block-to-block
    reliability: each odor's noise SD set so that the expected correlation between two trials across units is the
    flies' reliability r_o (sigma^2 = var over units x (1 - r_o) / r_o); per pseudo-fly three trials, read as the flies'
    (`block_average_distances`); the mean over `n_flies` pseudo-flies, as the flies' stage is their mean; `n_reps`
    replicates. Returns 'measured' / 'noise_removed' (pair Series, the mean over replicates) and their SD over
    replicates per pair ('measured_sd', 'noise_removed_sd')."""
    X = values[list(ODORS)].to_numpy(dtype=float).T                       # odors x units
    r = reliability.reindex(list(ODORS)).to_numpy(dtype=float)
    sig = np.sqrt(X.var(axis=1) * (1.0 - r) / r)
    out = {'measured': [], 'noise_removed': []}
    for _ in range(n_reps):
        m_acc, n_acc = np.zeros((len(ODORS), len(ODORS))), np.zeros((len(ODORS), len(ODORS)))
        for _f in range(n_flies):
            tr = X[None] + sig[None, :, None] * rng.standard_normal((3,) + X.shape)
            m, nr = block_average_distances(tr)
            m_acc += m
            n_acc += nr
        out['measured'].append(_pairs_from_square(m_acc / n_flies, ODORS))
        out['noise_removed'].append(_pairs_from_square(n_acc / n_flies, ODORS))
    res = {}
    for k, v in out.items():
        df = pd.concat(v, axis=1)
        res[k], res[f'{k}_sd'] = df.mean(axis=1), df.std(axis=1, ddof=1)
    return res


def flies_reliability(stage: str) -> pd.Series:
    """The flies' per-odor block-to-block reliability of a stage (`remy_reliability`, mean over flies)."""
    return remy_reliability(stage).mean(axis=0)


def claw_signal(change: pd.DataFrame, counts: pd.DataFrame, kc_of_claw: np.ndarray, share: float,
                scales: Tuple[float, float]) -> pd.DataFrame:
    """B3's version of what a claw's calcium reports: (1 - share) x its net drive's change / scales[0] + share x its KC's
    spike count / scales[1] (scales: the root mean square of each over all claws and odors in the wiring's APL-on run, so
    that share is the spike signal's share of the amplitude); share 0 is the drive alone."""
    if share == 0:
        return change
    ix = pd.Index(counts.index.get_level_values(H.KC_ID_LEVEL))
    k = counts.to_numpy(dtype=float)[ix.get_indexer(kc_of_claw)]
    if (ix.get_indexer(kc_of_claw) < 0).any():
        raise RuntimeError('claws whose KC is not in the counts')
    k = pd.DataFrame(k, index=change.index, columns=[H.odor_abbrev(c) for c in counts.columns])[list(change.columns)]
    return (1.0 - share) * change / scales[0] + share * k / scales[1]


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def stage_units(st: dict, stage: str) -> pd.DataFrame:
    """A run's stage as units x odors changes (window mean minus pre-odor mean), odor columns abbreviated."""
    x = st[f'{stage}_win'] - st[f'{stage}_pre']
    x.columns = [H.odor_abbrev(c) for c in x.columns]
    return x


def subset_mask(stage: str, units: pd.Index, meta: Optional[dict], subset: str) -> np.ndarray:
    """The units of `subset`: 'all', or 'shared' (those of the 28 shared glomeruli: a bouton's glomerulus, a claw's
    (`unit_meta`), a PN's own)."""
    if subset == 'all':
        return np.ones(len(units), dtype=bool)
    sh = list(shared_glomeruli())
    if stage == 'bouton':
        return np.isin(meta['bouton_glomerulus'], sh)
    if stage == 'claw':
        return np.isin(meta['claw_glomerulus'], sh)
    if stage == 'pn':
        return np.asarray(units.isin(sh))
    raise ValueError(stage)


def wiring_matrix(label: str) -> pd.DataFrame:
    """A wiring's claw x bouton matrix (E83's `wiring_parts`, which may simulate the base's APL-off run once)."""
    from experiments import E83_cross_wiring_fit as E83
    return E83.wiring_parts(label)['W']


def unit_meta(label: str, W: Optional[pd.DataFrame] = None) -> dict:
    """A wiring's units' metadata, cached under outputs/E90_model_stages/B/<label>_units.json: 'bouton_glomerulus' (in
    the wiring's bouton order), 'claw_glomerulus' and 'kc_of_claw' (in its claw order). `label` a part B wiring or a D1
    model's wiring."""
    path = local_dir('B') / f'{label}_units.json'
    if path.exists():
        d = json.loads(path.read_text())
        return {k: np.asarray(v) for k, v in d.items()}
    W = wiring_matrix(label) if W is None else W
    d = dict(bouton_glomerulus=[str(x) for x in W.columns.get_level_values('glomerulus')],
             claw_glomerulus=[str(x) for x in glomerulus_of_claws(W)],
             kc_of_claw=[int(x) for x in W.index.get_level_values(H.KC_ID_LEVEL)])
    path.write_text(json.dumps(d) + '\n')
    return {k: np.asarray(v) for k, v in d.items()}


def read_wiring(label: str, st_on: dict, st_off: dict, meta: dict, rng: np.random.Generator,
                with_noise: bool = True) -> List[dict]:
    """Every reading of one wiring's bouton and claw stages (APL on and off; SUBSETS; claws in CLAW_SPIKE_SHARES): rows of
    stage, subset, share, cond, reading, odor_a, odor_b, value (pair distance, 1 - Pearson across units), plus the noise
    check's truth / estimate pairs ('noisy_*' readings carry their replicate SD as `sd`)."""
    rows = []
    rel = {st: flies_reliability(st) for st in ('bouton', 'claw')}
    nflies = {'bouton': stage_samples()['bouton_measured'].shape[0], 'claw': stage_samples()['claw_measured'].shape[0]}
    kc_of_claw = np.asarray(meta['kc_of_claw'])
    for st in (st_on, st_off):
        if not np.array_equal(np.asarray(st['claw_win'].index.get_level_values(H.KC_ID_LEVEL)).astype(int), kc_of_claw.astype(int)):
            raise RuntimeError(f'{label}: the claws are not those of the units\' metadata')
        if len(st['bouton_win']) != len(meta['bouton_glomerulus']):
            raise RuntimeError(f'{label}: the boutons are not those of the units\' metadata')
    ch_on = stage_units(st_on, 'claw')
    scales = (_rms(ch_on.to_numpy()), _rms(st_on['spike_counts'].to_numpy(dtype=float)))
    for cond, st in (('on', st_on), ('off', st_off)):
        for stage in ('bouton', 'claw'):
            U = stage_units(st, stage)
            for subset in SUBSETS:
                m = subset_mask(stage, U.index, meta, subset)
                for share in (CLAW_SPIKE_SHARES if stage == 'claw' else (0.0,)):
                    V = claw_signal(U, st['spike_counts'], kc_of_claw, share, scales) if stage == 'claw' else U
                    V = _drop_silent(V.loc[m])
                    Z = session_z(V)
                    key = dict(wiring=label, stage=stage, subset=subset, share=share, cond=cond, n_units=int(V.shape[0]))
                    for rd, d in (('rate', odor_distances(V)), ('session_z', odor_distances(Z))):
                        rows += [dict(key, reading=rd, odor_a=a, odor_b=b, value=float(v), sd=np.nan) for (a, b), v in d.items()]
                    if with_noise:
                        nz = noisy_reading(Z, rel[stage], nflies[stage], rng)
                        for rd, k in (('noisy_nr', 'noise_removed'), ('noisy_measured', 'measured')):
                            rows += [dict(key, reading=rd, odor_a=a, odor_b=b, value=float(v), sd=float(nz[f'{k}_sd'][(a, b)]))
                                     for (a, b), v in nz[k].items()]
    return rows


def pn_rows(label: str, st: dict, cond: str) -> List[dict]:
    """The PN stage of a run on the shared glomeruli and on all of the wiring's: pair distances of the glomeruli's rate
    changes ('change'), of their window rates ('rate') and of the rates minus each glomerulus's lowest over the 16 odors
    ('minsub')."""
    rows = []
    win = st['pn_win'].copy()
    win.columns = [H.odor_abbrev(c) for c in win.columns]
    pre = st['pn_pre'].copy()
    pre.columns = win.columns
    for subset in SUBSETS:
        m = subset_mask('pn', win.index, None, subset)
        for rd, X in (('change', win - pre), ('rate', win), ('minsub', win.sub(win.min(axis=1), axis=0))):
            X = _drop_silent(X.loc[m])
            d = odor_distances(X)
            rows += [dict(wiring=label, stage='pn', subset=subset, share=0.0, cond=cond, n_units=int(X.shape[0]), reading=rd,
                          odor_a=a, odor_b=b, value=float(v), sd=np.nan) for (a, b), v in d.items()]
    return rows


def soma_rows(label: str, counts: Dict[str, pd.DataFrame], K: float) -> List[dict]:
    """B4: the somata of a wiring (its saved counts APL on and off) as raw counts and under the calibrated readout (K on
    the APL-on counts, as E83's `readout_K`), 1 - Pearson and 1 - cosine across KCs."""
    from pyolfsysm.readout import readout
    rows = []
    for cond, c in counts.items():
        c = H.abbrev_columns(c.astype(float))
        for ro, X in (('raw', c), ('cal', readout(c, 'saturating', K))):
            for meas in ('pearson', 'cosine'):
                d = odor_distances(X, meas)
                rows += [dict(wiring=label, stage='soma', subset='all', share=0.0, cond=cond, n_units=int(X.shape[0]),
                              reading=f'{ro}_{meas}', odor_a=a, odor_b=b, value=float(v), sd=np.nan) for (a, b), v in d.items()]
    return rows


def model_pairs(force: bool = False) -> pd.DataFrame:
    """Every wiring's readings of every stage (part B's saved stage matrices; `read_wiring`, `pn_rows`, `soma_rows`),
    cached in outputs/E90_model_stages/model_pairs.parquet; per wiring one generator (NOISE_SEED + the wiring's index)."""
    from experiments import E83_cross_wiring_fit as E83
    path = local_dir() / 'model_pairs.parquet'
    if path.exists() and not force:
        return pd.read_parquet(path)
    rows = []
    for i, w in enumerate(ALL_WIRINGS):
        st = {c: load_stage('B', f'{w}_{cc}') for c, cc in CONDS}
        if any(v is None for v in st.values()):
            print(f'{w}: no stage matrices yet', flush=True)
            continue
        t0 = time.time()
        rows += read_wiring(w, st['on'], st['off'], unit_meta(w), np.random.default_rng(NOISE_SEED + i))
        for c in ('on', 'off'):
            rows += pn_rows(w, st[c], c)
        K = E83.readout_K({'on': (st['on']['responses'], st['on']['spike_counts'])})
        rows += soma_rows(w, {c: st[c]['spike_counts'] for c in ('on', 'off')}, K)
        print(f'{w}: read in {time.time() - t0:.0f} s', flush=True)
    df = pd.DataFrame(rows)
    df.to_parquet(path)
    return df


KEY = ['stage', 'subset', 'share', 'cond', 'reading']


def populations(mp: pd.DataFrame) -> pd.DataFrame:
    """The population (mean over the wirings present) of every reading's pair distances: rows KEY + pair, value, n."""
    g = mp.groupby(KEY + ['odor_a', 'odor_b'])['value']
    return pd.DataFrame({'value': g.mean(), 'n_wirings': g.size()}).reset_index()


def _group_table(df: pd.DataFrame, by: Sequence[str]) -> pd.DataFrame:
    """Group means (GROUPS) of pair rows grouped by `by`."""
    out = []
    for k, sub in df.groupby(list(by)):
        s = pd.Series(sub.value.to_numpy(), index=pd.MultiIndex.from_arrays([sub.odor_a, sub.odor_b]))
        out.append(dict(zip(by, k if isinstance(k, tuple) else (k,)), **group_means(s)))
    return pd.DataFrame(out)


def noise_check(mp: pd.DataFrame) -> pd.DataFrame:
    """The check of the noise removal on the model's known values: per stage, subset, share and condition, the population's
    group mean of the noise-removed estimate ('noisy_nr') against the truth (the noise-free session z reading), and per
    wiring; passes for a stage if within NOISE_CHECK_TOL for every group, subset, share and condition in the population
    (noise_check.csv)."""
    rows = []
    d = mp[mp.stage.isin(['bouton', 'claw'])]
    for by_wiring in (False, True):
        keys = ['stage', 'subset', 'share', 'cond'] + (['wiring'] if by_wiring else [])
        src = d if by_wiring else populations(d)
        for k, sub in src.groupby(keys):
            kk = dict(zip(keys, k))
            tr = sub[sub.reading == 'session_z']
            es = sub[sub.reading == 'noisy_nr']
            ms = sub[sub.reading == 'noisy_measured']
            if not len(es):
                continue
            def ser(x):
                return pd.Series(x.value.to_numpy(), index=pd.MultiIndex.from_arrays([x.odor_a, x.odor_b]))
            gt, ge, gm = group_means(ser(tr)), group_means(ser(es)), group_means(ser(ms))
            for g in GROUPS:
                rows.append(dict(kk, wiring=kk.get('wiring', 'population'), group=g, truth=gt[g], estimate=ge[g],
                                 measured=gm[g], error=ge[g] - gt[g], error_measured=gm[g] - gt[g],
                                 within=bool(abs(ge[g] - gt[g]) <= NOISE_CHECK_TOL)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'noise_check.csv', index=False)
    return df


def check_passes(nc: pd.DataFrame) -> Dict[str, bool]:
    """Whether the noise removal passes its check, per stage (the population's rows, every group, subset, share and
    condition)."""
    p = nc[nc.wiring == 'population']
    return {st: bool(p[p.stage == st].within.all()) for st in ('bouton', 'claw')}


def flies_group(key: str) -> pd.DataFrame:
    """A fly stage's group rows (stages_flies.csv): group, mean, lo, hi."""
    d = pd.read_csv(out_dir() / 'stages_flies.csv')
    return d[(d.key == key) & (d.kind == 'group')].set_index('group')[['mean', 'lo', 'hi']]


def flies_step_group(step: str, claw_reading: str) -> pd.DataFrame:
    d = pd.read_csv(out_dir() / 'steps_flies.csv')
    return d[(d.step == step) & (d.claw == f'claw_{claw_reading}') & (d.kind == 'group')].set_index('group')[['mean', 'lo', 'hi']]


# the comparisons the verdict reads, each a model reading against the flies': the noise-free readings and the noisy
# reading with the noise removed against the flies' noise-removed values (if the stage's check passes), the noisy reading
# as measured against the flies as measured (always; alone if the check fails)
def reading_pairs(passed: bool) -> List[Tuple[str, str]]:
    if passed:
        return [('rate', 'noise_removed'), ('session_z', 'noise_removed'), ('noisy_nr', 'noise_removed')]
    return [('noisy_measured', 'measured')]


SOMA_STEP_READING = 'cal_pearson'
SOMA_STEP_BESIDE = ('cal_cosine', 'raw_pearson', 'raw_cosine')
FLY_STEP_OF = {('on', 'pearson'): 'on_ctrl', ('off', 'pearson'): 'off_tnt', ('on', 'cosine'): 'on_ctrl_cos',
               ('off', 'cosine'): 'off_tnt_cos'}


def model_steps(mp: pd.DataFrame) -> pd.DataFrame:
    """steps_model.csv: per wiring (and the population) the steps per pair: the somata (each soma reading, APL on and
    off) minus the claw stage with APL on (each claw reading, subset and share); rows step ('on' / 'off'), soma, subset,
    share, reading, wiring, odor_a, odor_b, value (`model_steps_of`)."""
    df = model_steps_of(mp)
    df.to_csv(out_dir() / 'steps_model.csv', index=False)
    return df


def effect_errors() -> Dict[str, float]:
    """The population's error in the APL effect per group (the 66-pair groups): the mean over the group's pairs of the
    pattern error (model minus data per-pair effect, each with its mean over the 120 pairs removed) of the eight wirings'
    population under the calibrated readout, against the main target (E83's `pattern_error`, 'z8')."""
    from experiments import E83_cross_wiring_fit as E83
    ds = [E83.pair_distances(E83.saved_fit(w)['saved'], 'cal') for w in ALL_WIRINGS]
    err = E83.pattern_error(E83.population(ds), E83.targets()['z8'])
    m = group_masks(err.index)
    return {g: float(err[m[g]].mean()) for g in GROUPS}


def part_bscore() -> pd.DataFrame:
    """Part B's tables: model_pairs (every reading, cached), stages_model.csv (per stage, subset, share, condition and
    reading, the population's and each wiring's group means beside the flies'), noise_check.csv, steps_model.csv,
    pn_stage.csv."""
    t0 = time.time()
    mp = model_pairs()
    nc = noise_check(mp)
    steps = model_steps(mp)
    gw = _group_table(mp, KEY + ['wiring'])
    gp = _group_table(populations(mp), KEY)
    gp['wiring'] = 'population'
    sm = pd.concat([gp, gw], ignore_index=True)
    sm.to_csv(out_dir() / 'stages_model.csv', index=False)
    pn_stage_table(mp)
    print(f'Bscore in {time.time() - t0:.0f} s; noise check: {check_passes(nc)}', flush=True)
    return sm


def pn_stage_table(mp: pd.DataFrame) -> pd.DataFrame:
    """B1: the model's PN stage on the 28 shared glomeruli against A3 (the PN dendrites' fly-mean entries and the ORN
    input): per reading the entries' correlation with the PN dendrites and with the ORN input, the per-odor totals and
    their CV over the 12 odors, and the group means of the pair distances (pn_stage.csv)."""
    st = load_stage('B', 'real_on')
    a3 = pd.read_csv(out_dir() / 'al_transfer.csv')
    ent = a3[(a3.section == 'entries') & (a3['map'] == 'linear')]
    ent = ent.assign(glomerulus=ent.key.str.split('|').str[0], odor=ent.key.str.split('|').str[1])
    P = ent[ent.stat == 'pn'].pivot(index='glomerulus', columns='odor', values='value')[list(ODORS)]
    O = ent[ent.stat == 'orn'].pivot(index='glomerulus', columns='odor', values='value')[list(ODORS)]
    gl = list(P.index)
    win = st['pn_win'].copy()
    win.columns = [H.odor_abbrev(c) for c in win.columns]
    pre = st['pn_pre'].copy()
    pre.columns = win.columns
    orn = H.abbrev_columns(H.apl16_orn_deltas().astype(float))
    rows = []
    for rd, X in (('change', win - pre), ('rate', win), ('minsub', win.sub(win.min(axis=1), axis=0)),
                  ('orn_input', orn)):
        X = X.loc[gl, list(ODORS)]
        tot = X.sum(axis=0)
        row = dict(reading=rd, r_entries_pn_dend=_corr(X.to_numpy().ravel(), P.to_numpy().ravel()),
                   r_entries_orn_dff=_corr(X.to_numpy().ravel(), O.to_numpy().ravel()),
                   r_entries_orn_input=_corr(X.to_numpy().ravel(), orn.loc[gl, list(ODORS)].to_numpy().ravel()),
                   cv_totals=_cv(tot), r_totals_pn_dend=_corr(tot, P.sum(axis=0)),
                   **{f'total_{o}': float(tot[o]) for o in ODORS},
                   **{f'dist_{g}': v for g, v in group_means(odor_distances(X)).items()})
        rows.append(row)
    for k, X in (('pn_dend', P), ('orn_dff', O)):
        tot = X.sum(axis=0)
        rows.append(dict(reading=k, cv_totals=_cv(tot), **{f'total_{o}': float(tot[o]) for o in ODORS},
                         **{f'dist_{g}': v for g, v in group_means(per_fly_distances(al_table('GH146' if k == 'pn_dend' else 'pebbled'), gl).mean(axis=0)).items()}))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'pn_stage.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# C: the verdict

def differs(model_pop: float, model_wirings: Sequence[float], flies: float, lo: float, hi: float) -> dict:
    """The rule for one stage or step, group and reading: |population - flies| >= VERDICT_DIFF, the population outside
    the flies' 90 % interval [lo, hi], and at least VERDICT_SAME_SIDE of the wirings on the population's side of the
    flies' value."""
    d = float(model_pop - flies)
    side = np.sign(d)
    n_same = int(sum(np.sign(float(w) - flies) == side for w in model_wirings)) if side != 0 else 0
    big = bool(abs(d) >= VERDICT_DIFF)
    outside = bool(model_pop < lo or model_pop > hi)
    same = bool(n_same >= VERDICT_SAME_SIDE)
    return dict(diff=d, big=big, outside=outside, n_same_side=n_same, same_side=same, differs=bool(big and outside and same))


def element_result(rows: pd.DataFrame) -> dict:
    """One stage's or step's result for a group over its readings: 'differs' (every reading differs), 'mixed' (some do,
    some do not), 'n_differ', 'n_readings', 'sign' (of the population minus the flies, where it differs)."""
    k = int(rows['differs'].sum())
    n = int(len(rows))
    signs = sorted(int(x) for x in set(np.sign(rows['diff'][rows['differs']])))
    return dict(differs=bool(n > 0 and k == n), mixed=bool(0 < k < n), n_differ=k, n_readings=n, signs=signs)


def group_verdict(res: Dict[str, dict], effect_err: float) -> str:
    """The verdict of one group from its elements' results ('bouton', 'claw', 'step_on', 'step_off'; `element_result`)
    and its error in the APL effect, by the entry's rules (fixed before any run): 'upstream' if the bouton or the claw
    stage differs; 'not resolved' if a stage's readings disagree; 'both steps' if both steps differ; 'the step with APL
    on' if that step differs and the step without APL does not (in every reading), 'the step without APL' the reverse;
    'not resolved' if a step's readings disagree, or nothing differs while the group's error in the effect is above
    VERDICT_EFFECT_ERR; 'none' otherwise."""
    if res['bouton']['differs'] or res['claw']['differs']:
        return 'upstream'
    if res['bouton']['mixed'] or res['claw']['mixed']:
        return 'not resolved'
    on, off = res['step_on'], res['step_off']
    if on['differs'] and off['differs']:
        return 'both steps'
    if on['differs'] and not off['mixed']:
        return 'the step with APL on'
    if off['differs'] and not on['mixed']:
        return 'the step without APL'
    if on['mixed'] or off['mixed']:
        return 'not resolved'
    if abs(effect_err) > VERDICT_EFFECT_ERR:
        return 'not resolved'
    return 'none'


# the subset each stage is read on for the verdict, as the recordings were made: Remy's PN boutons are GH146-labelled
# boutons, so the boutons of the GH146 glomeruli the model's input shares; her claws are KC claws, which receive from
# every PN, so all claws (the steps likewise). The other subset is reported beside, and the verdict with both subsets as
# readings is given too (verdict.json's 'both_subsets').
VERDICT_SUBSETS = {'bouton': ('shared',), 'claw': ('all',)}


def verdict_rows(sm: pd.DataFrame, steps: pd.DataFrame, passed: Dict[str, bool],
                 subsets: Optional[Dict[str, Sequence[str]]] = None) -> pd.DataFrame:
    """Every (element, group, reading) comparison the verdict reads (`differs`), the model's APL-on stages against the
    flies' (wild type) and the model's steps against the flies' (Yang's control flies and APL output blocked (TNT) minus
    the claws, Pearson), on `subsets` (VERDICT_SUBSETS by default)."""
    subsets = VERDICT_SUBSETS if subsets is None else subsets
    rows = []
    for stage in ('bouton', 'claw'):
        for mr, fr in reading_pairs(passed[stage]):
            fl = flies_group(f'{stage}_{fr}')
            sub = sm[(sm.stage == stage) & (sm.cond == 'on') & (sm.reading == mr) & sm.subset.isin(subsets[stage])]
            for (subset, share), x in sub.groupby(['subset', 'share']):
                pop = x[x.wiring == 'population'].iloc[0]
                wir = x[x.wiring != 'population']
                for g in GROUPS:
                    r = differs(pop[g], wir[g].to_numpy(), fl.loc[g, 'mean'], fl.loc[g, 'lo'], fl.loc[g, 'hi'])
                    rows.append(dict(element=stage, group=g, model_reading=mr, flies_reading=fr, subset=subset, share=share,
                                     model=float(pop[g]), flies=float(fl.loc[g, 'mean']), lo=float(fl.loc[g, 'lo']),
                                     hi=float(fl.loc[g, 'hi']), n_wirings=int(len(wir)), **r))
    gs = _group_table(steps, ['step', 'soma', 'subset', 'share', 'reading', 'wiring'])
    for step in ('on', 'off'):
        for mr, fr in reading_pairs(passed['claw']):
            fl = flies_step_group(FLY_STEP_OF[(step, 'pearson')], fr)
            sub = gs[(gs.step == step) & (gs.soma == SOMA_STEP_READING) & (gs.reading == mr) & gs.subset.isin(subsets['claw'])]
            for (subset, share), x in sub.groupby(['subset', 'share']):
                pop = x[x.wiring == 'population'].iloc[0]
                wir = x[x.wiring != 'population']
                for g in GROUPS:
                    r = differs(pop[g], wir[g].to_numpy(), fl.loc[g, 'mean'], fl.loc[g, 'lo'], fl.loc[g, 'hi'])
                    rows.append(dict(element=f'step_{step}', group=g, model_reading=mr, flies_reading=fr, subset=subset,
                                     share=share, model=float(pop[g]), flies=float(fl.loc[g, 'mean']), lo=float(fl.loc[g, 'lo']),
                                     hi=float(fl.loc[g, 'hi']), n_wirings=int(len(wir)), **r))
    return pd.DataFrame(rows)


def pn_verdict_rows(sm: pd.DataFrame) -> pd.DataFrame:
    """The PN stage against the PN dendrites by the same rule (reported; D2's condition): the model's PN rate changes on
    the shared glomeruli against the GH146 flies on the same glomeruli, APL on."""
    fl = flies_group('pn_dend28')
    rows = []
    for rd in ('change', 'rate', 'minsub'):
        x = sm[(sm.stage == 'pn') & (sm.cond == 'on') & (sm.subset == 'shared') & (sm.reading == rd)]
        pop = x[x.wiring == 'population'].iloc[0]
        wir = x[x.wiring != 'population']
        for g in GROUPS:
            r = differs(pop[g], wir[g].to_numpy(), fl.loc[g, 'mean'], fl.loc[g, 'lo'], fl.loc[g, 'hi'])
            rows.append(dict(element='pn', group=g, model_reading=rd, flies_reading='gh146_shared28', subset='shared',
                             share=0.0, model=float(pop[g]), flies=float(fl.loc[g, 'mean']), lo=float(fl.loc[g, 'lo']),
                             hi=float(fl.loc[g, 'hi']), n_wirings=int(len(wir)), **r))
    return pd.DataFrame(rows)


def part_c() -> dict:
    """The verdict per group (verdict.json, verdict_rows.csv): each group's elements over the readings that the noise
    check allows, its error in the effect, and its verdict; the step target's usability (step_target.json, part A); which
    of part D runs."""
    sm = pd.read_csv(out_dir() / 'stages_model.csv')
    steps = pd.read_csv(out_dir() / 'steps_model.csv')
    nc = pd.read_csv(out_dir() / 'noise_check.csv')
    passed = check_passes(nc)
    vr = verdict_rows(sm, steps, passed)
    both = {'bouton': SUBSETS, 'claw': SUBSETS}
    vb = verdict_rows(sm, steps, passed, both)
    pr = pn_verdict_rows(sm)
    pd.concat([vr.assign(verdict_subsets='primary'), vb.assign(verdict_subsets='both'), pr.assign(verdict_subsets='pn')],
              ignore_index=True).to_csv(out_dir() / 'verdict_rows.csv', index=False)
    ee = effect_errors()
    st = json.loads((out_dir() / 'step_target.json').read_text())
    out = dict(noise_check_passed=passed, effect_error=ee, groups={}, both_subsets={}, step_target_usable=bool(st['usable']),
               verdict_subsets={k: list(v) for k, v in VERDICT_SUBSETS.items()})
    pn_res = {g: element_result(pr[pr.group == g]) for g in GROUPS}
    for g in GROUPS:
        res = {el: element_result(vr[(vr.group == g) & (vr.element == el)]) for el in ('bouton', 'claw', 'step_on', 'step_off')}
        out['groups'][g] = dict(verdict=group_verdict(res, ee[g]), elements=res, pn=pn_res[g])
        rb = {el: element_result(vb[(vb.group == g) & (vb.element == el)]) for el in ('bouton', 'claw', 'step_on', 'step_off')}
        out['both_subsets'][g] = dict(verdict=group_verdict(rb, ee[g]), elements=rb)
    out['d1_due'] = bool(st['usable'])
    out['d2_due'] = bool(any(v['verdict'] == 'upstream' and v['elements']['bouton']['differs'] for v in out['groups'].values())
                         or any(v['differs'] for v in pn_res.values()))
    out['rules'] = dict(VERDICT_DIFF=VERDICT_DIFF, VERDICT_SAME_SIDE=VERDICT_SAME_SIDE, VERDICT_EFFECT_ERR=VERDICT_EFFECT_ERR,
                        NOISE_CHECK_TOL=NOISE_CHECK_TOL, soma_reading=SOMA_STEP_READING,
                        flies_steps='Yang control flies / APL output blocked (TNT) minus the claws, Pearson',
                        d2='due if a group is upstream with the bouton stage differing, or the PN stage differs for any group')
    _write_json(out_dir() / 'verdict.json', out)
    print(json.dumps({g: v['verdict'] for g, v in out['groups'].items()}, indent=1), flush=True)
    return out


# ---------------------------------------------------------------------------
# B5's readings and realism

def part_b5score() -> pd.DataFrame:
    """B5's stages read as part B's (`read_wiring`, `pn_rows`, `soma_rows`, K on its own APL-on counts) beside the real
    wiring at the 0.5 s pulse (stages_b5.csv: group means per stage, subset, share, condition and reading, for
    'real_pulse2' and 'real'), and the four bands of B5's point (b5_realism.json, E83's `realism`: E77's `realism_of`)."""
    from experiments import E83_cross_wiring_fit as E83
    st = {c: load_stage('B5', f'real_pulse2_{cc}') for c, cc in CONDS}
    if any(v is None for v in st.values()):
        raise RuntimeError('B5 has not run')
    rows = read_wiring('real_pulse2', st['on'], st['off'], unit_meta('real'), np.random.default_rng(NOISE_SEED + 100))
    for c in ('on', 'off'):
        rows += pn_rows('real_pulse2', st[c], c)
    o = {c: (st[c]['responses'], st[c]['spike_counts']) for c in ('on', 'off')}
    K = E83.readout_K(o)
    rows += soma_rows('real_pulse2', {c: st[c]['spike_counts'] for c in ('on', 'off')}, K)
    b5 = pd.DataFrame(rows)
    b5.to_parquet(local_dir() / 'b5_pairs.parquet')
    mp = pd.read_parquet(local_dir() / 'model_pairs.parquet')
    both = pd.concat([b5, mp[mp.wiring == 'real']], ignore_index=True)
    g = _group_table(both, KEY + ['wiring'])
    g.to_csv(out_dir() / 'stages_b5.csv', index=False)
    real = E83.realism(o)
    real0 = E83.realism(E83.saved_fit('real')['saved'])
    _write_json(out_dir() / 'b5_realism.json', dict(pulse2=real, pulse05_real=real0, K_pulse2=K,
                                                    K_pulse05=E83.saved_fit('real')['K']))
    print(json.dumps(real, indent=1, default=str), flush=True)
    return g


# ---------------------------------------------------------------------------
# D

D1_MODELS = ('fit1', 'P1_fc', 'P2_fc', 'M1_fc')            # E69's fit 1 and E72 part G's fits (E82's names)
D1_C2_WIRINGS = ('real', 'P1', 'P2', 'M1')                 # C2's family on the same four wirings


def d1_run(model: str) -> dict:
    """One firing-ceiling fit (E82_protocol_runs's `model_parts`: its parameters, refractory period and adaptation, its
    wiring and breadth) built as E82 builds it, run at its parameters APL on and off with APL's activity saved; its counts
    must equal its saved ones exactly; the stage matrices saved under outputs/E90_model_stages/D1/."""
    import gc
    from experiments import E82_protocol_runs as R
    m = R.model_parts(model)
    W = None if m['wiring'] is None else R.wiring_parts(m['wiring'])['W']
    sim, info = R.build_simulator(m['kc'], wPNKC=W, breadth=R.own_breadth(m), params=m['p'])
    cells = info['cells']
    mp_, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    mp_.kc.save_inh_sims = True
    Wd = info['wPNKC']
    d = local_dir('D1')
    res = {}
    for c, cc in CONDS:
        q = dict(m['p']) if c == 'on' else dict(m['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0)
        r, s_ = sim(q)
        rs, ss = m['saved'][c]
        nd = int((s_.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum())
        nr = int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum())
        res[c] = dict(counts_differ=nd, responses_differ=nr)
        if nd or nr:
            raise RuntimeError(f'D1 {model} {c}: {nd} counts and {nr} responses differ from its saved fit')
        st = read_stages(mp_, rv, Wd, order, [H.odor_abbrev(x) for x in r.columns], MODEL_WINDOW)
        _save_stage(d, f'{model}_{cc}', st, (r, s_))
        print(f'D1 {model} {c}: exact', flush=True)
    del sim, info, cells, mp_, rv
    gc.collect()
    return dict(res, wiring=m['wiring'] or 'real', source=list(m['source']), params=m['p'], kc=m['kc'])


def part_d1() -> pd.DataFrame:
    """D1: the firing-ceiling family's four fits read as C2's (identities first), their steps beside C2's family on the
    same four wirings, and both families' scores on the main target (d1_steps.csv, d1.json)."""
    from experiments import E82_protocol_runs as R
    from experiments import E83_cross_wiring_fit as E83
    t0 = time.time()
    path = out_dir() / 'd1.json'
    rec = json.loads(path.read_text()) if path.exists() else {}
    for mdl in D1_MODELS:
        if load_stage('D1', f'{mdl}_on') is None or mdl not in rec:
            rec[mdl] = d1_run(mdl)
            _write_json(path, rec)
    rows = []
    for i, mdl in enumerate(D1_MODELS):
        st = {c: load_stage('D1', f'{mdl}_{cc}') for c, cc in CONDS}
        rows += read_wiring(mdl, st['on'], st['off'], unit_meta(rec[mdl]['wiring']),
                            np.random.default_rng(NOISE_SEED + 200 + i))
        for c in ('on', 'off'):
            rows += pn_rows(mdl, st[c], c)
        o = {c: (st[c]['responses'], st[c]['spike_counts']) for c in ('on', 'off')}
        rows += soma_rows(mdl, {c: st[c]['spike_counts'] for c in ('on', 'off')}, E83.readout_K(o))
    fc = pd.DataFrame(rows)
    fc.to_parquet(local_dir() / 'd1_pairs.parquet')
    mp = pd.read_parquet(local_dir() / 'model_pairs.parquet')
    c2 = mp[mp.wiring.isin(D1_C2_WIRINGS)]
    out = []
    for fam, df in (('C2_S4', c2), ('fc_S4', fc)):
        steps = model_steps_of(df)
        pop = steps[steps.wiring == 'population']
        gt = _group_table(pop, ['step', 'soma', 'subset', 'share', 'reading'])
        gt.insert(0, 'family', fam)
        out.append(gt)
        stg = _group_table(populations(df), KEY)
        stg.insert(0, 'family', fam)
        stg.insert(1, 'step', 'stage')
        out.append(stg)
    tab = pd.concat(out, ignore_index=True)
    tab.to_csv(out_dir() / 'd1_steps.csv', index=False)
    # both families' populations of four on the main target (calibrated readout)
    tg = E83.targets()
    sc = {}
    for fam, runs in (('C2_S4', [E83.saved_fit(w)['saved'] for w in D1_C2_WIRINGS]),
                      ('fc_S4', [R.model_parts(m)['saved'] for m in D1_MODELS])):
        dist = E83.population([E83.pair_distances(o, 'cal') for o in runs])
        e = E83.effect_scores(dist, {'z8': tg['z8']})
        err = E83.pattern_error(dist, tg['z8'])
        mk = group_masks(err.index)
        sc[fam] = dict(e, **{f'err_{g}': float(err[mk[g]].mean()) for g in GROUPS})
    rec['main_target'] = sc
    rec['seconds'] = round(time.time() - t0, 1)
    _write_json(path, rec)
    _write_json(out_dir() / 'provenance_D1.json', dict(part='D1', script=H.script_path(__file__),
        model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), models=list(D1_MODELS)))
    return tab


def model_steps_of(mp: pd.DataFrame) -> pd.DataFrame:
    """`model_steps` without writing steps_model.csv."""
    claw = mp[(mp.stage == 'claw') & (mp.cond == 'on')]
    soma = mp[mp.stage == 'soma']
    rows = []
    for (w, sr, sc), s_ in soma.groupby(['wiring', 'reading', 'cond']):
        sv = pd.Series(s_.value.to_numpy(), index=pd.MultiIndex.from_arrays([s_.odor_a, s_.odor_b]))
        for (sub, sh, rd), c in claw[claw.wiring == w].groupby(['subset', 'share', 'reading']):
            cv = pd.Series(c.value.to_numpy(), index=pd.MultiIndex.from_arrays([c.odor_a, c.odor_b])).reindex(sv.index)
            st = sv - cv
            rows.append(pd.DataFrame(dict(step=sc, soma=sr, subset=sub, share=sh, reading=rd, wiring=w,
                                          odor_a=st.index.get_level_values(0), odor_b=st.index.get_level_values(1),
                                          value=st.to_numpy())))
    df = pd.concat(rows, ignore_index=True)
    pop = df.groupby(['step', 'soma', 'subset', 'share', 'reading', 'odor_a', 'odor_b'])['value'].mean().reset_index()
    pop['wiring'] = 'population'
    return pd.concat([df, pop], ignore_index=True)


D2_HANDLES = ('S1', 'S3')                 # E77's screen: the PN normalization (pn_inhadd) and the ORN exponent


def al_stage(extra_kws: Optional[dict] = None, orn: Optional[pd.DataFrame] = None) -> Dict[str, pd.DataFrame]:
    """The antennal-lobe stage alone of C2's family with an AL setting: the PN rates' window and pre-odor means
    (glomeruli x odors) of one fit_mb_model call at the bouton base's parameters (E3's; the PN layer does not depend on
    the KC layer's parameters, the wiring's claws or APL) with the bouton base's kwargs plus `extra_kws`, on the input
    `orn` (the lab's by default)."""
    from experiments import E10_joint_fit as J
    from experiments.E12_breadth_apl import e3_bouton_params
    from pyolfsysm.mb_model import fit_mb_model, get_sim_var
    orn = H.apl16_orn_deltas() if orn is None else orn
    p0 = e3_bouton_params()
    _, _, wPNKC, pd0 = fit_mb_model(orn, fixed_thr=p0['fixed_thr'], wAPLKC=p0['wAPLKC_scale'], wKCAPL=p0['wKCAPL_scale'],
                                    wAPLPN=p0['wAPLPN_scale'], wPNAPL=p0['wPNAPL_scale'],
                                    **dict(J.base_kws(None), **(extra_kws or {})), return_olfsysm_vars=True)
    mp_, rv = pd0['mp'], pd0['rv']
    t = step_times(mp_, rv)
    w = _windows(t, MODEL_WINDOW)
    pn = get_sim_var(rv.pn, 'pn_sims')
    gl = list(dict.fromkeys(wPNKC.columns.get_level_values('glomerulus')))
    out = {k: pd.DataFrame({H.odor_abbrev(c): np.asarray(pn[i])[:, w[k]].mean(axis=1) for i, c in enumerate(orn.columns)},
                           index=pd.Index(gl, name='glomerulus')) for k in ('win', 'pre')}
    return out


def pn_statistics(win: pd.DataFrame, pre: pd.DataFrame) -> dict:
    """The PN stage's statistics against A3 on the 28 shared glomeruli (rate changes): the group means of the pair
    distances and their rms difference from the PN dendrites', the CV of the per-odor totals, the correlation of the
    entries with the PN dendrites'."""
    a3 = pd.read_csv(out_dir() / 'al_transfer.csv')
    ent = a3[(a3.section == 'entries') & (a3['map'] == 'linear') & (a3.stat == 'pn')]
    P = ent.assign(g=ent.key.str.split('|').str[0], o=ent.key.str.split('|').str[1]).pivot(index='g', columns='o',
                                                                                         values='value')[list(ODORS)]
    gl = list(P.index)
    X = (win - pre).loc[gl, list(ODORS)]
    gm = group_means(odor_distances(X))
    fl = flies_group('pn_dend28')
    dif = {g: gm[g] - float(fl.loc[g, 'mean']) for g in TABLE_GROUPS}
    tot = X.sum(axis=0)
    return dict(**{f'dist_{g}': gm[g] for g in GROUPS}, **{f'diff_{g}': d for g, d in dif.items()},
                dist_rms_diff=float(np.sqrt(np.mean([d ** 2 for d in dif.values()]))), cv_totals=_cv(tot),
                cv_totals_flies=_cv(P.sum(axis=0)), r_entries_pn_dend=_corr(X.to_numpy().ravel(), P.to_numpy().ravel()),
                r_totals_pn_dend=_corr(tot, P.sum(axis=0)))


def part_d2() -> pd.DataFrame:
    """D2 (only if verdict.json says it is due): the PN stage (`al_stage`, `pn_statistics`) of E77's screen's
    normalization (S1: pn_inhadd with pn_inhsc raised as E3 / E32) and exponent (S3: the ORN input's positive deltas to a
    power) settings, identity points first (S1's must reproduce part B's PN window means exactly), beside the settings'
    saved main-target scores from E77's screen.csv; which setting is nearest the recorded stage (the rms over the five
    groups of the difference from the PN dendrites' distances; the totals' CV beside) and which the main target prefers
    (the lowest pattern part under the calibrated readout; the loss beside). d2_pn_settings.csv, d2.json."""
    from experiments import E77_input_screen as S
    v = json.loads((out_dir() / 'verdict.json').read_text())
    if not v.get('d2_due'):
        print('D2 is not due (verdict.json)', flush=True)
        return pd.DataFrame()
    t0 = time.time()
    specs = [x for x in S.point_specs() if x['handle'] in D2_HANDLES]
    screen = pd.read_csv(H.results_dir(S.EXPERIMENT) / 'screen.csv').set_index('point_id')
    ref = load_stage('B', 'real_on')
    rows = []
    for sp in specs:
        orn = S.orn_for(sp['orn'])
        st = al_stage(sp['kws'], orn)
        row = dict(point_id=sp['point_id'], handle=sp['handle'], setting=sp['setting'], value=sp['value'],
                   identity=bool(sp['identity']), **pn_statistics(st['win'], st['pre']))
        if sp['identity']:
            rw = ref['pn_win'].copy()
            rw.columns = [H.odor_abbrev(c) for c in rw.columns]
            row['identity_max_abs_diff'] = float(np.max(np.abs(st['win'][list(rw.columns)].to_numpy() - rw.to_numpy())))
        if sp['point_id'] in screen.index:
            r = screen.loc[sp['point_id']]
            for c in ('cal_loss', 'cal_pattern', 'cal_rho', 'cal_mean', 'raw_loss', 'raw_pattern', 'cal_err_ester', 'cal_err_alcohol'):
                if c in r.index:
                    row[c] = float(r[c])
        rows.append(row)
        print(f"D2 {sp['point_id']}: rms diff {row['dist_rms_diff']:.3f}, CV {row['cv_totals']:.3f}, "
              f"cal pattern {row.get('cal_pattern', np.nan):.4f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'd2_pn_settings.csv', index=False)
    res = {}
    for h, x in df.groupby('handle'):
        x = x.dropna(subset=['cal_pattern'])
        res[h] = dict(nearest_recorded=str(x.loc[x.dist_rms_diff.idxmin(), 'point_id']),
                      nearest_cv=str(x.loc[(x.cv_totals - x.cv_totals_flies).abs().idxmin(), 'point_id']),
                      preferred_by_target=str(x.loc[x.cal_pattern.idxmin(), 'point_id']),
                      preferred_by_loss=str(x.loc[x.cal_loss.idxmin(), 'point_id']))
        res[h]['same'] = res[h]['nearest_recorded'] == res[h]['preferred_by_target']
    _write_json(out_dir() / 'd2.json', dict(res, identity=[r.get('identity_max_abs_diff') for r in rows if r['identity']],
                                            seconds=round(time.time() - t0, 1)))
    _write_json(out_dir() / 'provenance_D2.json', dict(part='D2', script=H.script_path(__file__),
        model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), points=[x['point_id'] for x in specs]))
    return df


def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


TABLE_GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near')
TABLE_GROUP_HEAD = {'ester': 'ester (15)', 'alcohol': 'alcohol (6)', 'ester_alcohol': 'ester with alcohol (24)',
                    'other': 'other (21)', 'near': 'near (6)'}


def _flies_cell(fl: pd.DataFrame, g: str) -> str:
    return f"{_f(fl.loc[g, 'mean'])} [{_f(fl.loc[g, 'lo'], '.2f')}, {_f(fl.loc[g, 'hi'], '.2f')}]"


def write_variants() -> pd.DataFrame:
    """variants.csv (one row per saved run, the archive's params.csv source): tag, wiring, part, the absolute parameters,
    the protocol, the response rates and the counts' identity."""
    from experiments import E83_cross_wiring_fit as E83
    rows = []
    idn = json.loads((out_dir() / 'identity.json').read_text()) if (out_dir() / 'identity.json').exists() else {}
    for w in ALL_WIRINGS:
        st = load_stage('B', f'{w}_on')
        if st is None:
            continue
        f = E83.saved_fit(w)
        off = load_stage('B', f'{w}_acute_off')
        rows.append(dict(tag=w, variant=w, part='B', wiring=w, ref_experiment=f['folder'], ref_fit=f['tag'], k=2.0,
                         **{k: f['p'][k] for k in E83.PARAM_KEYS if k in f['p']}, kc_kc_tau=0.05, time_stim_end=0.5,
                         rate_on=float(st['responses'].to_numpy(dtype=float).mean()),
                         rate_acute_off=float(off['responses'].to_numpy(dtype=float).mean()),
                         identity_exact=bool(idn.get(f'run_{w}_on', {}).get('exact') and idn.get(f'run_{w}_off', {}).get('exact'))))
    b5 = out_dir() / 'b5_point.json'
    if b5.exists():
        b = json.loads(b5.read_text())
        st = load_stage('B5', 'real_pulse2_on')
        off = load_stage('B5', 'real_pulse2_acute_off')
        rows.append(dict(tag='real_pulse2', variant='real_pulse2', part='B5', wiring='real', ref_experiment='E82_flies_protocol',
                         ref_fit='B1_C2', k=2.0, **b['params'], kc_kc_tau=0.05, time_stim_end=PULSE, time_end=PULSE + 0.25,
                         rate_on=float(st['responses'].to_numpy(dtype=float).mean()),
                         rate_acute_off=float(off['responses'].to_numpy(dtype=float).mean()), identity_exact=bool(b.get('exact'))))
    d1 = out_dir() / 'd1.json'
    if d1.exists():
        rec = json.loads(d1.read_text())
        for mdl in D1_MODELS:
            if mdl not in rec:
                continue
            st = load_stage('D1', f'{mdl}_on')
            off = load_stage('D1', f'{mdl}_acute_off')
            r = rec[mdl]
            rows.append(dict(tag=mdl, variant=mdl, part='D1', wiring=r['wiring'], ref_experiment=r['source'][0],
                             ref_fit=r['source'][1], k=2.0, **r['params'], **r['kc'], kc_kc_tau=0.05, time_stim_end=0.5,
                             rate_on=float(st['responses'].to_numpy(dtype=float).mean()),
                             rate_acute_off=float(off['responses'].to_numpy(dtype=float).mean()),
                             identity_exact=bool(r['on']['counts_differ'] == 0 and r['off']['counts_differ'] == 0)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'variants.csv', index=False)
    return df


def comparison_table() -> None:
    """comparison_table.md: B0's identities, the stage table and the steps (flies and model, by group), the noise check,
    the verdict, the step target, B5, the PN stage and D1 if run; and the figures."""
    import subprocess
    o = out_dir()
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=H.HERE, capture_output=True, text=True).stdout.strip()
    L = [f"# E90: the model's stages against the flies' stages", '',
         f'Generated by experiments/E90_model_stages.py at {commit}. Distances are 1 - Pearson r across units unless '
         'stated; groups of the 66 megamat pairs; flies\' values with 90 % intervals from resampling each stage\'s flies '
         f'(Yang\'s recordings), {N_BOOT} resamples. The model is C2\'s family (the breadth-2 APL-effect fit, 6.8 %, the '
         '0.5 s pulse) on E72 part B\'s eight wirings at their own fits; "population" is the mean over the eight of the '
         'pair distances.', '']
    idn = json.loads((o / 'identity.json').read_text())
    ci = idn['claw_identity']
    L += ['## B0: identities', '',
          f"- Every wiring's run with the rates saved reproduces its saved fit exactly, APL on and off: "
          f"{all(v.get('exact') for k, v in idn.items() if k.startswith('run_'))} "
          f"({sum(1 for k in idn if k.startswith('run_'))} runs, counts and responses differing: "
          f"{sum(v.get('counts_differ', 0) + v.get('responses_differ', 0) for k, v in idn.items() if k.startswith('run_'))}).",
          f"- Claw inputs (wiring x bouton rate minus APL's inhibition at the step before, clipped at 0) against olfsysm's "
          f"claw_sims on the real wiring ({', '.join(ci['odors'])}; {ci['n_claws']} claws, {ci['steps_compared']} odor-steps): "
          f"{ci['claw_values_differ']} values differ (largest |difference| {ci['claw_max_abs_diff']:g}); the run's counts "
          f"equal C2's saved counts of those odors ({ci['counts_differ']} differ); one bouton per claw: {ci['one_bouton_per_claw']}.",
          '']
    fl_keys = [('orn', 'ORNs (olfdb, 9 flies, 36 glomeruli)'), ('orn_pebbled28', 'ORNs, 28 shared glomeruli'),
               ('pn_dend', 'PN dendrites (GH146, 7 flies, 29 glomeruli)'), ('pn_dend28', 'PN dendrites, 28 shared glomeruli'),
               ('bouton_measured', 'PN boutons as measured (5 flies)'), ('bouton_noise_removed', 'PN boutons, noise removed'),
               ('claw_measured', 'KC claws as measured (4 flies)'), ('claw_noise_removed', 'KC claws, noise removed'),
               ('soma_remy', "KC somata, Remy's wild type (18 flies)"), ('soma_ctrl', 'KC somata, control flies (11 recordings)'),
               ('soma_tnt', 'KC somata, APL output blocked (TNT) (8 recordings)'), ('soma_ctrl_cos', 'KC somata, control flies, cosine'),
               ('soma_tnt_cos', 'KC somata, APL output blocked (TNT), cosine')]
    L += ['## A1 and B: the stage table by group', '', '### The flies', '',
          '| stage | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) + ' |',
          '|---|' + '---:|' * len(TABLE_GROUPS)]
    for k, lab in fl_keys:
        fl = flies_group(k)
        L.append(f'| {lab} | ' + ' | '.join(_flies_cell(fl, g) for g in TABLE_GROUPS) + ' |')
    sm = pd.read_csv(o / 'stages_model.csv')
    pop = sm[sm.wiring == 'population']

    def mrow(stage, subset, share, cond, reading):
        x = pop[(pop.stage == stage) & (pop.subset == subset) & (pop.share == share) & (pop.cond == cond) & (pop.reading == reading)]
        w = sm[(sm.wiring != 'population') & (sm.stage == stage) & (sm.subset == subset) & (sm.share == share)
               & (sm.cond == cond) & (sm.reading == reading)]
        if not len(x):
            return None
        return ' | '.join(f"{_f(x.iloc[0][g])} ({_f(w[g].min(), '.2f')}-{_f(w[g].max(), '.2f')})" for g in TABLE_GROUPS)
    L += ['', '### The model (population of eight; the wirings\' range in parentheses)', '',
          '| stage | subset | share of the spike signal | APL | reading | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) + ' |',
          '|---|---|---:|---|---|' + '---:|' * len(TABLE_GROUPS)]
    for stage, subsets, shares, readings in (('pn', ('shared', 'all'), (0.0,), ('change', 'rate', 'minsub')),
                                             ('bouton', ('shared', 'all'), (0.0,), MODEL_READINGS),
                                             ('claw', ('all', 'shared'), CLAW_SPIKE_SHARES, MODEL_READINGS),
                                             ('soma', ('all',), (0.0,), ('cal_pearson', 'cal_cosine', 'raw_pearson', 'raw_cosine'))):
        for subset in subsets:
            for share in shares:
                for cond in ('on', 'off'):
                    for rd in readings:
                        r = mrow(stage, subset, share, cond, rd)
                        if r is not None:
                            L.append(f'| {stage} | {subset} | {share:g} | {cond} | {rd} | {r} |')
    # steps
    sf = pd.read_csv(o / 'steps_flies_summary.csv')
    L += ['', '## A2: the flies\' steps (somata minus claws)', '',
          '| step | claws | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) +
          ' | ester minus alcohol [90 %] | r with the APL effect | r with the common error | halves agree | floor | SD over pairs |',
          '|---|---|' + '---:|' * (len(TABLE_GROUPS) + 6)]
    stp = pd.read_csv(o / 'steps_flies.csv')
    for _, r in sf.iterrows():
        x = stp[(stp.step == r.step) & (stp.claw == r.claw) & (stp.kind == 'group')].set_index('group')
        L.append(f'| {r.step} | {r.claw[5:]} | ' + ' | '.join(_flies_cell(x, g) for g in TABLE_GROUPS) +
                 f" | {_f(r.ester_minus_alcohol)} [{_f(r.ester_minus_alcohol_lo, '.2f')}, {_f(r.ester_minus_alcohol_hi, '.2f')}]"
                 f' | {_f(r.r_effect)} | {_f(r.r_common_error)} | {_f(r.halves_common)} | {_f(r.floor)} | {_f(r.sd_over_pairs)} |')
    L += ['', 'on_remy: Remy\'s wild-type somata minus the claws; on_ctrl / off_tnt: Yang\'s control flies / APL output '
          'blocked (TNT) minus the claws (Pearson; _cos: the somata as cosine). Halves: Pearson r over the 66 pairs between '
          'the steps from the two halves of the claw flies (three splits of two and two), all of the somata\'s flies in '
          'both; floor: the rms over the 66 pairs (means removed) of a resampled step minus the measured one.', '']
    st = json.loads((o / 'step_target.json').read_text())
    L += ['## The step target', '', f"Rule: {st['rule']}.", '']
    for k, v in st['readings'].items():
        L.append(f"- claws {k}: halves agree {v['halves_agreement']:.3f} (somata also split: {v['halves_agreement_split_somata']:.3f}), "
                 f"floor {v['floor']:.4f} against half its SD {0.5 * v['sd_over_pairs']:.4f}: usable {v['usable']}.")
    L += [f"- Usable: {st['usable']}.", '']
    # model steps
    sg = _group_table(pd.read_csv(o / 'steps_model.csv'), ['step', 'soma', 'subset', 'share', 'reading', 'wiring'])
    L += ['## B: the model\'s steps (somata minus its claw stage with APL on)', '',
          '| step | somata | claws (subset, share, reading) | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) + ' |',
          '|---|---|---|' + '---:|' * len(TABLE_GROUPS)]
    for step in ('on', 'off'):
        for soma in (SOMA_STEP_READING,) + SOMA_STEP_BESIDE:
            for sub, sh, rd in (('all', 0.0, 'session_z'), ('all', 0.0, 'rate'), ('all', 0.0, 'noisy_nr'), ('all', 0.25, 'session_z'),
                                ('all', 0.5, 'session_z'), ('shared', 0.0, 'session_z')):
                if soma != SOMA_STEP_READING and (sub, sh, rd) != ('all', 0.0, 'session_z'):
                    continue
                x = sg[(sg.step == step) & (sg.soma == soma) & (sg.subset == sub) & (sg.share == sh) & (sg.reading == rd)]
                pw, ww = x[x.wiring == 'population'], x[x.wiring != 'population']
                if not len(pw):
                    continue
                L.append(f'| {step} | {soma} | {sub}, {sh:g}, {rd} | ' + ' | '.join(
                    f"{_f(pw.iloc[0][g])} ({_f(ww[g].min(), '.2f')}-{_f(ww[g].max(), '.2f')})" for g in TABLE_GROUPS) + ' |')
    # noise check
    nc = pd.read_csv(o / 'noise_check.csv')
    pnc = nc[nc.wiring == 'population']
    L += ['', '## The check of the noise removal', '',
          f'On the model\'s known values (the session z reading), three trials per odor at the flies\' per-odor block-to-block '
          f'reliability (boutons {remy_reliability("bouton")[list(ODORS)].to_numpy().mean():.3f}, claws '
          f'{remy_reliability("claw")[list(ODORS)].to_numpy().mean():.3f} on average), as many pseudo-flies as the stage has '
          f'flies, {N_NOISE_REPS} replicates; tolerance {NOISE_CHECK_TOL}.', '',
          '| stage | largest abs error of the noise-removed estimate (population, every group, subset, share, APL on / off) | '
          'largest over single wirings | as measured: mean error (attenuation) | passes |', '|---|---:|---:|---:|---|']
    for stg in ('bouton', 'claw'):
        a, b = pnc[pnc.stage == stg], nc[(nc.stage == stg) & (nc.wiring != 'population')]
        L.append(f'| {stg} | {_f(a.error.abs().max(), ".4f")} | {_f(b.error.abs().max(), ".4f")} | '
                 f'{_f(a.error_measured.mean(), ".3f")} | {bool(a.within.all())} |')
    # verdict
    v = json.loads((o / 'verdict.json').read_text())
    vr = pd.read_csv(o / 'verdict_rows.csv')
    L += ['', '## C: the verdict', '',
          f"Readings: the noise check {v['noise_check_passed']}; subsets {v['verdict_subsets']} (the steps on the claws' "
          f"subset); the model's steps with the somata under the calibrated readout (Pearson) against Yang's control flies "
          f"and APL output blocked (TNT) minus the claws.", '',
          '| group | verdict | boutons | claws | step with APL on | step without APL | PN stage (beside) | error in the effect '
          '| verdict with both subsets as readings |', '|---|---|---|---|---|---|---|---:|---|']

    def el(e):
        return f"{e['n_differ']}/{e['n_readings']}" + (' differs' if e['differs'] else ' mixed' if e['mixed'] else '')
    for g in GROUPS:
        x = v['groups'][g]
        L.append(f"| {GROUP_NAMES[g]} | {x['verdict']} | {el(x['elements']['bouton'])} | {el(x['elements']['claw'])} | "
                 f"{el(x['elements']['step_on'])} | {el(x['elements']['step_off'])} | {el(x['pn'])} | "
                 f"{_f(v['effect_error'][g], '+.3f')} | {v['both_subsets'][g]['verdict']} |")
    L += ['', 'Cells: readings in which the element differs (all of the rule\'s three conditions) / readings read.', '',
          '### The comparisons behind it (primary subsets; model population against the flies)', '',
          '| element | group | model reading | flies | share | model | flies [90 %] | wirings on the side | differs |',
          '|---|---|---|---|---:|---:|---:|---:|---|']
    vp = vr[vr.verdict_subsets.isin(['primary', 'pn'])]
    for _, r in vp[vp.group.isin(TABLE_GROUPS)].iterrows():
        L.append(f"| {r.element} | {r.group} | {r.model_reading} | {r.flies_reading} | {r.share:g} | {_f(r.model)} | "
                 f"{_f(r.flies)} [{_f(r.lo, '.2f')}, {_f(r.hi, '.2f')}] | {int(r.n_same_side)}/{int(r.n_wirings)} | {bool(r.differs)} |")
    # PN stage
    ps = pd.read_csv(o / 'pn_stage.csv')
    L += ['', '## B1: the PN stage on the 28 shared glomeruli (outcome g)', '',
          '| reading | r of the entries with the PN dendrites | with the ORN dF/F | with the ORN input | CV of the per-odor totals | '
          'r of the totals with the PN dendrites\' | ' + ' | '.join(f'distance, {TABLE_GROUP_HEAD[g]}' for g in TABLE_GROUPS) + ' |',
          '|---|' + '---:|' * (5 + len(TABLE_GROUPS))]
    for _, r in ps.iterrows():
        L.append(f"| {r.reading} | {_f(r.get('r_entries_pn_dend'))} | {_f(r.get('r_entries_orn_dff'))} | "
                 f"{_f(r.get('r_entries_orn_input'))} | {_f(r.cv_totals)} | {_f(r.get('r_totals_pn_dend'))} | " +
                 ' | '.join(_f(r[f'dist_{g}']) for g in TABLE_GROUPS) + ' |')
    a3 = pd.read_csv(o / 'al_transfer.csv')
    sm3 = a3[a3.section == 'summary']
    L += ['', 'A3 (the flies, linear dF/F and through the saturating map, the bound for the indicator): ' + '; '.join(
        f"{r['map']} {r.stat} {_f(r.value)} [{_f(r.lo, '.2f')}, {_f(r.hi, '.2f')}]" for _, r in sm3.iterrows()) + '.', '']
    # B5
    if (o / 'stages_b5.csv').exists():
        b5 = pd.read_csv(o / 'stages_b5.csv')
        rz = json.loads((o / 'b5_realism.json').read_text())
        L += ['## B5: the flies\' 2 s pulse (E82\'s re-held B1_C2 point), stages over 0.5-2 s, against the 0.5 s pulse (real wiring)', '',
              '| stage | subset | APL | reading | pulse | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) + ' |',
              '|---|---|---|---|---|' + '---:|' * len(TABLE_GROUPS)]
        for stage, sub, rd in (('pn', 'shared', 'change'), ('bouton', 'shared', 'session_z'), ('bouton', 'shared', 'rate'),
                               ('claw', 'all', 'session_z'), ('claw', 'all', 'rate'), ('soma', 'all', 'cal_pearson'),
                               ('soma', 'all', 'raw_cosine')):
            for cond in ('on', 'off'):
                for w, lab in (('real', '0.5 s'), ('real_pulse2', '2 s')):
                    x = b5[(b5.wiring == w) & (b5.stage == stage) & (b5.subset == sub) & (b5.share == 0.0) & (b5.cond == cond) & (b5.reading == rd)]
                    if len(x):
                        L.append(f'| {stage} | {sub} | {cond} | {rd} | {lab} | ' + ' | '.join(_f(x.iloc[0][g]) for g in TABLE_GROUPS) + ' |')
        L += ['', '### Realism of B5\'s point (the four bands, the rule of 2026-11-10 16:10; the earlier rule beside)', '',
              '| statistic | 2 s pulse (B5) | 0.5 s pulse (C2, real wiring) | band |', '|---|---:|---:|---|',
              f"| APL-on response rate | {_f(rz['pulse2']['rate_on'], '.4f')} | {_f(rz['pulse05_real']['rate_on'], '.4f')} | 0.045-0.085 |",
              f"| APL-off / APL-on ratio | {_f(rz['pulse2']['ratio'])} | {_f(rz['pulse05_real']['ratio'])} | 1.6-2.9 |",
              f"| response rate ρ, APL on | {_f(rz['pulse2']['rate_rho_on'])} | {_f(rz['pulse05_real']['rate_rho_on'])} | 0.37 or more |",
              f"| odors per KC, distance, APL on | {_f(rz['pulse2']['opk_on_tv'])} | {_f(rz['pulse05_real']['opk_on_tv'])} | 0.155 |",
              f"| odors per KC, distance, APL off | {_f(rz['pulse2']['opk_off_tv'])} | {_f(rz['pulse05_real']['opk_off_tv'])} | 0.223 |",
              f"| largest spike count, APL on / off | {rz['pulse2']['max_count_on']} / {rz['pulse2']['max_count_off']} | "
              f"{rz['pulse05_real']['max_count_on']} / {rz['pulse05_real']['max_count_off']} | |",
              f"| flag | {rz['pulse2']['realism_flag'] or 'none'} | {rz['pulse05_real']['realism_flag'] or 'none'} | |",
              f"| flag, earlier rule (Yang's z) | {rz['pulse2']['realism_flag_zcrit'] or 'none'} | {rz['pulse05_real']['realism_flag_zcrit'] or 'none'} | |", '']
    if (o / 'd1_steps.csv').exists():
        d1 = pd.read_csv(o / 'd1_steps.csv')
        rec = json.loads((o / 'd1.json').read_text())
        L += ['## D1: the firing-ceiling family (fit 1, P1, P2, M1) beside C2\'s on the same four wirings', '',
              '| family | element | somata | claws reading | ' + ' | '.join(TABLE_GROUP_HEAD[g] for g in TABLE_GROUPS) + ' |',
              '|---|---|---|---|' + '---:|' * len(TABLE_GROUPS)]
        for fam in ('C2_S4', 'fc_S4'):
            for step in ('on', 'off'):
                x = d1[(d1.family == fam) & (d1.step == step) & (d1.soma == SOMA_STEP_READING) & (d1.subset == 'all')
                       & (d1.share == 0.0) & (d1.reading == 'session_z')]
                if len(x):
                    L.append(f'| {fam} | step {step} | {SOMA_STEP_READING} | all, session_z | ' + ' | '.join(_f(x.iloc[0][g]) for g in TABLE_GROUPS) + ' |')
            for stage, sub in (('bouton', 'shared'), ('claw', 'all')):
                x = d1[(d1.family == fam) & (d1.step == 'stage') & (d1.stage == stage) & (d1.subset == sub) & (d1.share == 0.0)
                       & (d1.cond == 'on') & (d1.reading == 'session_z')]
                if len(x):
                    L.append(f'| {fam} | {stage} (APL on) | | {sub}, session_z | ' + ' | '.join(_f(x.iloc[0][g]) for g in TABLE_GROUPS) + ' |')
        mt = rec.get('main_target', {})
        if mt:
            L += ['', '| family (population of four) | APL-effect loss (z8) | pattern part | ρ on the per-pair APL effect | mean APL effect | '
                  'error, ester | error, alcohol |', '|---|---:|---:|---:|---:|---:|---:|']
            for fam, r in mt.items():
                L.append(f"| {fam} | {_f(r['z8_loss'], '.4f')} | {_f(r['z8_pattern'], '.4f')} | {_f(r['z8_rho'])} | "
                         f"{_f(r['mean_effect'], '+.3f')} | {_f(r['err_ester'], '+.3f')} | {_f(r['err_alcohol'], '+.3f')} |")
        L.append('')
    if (o / 'd2_pn_settings.csv').exists():
        d2 = pd.read_csv(o / 'd2_pn_settings.csv')
        dj = json.loads((o / 'd2.json').read_text())
        L += ['## D2: the PN stage of E77\'s screen\'s normalization (S1) and exponent (S3) settings beside their main-target '
              'scores', '', '| point | rms difference from the PN dendrites\' group distances | ' +
              ' | '.join(f'distance, {TABLE_GROUP_HEAD[g]}' for g in TABLE_GROUPS) +
              ' | CV of the totals (PN dendrites ' + _f(d2.cv_totals_flies.iloc[0]) + ') | r of the entries | pattern part (cal) '
              '| APL-effect loss (cal) | ρ on the per-pair APL effect (cal) |', '|---|' + '---:|' * (len(TABLE_GROUPS) + 6)]
        for _, r in d2.iterrows():
            L.append(f'| {r.point_id} | {_f(r.dist_rms_diff)} | ' + ' | '.join(_f(r[f"dist_{g}"]) for g in TABLE_GROUPS) +
                     f' | {_f(r.cv_totals)} | {_f(r.r_entries_pn_dend)} | {_f(r.get("cal_pattern"), ".4f")} | '
                     f'{_f(r.get("cal_loss"), ".4f")} | {_f(r.get("cal_rho"))} |')
        L += [''] + [f"- {h}: nearest the recorded stage {x['nearest_recorded']} (nearest on the totals' CV {x['nearest_cv']}); "
                     f"preferred by the main target's pattern part {x['preferred_by_target']} (by the loss {x['preferred_by_loss']}); "
                     f"the same: {x['same']}." for h, x in dj.items() if h in D2_HANDLES] + ['']
    (o / 'comparison_table.md').write_text('\n'.join(L) + '\n')
    write_variants()
    figures()


def figures() -> None:
    """E90_stages.png (the group means by stage, flies with intervals and the model's population with its wirings, ester
    and alcohol groups and the other two) and E90_steps.png (per pair: the flies' step in wild-type flies against the
    APL effect; the model's steps against the flies')."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    o = out_dir()
    sm = pd.read_csv(o / 'stages_model.csv')
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.6), sharey=True)
    stages = [('ORNs', 'orn', ('pn', 'shared', 0.0, 'on', 'change'), None),
              ('PN dendrites', 'pn_dend28', ('pn', 'shared', 0.0, 'on', 'change'), 'model PNs'),
              ('PN boutons', 'bouton_noise_removed', ('bouton', 'shared', 0.0, 'on', 'session_z'), None),
              ('KC claws', 'claw_noise_removed', ('claw', 'all', 0.0, 'on', 'session_z'), None),
              ('somata, APL on', 'soma_ctrl', ('soma', 'all', 0.0, 'on', 'cal_pearson'), None),
              ('somata, APL off', 'soma_tnt', ('soma', 'all', 0.0, 'off', 'cal_pearson'), None)]
    for ax, g in zip(axes, ('ester', 'alcohol', 'ester_alcohol', 'other')):
        x = np.arange(len(stages))
        fm = [flies_group(k).loc[g] for _, k, _, _ in stages]
        ax.errorbar(x - 0.08, [f['mean'] for f in fm], yerr=[[f['mean'] - f['lo'] for f in fm], [f['hi'] - f['mean'] for f in fm]],
                    fmt='o-', color='k', lw=1.4, capsize=3, label='flies (90 % interval)')
        for i, (_, _, (st, sub, sh, cond, rd), _) in enumerate(stages):
            if i == 0:
                continue
            sel = sm[(sm.stage == st) & (sm.subset == sub) & (sm.share == sh) & (sm.cond == cond) & (sm.reading == rd)]
            w = sel[sel.wiring != 'population'][g].to_numpy()
            ax.plot(np.full(len(w), i + 0.08), w, '.', color='#d9480f', alpha=0.5, ms=5)
        mv = [np.nan] + [float(sm[(sm.stage == st) & (sm.subset == sub) & (sm.share == sh) & (sm.cond == cond) &
                                  (sm.reading == rd) & (sm.wiring == 'population')][g].iloc[0]) for _, _, (st, sub, sh, cond, rd), _ in stages[1:]]
        ax.plot(x + 0.08, mv, 's--', color='#d9480f', lw=1.2, ms=5, label='model, population of 8 (dots: wirings)')
        ax.set_xticks(x)
        ax.set_xticklabels([s[0] for s in stages], rotation=35, ha='right', fontsize=8)
        ax.set_title(GROUP_NAMES[g], fontsize=10)
        ax.grid(True, color='0.9', lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('distance between two odors\n(1 - Pearson r across units), group mean')
    axes[0].legend(fontsize=7.5, frameon=False, loc='lower right')
    fig.suptitle('Stages by group: flies (ORNs and PN dendrites per fly; boutons and claws noise removed; somata: control '
                 'flies / APL output blocked (TNT)) and the model (PN rate changes, session z of boutons and claws, somata '
                 'under the calibrated readout)', fontsize=9)
    fig.tight_layout()
    fig.savefig(o / 'E90_stages.png', dpi=150)
    plt.close(fig)
    # the steps per pair
    S = stage_samples()
    stp = pd.read_csv(o / 'steps_flies.csv')
    pr = stp[stp.kind == 'pair']
    eff = _target_effect()
    ms = pd.read_csv(o / 'steps_model.csv')
    ms = ms[(ms.wiring == 'population') & (ms.soma == SOMA_STEP_READING) & (ms.subset == 'all') & (ms.share == 0.0)
            & (ms.reading == 'session_z')]
    masks = group_masks(megamat_pairs())
    col = np.where(masks['ester'], '#2b6cb0', np.where(masks['alcohol'], '#d9480f', np.where(masks['ester_alcohol'], '#6b46c1', '0.6')))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    def series(df, step, claw):
        x = df[(df.step == step) & (df.claw == claw)]
        return pd.Series(x['mean'].to_numpy(), index=pd.MultiIndex.from_arrays([x.odor_a, x.odor_b])).reindex(megamat_pairs())
    on_w = series(pr, 'on_remy', 'claw_noise_removed')
    ax = axes[0]
    ax.scatter(eff.to_numpy(), on_w.to_numpy(), c=col, s=18)
    ax.set_xlabel("per-pair APL effect, KCs (control flies minus APL output blocked (TNT)), cosine")
    ax.set_ylabel("step from claws to somata, wild type (Remy's somata minus claws, noise removed)")
    ax.set_title(f'flies: r = {_corr(eff, on_w):.2f}', fontsize=10)
    for ax, (step, flstep) in zip(axes[1:], (('on', 'on_ctrl'), ('off', 'off_tnt'))):
        fv = series(pr, flstep, 'claw_noise_removed')
        mv = ms[ms.step == step]
        mv = pd.Series(mv.value.to_numpy(), index=pd.MultiIndex.from_arrays([mv.odor_a, mv.odor_b])).reindex(megamat_pairs())
        ax.scatter(fv.to_numpy(), mv.to_numpy(), c=col, s=18)
        lim = [min(np.nanmin(fv), np.nanmin(mv)) - 0.05, max(np.nanmax(fv), np.nanmax(mv)) + 0.05]
        ax.plot(lim, lim, ':', color='0.5', lw=1)
        ax.set_xlabel(f"flies: {'control flies' if step == 'on' else 'APL output blocked (TNT)'} minus claws (noise removed)")
        ax.set_ylabel(f"model: somata APL {step} (calibrated readout) minus claws (session z), population")
        ax.set_title(f"step {'with APL on' if step == 'on' else 'without APL'}: r = {_corr(fv, mv):.2f}", fontsize=10)
    for ax in axes:
        ax.grid(True, color='0.9', lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].scatter([], [], c='#2b6cb0', label='ester group')
    axes[0].scatter([], [], c='#d9480f', label='alcohol group')
    axes[0].scatter([], [], c='#6b46c1', label='ester with alcohol')
    axes[0].scatter([], [], c='0.6', label='other pairs')
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.suptitle('The step from claws to somata, per megamat pair (1 - Pearson r)', fontsize=10)
    fig.tight_layout()
    fig.savefig(o / 'E90_steps.png', dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'B0', 'B', 'B5', 'B5score', 'Bscore', 'C', 'D1', 'D2', 'table'])
    ap.add_argument('--wirings', nargs='*', default=None)
    a = ap.parse_args()
    {'A': part_a, 'B0': part_b0, 'B': lambda: part_b(a.wirings), 'B5': part_b5, 'B5score': part_b5score,
     'Bscore': part_bscore, 'C': part_c,
     'D1': part_d1, 'D2': part_d2, 'table': comparison_table}[a.part]()


if __name__ == '__main__':
    main()
