#!/usr/bin/env python3
"""E119 part 0: Remy's per-ROI megamat recordings (PN boutons, KC claws) reduced to ROI x trial window means, each ROI's
pre-odor noise SD and a lock file, with the checks that they reproduce Remy's own per-fly distance matrices and the
decoded dashboard matrices.

Reads the byte copies in <shared>/data/remy_per_roi/source/ (made by `--copy` at viki's user's instruction of 2026-11-20,
Remy's data shared with her OK; source/README.md lists every file with its sha256) of Remy's analysis outputs on tensor,
<REMY_ROOT>/data_by_acq/ (read-only; read directly if the copy is absent):
- kc_claws/volumetric_with_neuropil_correction/<acq>/trials/xrds_bc_trials.nc (4 flies) and
  pn_boutons/volumetric_without_neuropil_correction/<acq>/trials/xrds_bc_trials.nc (6 flies, of which 2022-10-21 fly02
  is excluded from every analysis, EXCLUDED, as Remy's dashboards leave it out; checks 0 and 1 run on it): trials (17 odors x 3) x
  ROIs x time (-5 to 19.8 s, 0.2 s frames), each trial's baseline (the median of -5 to 0 s) subtracted, variables F,
  Fneu, Fc, spks, F_zscore and Fc_zscore;
- <acq>/respvec/xrds_mean_peak.nc: Remy's own per-trial window means (0.5 to 2.5 s), check 0;
- <stage folder>/aggregate/xrds_stim_rdm_concat.nc: Remy's per-fly 17 x 17 correlation distances (1 - Pearson r across
  ROIs of the trial-mean window means), check 1;
- <shared>/data/remy_rdm_timeseries/ (decoded from Remy's dashboards): the per-fly per-frame matrices, check 2.

Writes <shared>/data/remy_per_roi/ (README.md there):
- window_means.parquet: one row per (stage, acq, roi, trial) with the trial's odor and repeat, and per variable
  (VARIABLES) and window (WINDOWS) the mean over the window's frames (closed windows, as Remy's respvec and the
  dashboards' window_csvs.py select frames), columns '<variable>__<window>';
- noise_sd.parquet: per (stage, acq, roi) and variable, the SD of the frames within the 51 pre-odor windows (PRE_NOISE),
  each window's mean removed (the baseline-noise unit of the realism bands' data side,
  notes/response_criterion_decision.md);
- frame_rdms.npz: per stage and acq the per-frame trial-mean correlation distances (frames x 17 x 17) of the variable
  Remy's figures use (DIST_VAR), with the frame times;
- lock.json: the files read (path, size, mtime, sha256), the flies, ROI counts, frame times, odor names and their trial
  order, the windows and the checks' results; checks.csv; provenance.json.

The script stops (exit 1, nothing but the checks written) if a check fails.

Part 0b (`--part0b`; the same fly, 2023-04-11 fly 3, megamat1_redo): its claws and somata cut from the copy of Remy's
suite2p output (source/same_fly/; copied with `--copy` from matrix at viki's user's instruction of 2026-11-20), gated on
reproducing the decoded same-fly series (data/remy_rdm_timeseries/same_fly/) at r >= 0.98 and Remy's claws-against-somata
agreement of 0.86; writes same_fly_checks.json and, if the gate passes, same_fly_trials.npz.

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u data_prep/remy_per_roi_data.py --copy      # the copies, once
    nice -n 10 ./venv/bin/python -u data_prep/remy_per_roi_data.py [--out DIR] [--no-sha]
    nice -n 10 ./venv/bin/python -u data_prep/remy_per_roi_data.py --part0b
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

REMY_ROOT = Path('/local/tensor/Remy/odor_space_collab/analysis_outputs')
# the same fly's suite2p output and acquisition files (part 0b), read-only
SAME_FLY_SRC = Path('/local/matrix/Remy-Data/projects/odor_space_collab/processed_data/2023-04-11/3/megamat1_redo')
SAME_FLY_NAME = '2023-04-11__fly03__megamat1_redo'
SAME_FLY_COMBINED = ('F.npy', 'Fneu.npy', 'stat.npy', 'iscell.npy', 'iscell_calyx.npy', 'iscell_soma.npy', 'ops.npy',
                     'ops.yaml')
SAME_FLY_ACQ_FILES = ('timestamps.npy', 'timestamps0.npy', 'stim_list.json', 'pin_odor_list.json',
                      'thorimage_settings.yaml', 'Experiment.xml')
# the bouton fly Remy's dashboards and bouton analysis leave out (it matches no dashboard fly; its median cellprob is
# 0.03): copied, never analyzed (viki's user, 2026-11-20)
EXCLUDED = {('bouton', '2022-10-21__fly02__megamat1_001__dsub3'):
            "left out of Remy's bouton analysis (matches no dashboard fly); not analyzed (viki's user, 2026-11-20)"}
STAGE_DIRS = {'bouton': 'pn_boutons/volumetric_without_neuropil_correction',
              'claw': 'kc_claws/volumetric_with_neuropil_correction'}
VARIABLES = ('F', 'Fc', 'F_zscore', 'Fc_zscore', 'spks')
# the variable Remy's figures and the dashboards use per stage (E90's REMY_VARIANT)
DIST_VAR = {'bouton': 'F_zscore', 'claw': 'Fc_zscore'}
# the variable read in baseline-noise units per stage (the boutons were not neuropil-corrected)
NOISE_VAR = {'bouton': 'F', 'claw': 'Fc'}
# the windows (s after odor onset), closed: Remy's respvec window, E90's STAGE_WINDOW, the dashboards' early / mid / late
# windows and the pre-odor window (the onset frame left out)
WINDOWS = {'remy': (0.5, 2.5), 'e90': (0.5, 3.0), 'early': (0.5, 1.0), 'mid': (1.0, 2.0), 'late': (2.0, 3.0),
           'pre': (-5.0, -0.2), 'base': (-5.0, 0.0)}
# Remy's respvec subtracts the mean over 'base' (-5 to 0 s, closed) where its file says so (respvec.baseline_method
# 'mean': the boutons), and nothing otherwise (the claws); her stim RDMs are made from those respvec values
PRE_NOISE = (-5.0, -0.2)                # the pre-odor frames whose SD (each trial's mean removed) is the noise unit
EXACT_TOL = 1e-6                        # checks 0 and 1: Remy's window means and stim RDMs reproduced to this
DASH_WINDOW = (0.5, 3.0)                # check 2: the window over which the per-frame matrices are averaged
DASH_R_MIN = 0.98                       # check 2: per fly, Pearson r over the 136 pairs with the dashboard's
SAME_FLY_R_MIN = 0.98                   # part 0b's gate (the entry)
SAME_FLY_AGREEMENT = 0.86               # Remy's claws-against-somata agreement (1-2 s), part 0b's gate
SAME_FLY_AGREEMENT_TOL = 0.02


def out_dir() -> Path:
    return H.shared_dir() / 'data' / 'remy_per_roi'


def odor_of_stim(stim: str) -> str:
    """'ms @ -3.0' -> 'ms'."""
    return str(stim).split('@')[0].strip()


def acquisitions(root: Path = REMY_ROOT) -> pd.DataFrame:
    """One row per stage and acquisition (stage, acq: the folder name, path of xrds_bc_trials.nc), sorted."""
    rows = []
    for stage, sub in STAGE_DIRS.items():
        for p in sorted((root / 'data_by_acq' / sub).glob('*/trials/xrds_bc_trials.nc')):
            rows.append(dict(stage=stage, acq=p.parent.parent.name, path=p))
    return pd.DataFrame(rows)


def frame_mask(t: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    """The frames of a closed window."""
    t = np.asarray(t, dtype=float)
    return (t >= window[0] - 1e-6) & (t <= window[1] + 1e-6)


def window_means(x: np.ndarray, t: np.ndarray, windows: Dict[str, Tuple[float, float]] = WINDOWS
                 ) -> Dict[str, np.ndarray]:
    """{window: trials x ROIs} means of x (trials x ROIs x frames) over each window's frames."""
    x = np.asarray(x, dtype=float)
    return {k: x[..., frame_mask(t, w)].mean(axis=-1) for k, w in windows.items()}


def noise_sd(x: np.ndarray, t: np.ndarray, window: Tuple[float, float] = PRE_NOISE) -> np.ndarray:
    """Per ROI the SD of the pre-odor frames of all trials, each trial's pre-odor mean removed (ROIs; ddof 1 per trial's
    removed mean)."""
    p = np.asarray(x, dtype=float)[..., frame_mask(t, window)]          # trials x ROIs x frames
    d = p - p.mean(axis=-1, keepdims=True)
    n_tr, _, n_fr = p.shape
    return np.sqrt((d ** 2).sum(axis=(0, 2)) / (n_tr * (n_fr - 1)))


def trial_repeats(stims: Sequence[str]) -> np.ndarray:
    """The repeat (0, 1, 2, ...) of every trial: its stimulus's occurrence number in the trial order."""
    seen: Dict[str, int] = {}
    out = []
    for s in stims:
        out.append(seen.get(s, 0))
        seen[s] = seen.get(s, 0) + 1
    return np.asarray(out, dtype=int)


def _pearson_dist(M: np.ndarray) -> np.ndarray:
    """1 - Pearson r between the rows of M (odors x ROIs)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        return 1.0 - np.corrcoef(M)


def odor_means(wm: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """The trial means per odor of trial-first values (trials x ...) -> (odors x ...)."""
    o = np.asarray(odors_of_trials)
    return np.stack([np.asarray(wm, dtype=float)[o == k].mean(axis=0) for k in odors])


def stim_rdm(wm: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """1 - Pearson r across ROIs of the trial-mean window means (trials x ROIs) of every two odors (odors x odors)."""
    return _pearson_dist(odor_means(wm, odors_of_trials, odors))


def frame_rdms(x: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """Per frame, 1 - Pearson r across ROIs of the trial-mean responses (x trials x ROIs x frames -> frames x odors x
    odors)."""
    m = odor_means(x, odors_of_trials, odors)                          # odors x ROIs x frames
    z = m - m.mean(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = z / np.linalg.norm(z, axis=1, keepdims=True)
    return 1.0 - np.einsum('aut,but->tab', z, z)


def _sha256(path: Path, block: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(block), b''):
            h.update(b)
    return h.hexdigest()


def file_entry(path: Path, sha: bool = True) -> dict:
    st = Path(path).stat()
    return dict(path=str(path), size=int(st.st_size), mtime=datetime.fromtimestamp(st.st_mtime).isoformat(timespec='seconds'),
                sha256=_sha256(path) if sha else None)


def reduce_acq(stage: str, acq: str, path: Path, sha: bool = True) -> dict:
    """One acquisition's window means (long DataFrame), noise SDs, frame RDMs (DIST_VAR), its stim RDMs per variable
    (Remy's window), Remy's own window means checked (check 0), and its lock entry."""
    import xarray as xr
    path = Path(path)
    ds = xr.open_dataset(path)
    try:
        t = np.round(ds['time'].to_numpy().astype(float), 6)
        stims = [str(s) for s in ds['stim'].to_numpy()]
        odors_t = [odor_of_stim(s) for s in stims]
        odors = sorted(set(odors_t))
        rep = trial_repeats(stims)
        if 'stim_occ' in ds and not np.array_equal(rep, ds['stim_occ'].to_numpy().astype(int)):
            raise RuntimeError(f'{acq}: the repeats are not the file\'s stim_occ')
        rois = ds['cells'].to_numpy().astype(int)
        rp = path.parent.parent / 'respvec' / 'xrds_mean_peak.nc'
        remy_base = False
        if rp.exists():
            with xr.open_dataset(rp) as mp0:
                remy_base = str(mp0.attrs.get('respvec.baseline_method', '')) == 'mean'
        cols, noise, rdms, remy_resp = {}, {}, {}, {}
        frames = None
        for v in VARIABLES:
            x = ds[v].to_numpy().astype(float)                          # trials x ROIs x frames
            wm_v = window_means(x, t)
            for w, a in wm_v.items():
                cols[f'{v}__{w}'] = a.ravel()
            noise[v] = noise_sd(x, t)
            remy_resp[v] = wm_v['remy'] - (wm_v['base'] if remy_base else 0.0)
            rdms[v] = stim_rdm(remy_resp[v], odors_t, odors)
            if v == DIST_VAR[stage]:
                frames = frame_rdms(x, odors_t, odors)
            del x
        attrs = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in ds.attrs.items()}
        n_tr = len(stims)
    finally:
        ds.close()
    # check 0: Remy's own per-trial window means (respvec/xrds_mean_peak.nc)
    c0 = {}
    if rp.exists():
        mp = xr.open_dataset(rp)
        try:
            if not np.array_equal(mp['cells'].to_numpy().astype(int), rois):
                raise RuntimeError(f'{acq}: respvec ROIs differ')
            for v in VARIABLES:
                c0[v] = float(np.nanmax(np.abs(remy_resp[v] - mp[v].to_numpy().astype(float))))
        finally:
            mp.close()
    idx = pd.MultiIndex.from_product([[stage], [acq], range(n_tr), rois], names=['stage', 'acq', 'trial', 'roi'])
    wm = pd.DataFrame(cols, index=idx).reset_index()
    wm.insert(4, 'odor', np.repeat(odors_t, len(rois)))
    wm.insert(5, 'repeat', np.repeat(rep, len(rois)))
    ns = pd.DataFrame(noise)
    ns.insert(0, 'roi', rois)
    ns.insert(0, 'acq', acq)
    ns.insert(0, 'stage', stage)
    lock = dict(stage=stage, acq=acq, file=file_entry(path, sha), respvec=file_entry(rp, sha) if rp.exists() else None,
                n_rois=int(len(rois)), n_trials=n_tr, frame_times=[float(t[0]), float(t[-1]), float(np.median(np.diff(t)))],
                n_frames=int(len(t)), remy_respvec_baseline_subtracted=remy_base, odors=odors, trial_odors=odors_t, attrs=attrs)
    return dict(window_means=wm, noise_sd=ns, frames=frames, t=t, rdms=rdms, odors=odors, lock=lock, check0=c0)


def remy_aggregate(stage: str, root: Path = REMY_ROOT):
    """Remy's per-fly stim RDMs (xarray Dataset, acq x 17 x 17) and the acq folder name of each of its flies."""
    import xarray as xr
    ds = xr.open_dataset(root / 'data_by_acq' / STAGE_DIRS[stage] / 'aggregate' / 'xrds_stim_rdm_concat.nc')
    names = [str(t).replace('/', '__').replace('__', '__fly0', 1) for t in ds['title_str'].to_numpy()]
    return ds, names


def _acq_key(acq: str) -> Tuple[str, int, str]:
    """('2023-04-11', 3, 'megamat1') from '2023-04-11__fly03__megamat1' (any '__dsub3' suffix kept in the name)."""
    p = acq.split('__')
    return p[0], int(p[1].replace('fly', '')), '__'.join(p[2:])


def check_remy_rdms(stage: str, rdms: Dict[str, Dict[str, np.ndarray]], odors: Sequence[str],
                    root: Path = REMY_ROOT) -> pd.DataFrame:
    """Check 1: our stim RDMs (Remy's window, every variable; {acq: {variable: odors x odors}}) against her aggregate
    xrds_stim_rdm_concat.nc, fly by fly (matched by date, fly and recording name). Rows: stage, acq, variable,
    max_abs_diff, passed."""
    ds, _ = remy_aggregate(stage, root)
    rows = []
    try:
        lab = [odor_of_stim(s) for s in ds['stim_row'].to_numpy()]
        ix = [lab.index(o) for o in odors]
        keys = list(zip([str(x) for x in ds['date_imaged'].to_numpy()], [int(x) for x in ds['fly_num'].to_numpy()],
                        [str(x) for x in ds['thorimage_name'].to_numpy()]))
        for acq, per_var in rdms.items():
            k = _acq_key(acq)
            if k not in keys:
                rows.append(dict(stage=stage, acq=acq, variable='', max_abs_diff=np.nan, passed=False,
                                 note='not in the aggregate'))
                continue
            i = keys.index(k)
            for v, R in per_var.items():
                theirs = ds[v].to_numpy()[i][np.ix_(ix, ix)]
                d = float(np.nanmax(np.abs(R - theirs)))
                rows.append(dict(stage=stage, acq=acq, variable=v, max_abs_diff=d, passed=bool(d <= EXACT_TOL), note=''))
        for k in keys:
            if f'{k[0]}__fly{k[1]:02d}__{k[2]}' not in rdms:
                rows.append(dict(stage=stage, acq=f'{k[0]}__fly{k[1]:02d}__{k[2]}', variable='', max_abs_diff=np.nan,
                                 passed=False, note='in the aggregate, not read'))
    finally:
        ds.close()
    return pd.DataFrame(rows)


def dashboard_window_means(stage: str, window: Tuple[float, float] = DASH_WINDOW) -> Tuple[np.ndarray, List[str]]:
    """The decoded dashboard's per-fly `stim` matrices (DIST_VAR) averaged over the window's frames (flies x 17 x 17,
    1 - r), and their odors (E90's `remy_window_means`)."""
    from experiments import E90_model_stages as E90
    a, t, odors = E90.remy_matrices(stage, 'stim')
    return a[:, frame_mask(t, window)].mean(axis=1), odors


def match_flies(ours: Dict[str, np.ndarray], theirs: np.ndarray) -> pd.DataFrame:
    """Pearson r over the upper-triangle pairs between each of our flies and each dashboard fly, and the one-to-one
    match maximizing the summed r (Hungarian on -r). Rows: acq, dash_fly, r, r_second (our fly's next best), matched."""
    from scipy.optimize import linear_sum_assignment
    acqs = list(ours)
    n = theirs.shape[1]
    iu = np.triu_indices(n, 1)
    R = np.array([[np.corrcoef(ours[a][iu], theirs[j][iu])[0, 1] for j in range(theirs.shape[0])] for a in acqs])
    ri, ci = linear_sum_assignment(-np.nan_to_num(R, nan=-2.0))
    rows = []
    for i, a in enumerate(acqs):
        order = np.sort(R[i])[::-1]
        m = ci[list(ri).index(i)] if i in ri else -1
        rows.append(dict(acq=a, dash_fly=int(m) + 1 if m >= 0 else -1, r=float(R[i, m]) if m >= 0 else np.nan,
                         r_second=float(order[1]) if len(order) > 1 else np.nan, matched=bool(m >= 0)))
    return pd.DataFrame(rows)


def check_dashboard(stage: str, frames: Dict[str, np.ndarray], t: np.ndarray, odors: Sequence[str],
                    excluded: Sequence[str] = ()) -> pd.DataFrame:
    """Check 2: our per-frame RDMs' means over DASH_WINDOW against the decoded dashboard's, flies matched one to one
    (`match_flies`). Passes if every dashboard fly is matched at r >= DASH_R_MIN by one of our flies that is not
    `excluded`, and every excluded fly is unmatched with no dashboard fly at r >= DASH_R_MIN (the dashboard has 5 of
    the 6 bouton flies). Rows: stage, acq, dash_fly, r, r_second, r_best, matched, excluded, passed."""
    theirs, dodors = dashboard_window_means(stage)
    ix = [list(dodors).index(o) for o in odors]
    theirs = theirs[:, ix][:, :, ix]
    m = frame_mask(t, DASH_WINDOW)
    ours = {a: np.nanmean(f[m], axis=0) for a, f in frames.items()}
    mf = match_flies(ours, theirs)
    iu = np.triu_indices(len(odors), 1)
    mf['r_best'] = [max(np.corrcoef(ours[a][iu], th[iu])[0, 1] for th in theirs) for a in mf['acq']]
    mf.insert(0, 'stage', stage)
    mf['excluded'] = mf['acq'].isin(list(excluded))
    mf['passed'] = np.where(mf['excluded'], ~mf['matched'] & (mf['r_best'] < DASH_R_MIN),
                            mf['matched'] & (mf['r'] >= DASH_R_MIN))
    if int((mf['matched'] & ~mf['excluded']).sum()) < theirs.shape[0]:
        mf['passed'] = False
    return mf


README = """# Remy's per-ROI megamat recordings, reduced (E119 part 0)

Made by `data_prep/remy_per_roi_data.py` (model_experiment_queue; provenance.json gives the commit) from the copies in
`source/` (source/README.md: sources, sizes, sha256; copied at viki's user's instruction of 2026-11-20, Remy's data shared
with her OK) of Remy's analysis outputs on tensor, `{root}/data_by_acq/`:
- PN boutons: `pn_boutons/volumetric_without_neuropil_correction/<acq>/trials/xrds_bc_trials.nc`, {n_bouton} flies
  analyzed;
- KC claws: `kc_claws/volumetric_with_neuropil_correction/<acq>/trials/xrds_bc_trials.nc`, {n_claw} flies.

Excluded from every analysis: the bouton fly 2022-10-21 fly02 (`megamat1_001__dsub3`). Remy's dashboards and bouton
analysis leave it out: it matches none of the dashboard's five bouton flies (check 2), and its median cellprob is 0.03.
It is copied in `source/` and listed in `lock.json` (`excluded`), and checks 0 and 1 are run on it. It has no rows in
`window_means.parquet`, `noise_sd.parquet` or `frame_rdms.npz`.

Each file is trials (17 megamat odors x 3 repeats, all at 10^-3) x ROIs x frames (-5 to 19.8 s after odor onset, 0.2 s),
each trial's baseline (the median of -5 to 0 s) subtracted, with variables F, Fneu, Fc (F - 0.7 Fneu for the claws), spks,
F_zscore and Fc_zscore (each ROI z-scored over the session). Remy's figures use F_zscore for the boutons and Fc_zscore for
the claws. The 200 MB files stay where they are; this folder holds what E119 reads.

Files:
- `window_means.parquet`: one row per stage, acq, trial and ROI (`roi` is Remy's suite2p cell number), with the trial's
  `odor` (abbreviation) and `repeat` (0 to 2: the odor's occurrence number in the trial order), and for every variable and
  window a column `<variable>__<window>`: the mean over the window's frames. Windows are closed, in s after onset:
  {windows}. `remy` is Remy's respvec window, `e90` is E90's stage window, and `pre` is the pre-odor window without the
  onset frame.
- `noise_sd.parquet`: per stage, acq and ROI, for every variable, the SD of the frames in the 51 pre-odor windows
  ({pre}), each window's own mean removed. A response in baseline-noise units is (window mean - pre mean) / this SD,
  as the realism bands' data side reads Yang's KCs.
- `frame_rdms.npz`: per stage and acq (`<stage>__<acq>`), frames x 17 x 17 correlation distances (1 - Pearson r across
  ROIs of the trial-mean responses at each frame) of F_zscore (boutons) or Fc_zscore (claws); `t` the frame times;
  `odors` the order.
- `lock.json`: every file read (path, size, mtime, sha256), each fly's ROI count, frame times, odors and trial order,
  the files' attributes, and the checks.
- `checks.csv`: the checks below, one row per fly and variable.
- `provenance.json`: script, commits, when, settings.

Checks (all passed, or the script writes nothing but the checks):
0. Remy's own per-trial window means (`respvec/xrds_mean_peak.nc`, 0.5 to 2.5 s) are reproduced to {tol:g} for every
   variable. For the boutons her respvec subtracts the mean over -5 to 0 s (its `respvec.baseline_method` is 'mean'); for
   the claws it subtracts nothing. `lock.json` records which (`remy_respvec_baseline_subtracted`).
1. Remy's per-fly stim distance matrices (`aggregate/xrds_stim_rdm_concat.nc`: 1 - Pearson r across ROIs of the
   trial-mean 0.5 to 2.5 s window means) are reproduced to {tol:g} for every variable and fly.
2. The decoded dashboard matrices (`data/remy_rdm_timeseries/`, per-frame `stim` matrices, averaged over
   {dash[0]} to {dash[1]} s) are reproduced from our per-frame matrices: flies matched one to one, Pearson r over the 136
   pairs of at least {dash_r} for every dashboard fly. The dashboard's claw frames are 0.25 s apart, ours 0.2 s, so
   this check is not exact. The dashboard has {n_dash_bouton} of the {n_bouton} bouton flies.

Not here: the same fly's claws and somata (megamat1_redo, Remy's suite2p output on matrix; E119 part 0b), which wait
until that source is cleared.
"""


def write_readme(out: Path, n: Dict[str, int], n_dash_bouton: int) -> None:
    (out / 'README.md').write_text(README.format(root=REMY_ROOT, n_bouton=n.get('bouton', 0), n_claw=n.get('claw', 0),
        windows=', '.join(f'`{k}` {a:g} to {b:g}' for k, (a, b) in WINDOWS.items()), pre=f'{PRE_NOISE[0]:g} to {PRE_NOISE[1]:g} s',
        tol=EXACT_TOL, dash=DASH_WINDOW, dash_r=DASH_R_MIN, n_dash_bouton=n_dash_bouton))


# ---------------------------------------------------------------------------
# the copy of the sources (viki's user, 2026-11-20: Remy's data, shared with her OK)

def source_copy_dir(out: Optional[Path] = None) -> Path:
    return (out or out_dir()) / 'source'


def source_files(root: Path = REMY_ROOT, same_fly: Path = SAME_FLY_SRC) -> List[Tuple[Path, str]]:
    """(source path, path under source/) of every file the copy holds: per stage every acquisition's
    trials/xrds_bc_trials.nc and respvec/ files and the stage's aggregate/ files (relative to <root>/data_by_acq/); the
    same fly's suite2p combined/ files and acquisition files (under same_fly/<SAME_FLY_NAME>/). The copy of
    <root>/data_by_acq/ is source/data_by_acq/, so `acquisitions(source_copy_dir())` reads the copy."""
    out = []
    base = root / 'data_by_acq'
    for sub in STAGE_DIRS.values():
        for p in sorted((base / sub).glob('*/trials/xrds_bc_trials.nc')):
            out.append((p, 'data_by_acq/' + str(p.relative_to(base))))
            for q in sorted((p.parent.parent / 'respvec').glob('*')):
                out.append((q, 'data_by_acq/' + str(q.relative_to(base))))
        for q in sorted((base / sub / 'aggregate').glob('*')):
            if q.is_file():
                out.append((q, 'data_by_acq/' + str(q.relative_to(base))))
    comb = same_fly / 'source_extraction_s2p' / 'suite2p' / 'combined'
    for f in SAME_FLY_COMBINED:
        if (comb / f).exists():
            out.append((comb / f, f'same_fly/{SAME_FLY_NAME}/suite2p_combined/{f}'))
    for f in SAME_FLY_ACQ_FILES:
        if (same_fly / f).exists():
            out.append((same_fly / f, f'same_fly/{SAME_FLY_NAME}/{f}'))
    return out


SOURCES_README = """# Copies of Remy's per-ROI megamat recordings (E119)

Copied {when} by `data_prep/remy_per_roi_data.py --copy` at viki's user's instruction of 2026-11-20 ("you can read
remy's bouton and claw data for all those flies and put in the shared folder. she's fine with it and it's all the lab's
data anyway"): Remy's data, shared with her OK. The sources were only read; every file is a byte copy (sha256 below).

Sources:
- `data_by_acq/...`: {root}/data_by_acq/ (Remy's analysis outputs on tensor): every acquisition's
  `trials/xrds_bc_trials.nc` and `respvec/`, and each stage's `aggregate/`.
- `same_fly/{same}/`: {same_src}/ (Remy's processed data on matrix), 2023-04-11 fly 3, megamat1_redo, the recording of
  the same fly's claws and somata: suite2p's `combined/` files (as `suite2p_combined/`) and the acquisition files.
  Not the raw movie.

Excluded from every analysis (copied all the same): {excluded}.

| file | source | bytes | sha256 |
|---|---|---:|---|
{rows}
"""


def copy_sources(out: Optional[Path] = None, root: Path = REMY_ROOT, same_fly: Path = SAME_FLY_SRC) -> pd.DataFrame:
    """Byte-copies `source_files` into source/ (shutil.copyfile: no mode or owner carried over; the sources are only
    read), skipping a file already copied with the same size and sha256, and writes source/README.md and
    source/files.csv (file, source, bytes, sha256)."""
    import shutil
    dest = source_copy_dir(out)
    rows = []
    for src, rel in source_files(root, same_fly):
        d = dest / rel
        d.parent.mkdir(parents=True, exist_ok=True)
        sha = _sha256(src)
        if not (d.exists() and d.stat().st_size == src.stat().st_size and _sha256(d) == sha):
            shutil.copyfile(src, d)
            if _sha256(d) != sha:
                raise RuntimeError(f'{rel}: the copy differs from its source')
        rows.append(dict(file=rel, source=str(src), bytes=int(src.stat().st_size), sha256=sha))
        print(f'copied {rel}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(dest / 'files.csv', index=False)
    when = datetime.now().astimezone().isoformat(timespec='seconds')
    (dest / 'README.md').write_text(SOURCES_README.format(
        when=when, root=root, same=SAME_FLY_NAME, same_src=same_fly,
        excluded='; '.join(f'{a} ({s}): {why}' for (s, a), why in EXCLUDED.items()),
        rows='\n'.join(f'| `{r.file}` | `{r.source}` | {r.bytes} | `{r.sha256}` |' for r in df.itertuples())))
    return df


# ---------------------------------------------------------------------------
# part 0b: the same fly's claws and somata (megamat1_redo)

def suite2p_traces(folder: Path) -> dict:
    """F, Fneu (ROIs x frames) and iscell (ROIs x 2) of a suite2p 'combined' folder (F.npy, Fneu.npy, iscell.npy),
    restricted to iscell's first column."""
    folder = Path(folder)
    F = np.load(folder / 'F.npy')
    Fn = np.load(folder / 'Fneu.npy')
    ic = np.load(folder / 'iscell.npy')
    keep = np.asarray(ic)[:, 0].astype(bool)
    return dict(F=F[keep], Fneu=Fn[keep], cells=np.flatnonzero(keep))


def cut_trials(traces: np.ndarray, frame_times: np.ndarray, onsets: Sequence[float],
               t: Optional[np.ndarray] = None, baseline: Tuple[float, float] = (-5.0, 0.0)) -> np.ndarray:
    """Cut traces (ROIs x frames, at `frame_times`) into trials (trials x ROIs x len(t)) around the odor `onsets`, by
    linear interpolation onto t (default Remy's -5 to 19.8 s at 0.2 s), each trial's baseline (the median over the
    closed `baseline` window of t) subtracted, as Remy's xrds_bc_trials."""
    if t is None:
        t = np.round(np.arange(-5.0, 19.81, 0.2), 6)
    X = np.asarray(traces, dtype=float)
    ft = np.asarray(frame_times, dtype=float)
    b = frame_mask(t, baseline)
    out = np.empty((len(onsets), X.shape[0], len(t)))
    for k, o in enumerate(onsets):
        tt = o + np.asarray(t)
        if tt[0] < ft[0] or tt[-1] > ft[-1]:
            raise ValueError(f'trial {k}: its window leaves the recording')
        y = np.stack([np.interp(tt, ft, x) for x in X])
        out[k] = y - np.median(y[:, b], axis=1, keepdims=True)
    return out


def session_zscore(x: np.ndarray) -> np.ndarray:
    """Each ROI's trials (trials x ROIs x frames) z-scored over all of its frames (as Remy's *_zscore of a trial-cut
    session)."""
    x = np.asarray(x, dtype=float)
    m = x.mean(axis=(0, 2), keepdims=True)
    s = x.std(axis=(0, 2), keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (x - m) / s


def trial_pair_rdms(x: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """Per frame, 1 - the mean over the trial pairs of two odors of the single-trial Pearson r across ROIs (the same-fly
    series' estimator, 'mean over the 9 trial pairs'); on the diagonal the mean over the pairs of different trials
    (frames x odors x odors)."""
    x = np.asarray(x, dtype=float)
    z = x - x.mean(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = z / np.linalg.norm(z, axis=1, keepdims=True)
    C = np.einsum('iut,jut->tij', z, z)                                 # frames x trials x trials
    o = np.asarray(odors_of_trials)
    out = np.empty((x.shape[2], len(odors), len(odors)))
    for a, oa in enumerate(odors):
        ia = np.flatnonzero(o == oa)
        for b, ob in enumerate(odors):
            ib = np.flatnonzero(o == ob)
            blk = C[:, ia][:, :, ib]
            if a == b:
                off = ~np.eye(len(ia), dtype=bool)
                out[:, a, b] = 1.0 - blk[:, off].mean(axis=1)
            else:
                out[:, a, b] = 1.0 - blk.mean(axis=(1, 2))
    return out


def same_fly_checks(claw_frames: np.ndarray, soma_frames: np.ndarray, t: np.ndarray, odors: Sequence[str],
                    decoded_dir: Optional[Path] = None) -> dict:
    """Part 0b's gate: our same-fly per-frame matrices (`trial_pair_rdms`, frames x odors x odors at times t) against the
    decoded same-fly series (data/remy_rdm_timeseries/same_fly/), each at r >= SAME_FLY_R_MIN over the pairs of every
    frame both have (mean over frames), and the claws-against-somata agreement over 1 to 2 s (Pearson r of the window
    means' upper triangles) within SAME_FLY_AGREEMENT_TOL of Remy's 0.86."""
    d = Path(decoded_dir) if decoded_dir is not None else H.shared_dir() / 'data' / 'remy_rdm_timeseries' / 'same_fly'
    dodors = json.loads((d.parent / 'odors.json').read_text())
    ix = [list(dodors).index(o) for o in odors]
    td = np.load(d / 'times_samefly.npy')
    n = len(odors)
    iu = np.triu_indices(n, 1)
    res = {}
    for name, ours in (('claws', claw_frames), ('somata', soma_frames)):
        th = np.load(d / f'kc_{name}_samefly_corrdist.npy')[:, ix][:, :, ix]
        rs = []
        for k, tk in enumerate(td):
            j = np.flatnonzero(np.abs(np.asarray(t) - tk) < 1e-6)
            if len(j) and np.isfinite(ours[j[0]][iu]).all():
                rs.append(np.corrcoef(ours[j[0]][iu], th[k][iu])[0, 1])
        res[f'{name}_r'] = float(np.mean(rs)) if rs else np.nan
        res[f'{name}_n_frames'] = len(rs)
    m = frame_mask(t, (1.0, 2.0))
    agree = float(np.corrcoef(np.nanmean(claw_frames[m], axis=0)[iu], np.nanmean(soma_frames[m], axis=0)[iu])[0, 1])
    res['agreement_1_2s'] = agree
    res['passed'] = bool(res['claws_r'] >= SAME_FLY_R_MIN and res['somata_r'] >= SAME_FLY_R_MIN
                         and abs(agree - SAME_FLY_AGREEMENT) <= SAME_FLY_AGREEMENT_TOL)
    return res


SAME_FLY_DECODED_T = np.round(np.arange(-5.0, 19.76, 0.25), 6)     # the decoded same-fly series' frame grid
NEUCOEFF = 0.7


def same_fly_sessions(src: Path) -> dict:
    """The same fly's continuous traces from the copy (src = source/same_fly/<SAME_FLY_NAME>): per kind ('claw' from
    iscell_calyx, 'soma' from iscell_soma) Fc = F - NEUCOEFF Fneu (ROIs x volumes) and its z-score over the whole
    session per ROI; the volume times (timestamps' stack_times), the odor onsets (olf_ict) and the trial odors
    (stim_list_flatstr)."""
    comb = src / 'suite2p_combined'
    F = np.load(comb / 'F.npy')
    Fn = np.load(comb / 'Fneu.npy')
    ts = np.load(src / 'timestamps.npy', allow_pickle=True).item()
    stims = json.loads((src / 'stim_list.json').read_text())['stim_list_flatstr']
    vt = np.asarray(ts['stack_times'], dtype=float)
    onsets = np.asarray(ts['olf_ict'], dtype=float)
    if F.shape[1] != len(vt) or len(onsets) != len(stims):
        raise RuntimeError(f'{F.shape[1]} volumes for {len(vt)} stack times; {len(onsets)} onsets for {len(stims)} trials')
    out = dict(t=vt, onsets=onsets, odors_of_trials=[odor_of_stim(s) for s in stims])
    for kind, f in (('claw', 'iscell_calyx.npy'), ('soma', 'iscell_soma.npy')):
        keep = np.load(comb / f)[:, 0].astype(bool)
        Fc = F[keep] - NEUCOEFF * Fn[keep]
        out[f'{kind}_Fc'] = Fc
        out[f'{kind}_Fc_zscore'] = (Fc - Fc.mean(axis=1, keepdims=True)) / Fc.std(axis=1, keepdims=True)
        out[f'{kind}_cells'] = np.flatnonzero(keep)
    return out


def part_0b(out: Optional[Path] = None) -> dict:
    """Part 0b: the same fly's claw and soma trials cut from the copied suite2p output (`same_fly_sessions`,
    `cut_trials`), the gate (`same_fly_checks` on the decoded series' 0.25 s grid), and, if it passes,
    same_fly_trials.npz: per kind the trials (repeats x ROIs x 17 odors) in baseline-noise units ('<kind>_noise': the
    0.5-3 s window's change over the pre-odor frame SD of Fc, as part 0) and as the change of the session z ('<kind>_z');
    same_fly_checks.json either way."""
    out = out or out_dir()
    src = source_copy_dir(out) / 'same_fly' / SAME_FLY_NAME
    S = same_fly_sessions(src)
    o_t = S['odors_of_trials']
    odors = sorted(set(o_t))
    rep = trial_repeats(o_t)
    fr = {k: trial_pair_rdms(cut_trials(S[f'{k}_Fc_zscore'], S['t'], S['onsets'], SAME_FLY_DECODED_T), o_t, odors)
          for k in ('claw', 'soma')}
    chk = same_fly_checks(fr['claw'], fr['soma'], SAME_FLY_DECODED_T, odors)
    chk.update(n_claw=int(len(S['claw_cells'])), n_soma=int(len(S['soma_cells'])), neucoeff=NEUCOEFF,
               generated=datetime.now().astimezone().isoformat(timespec='seconds'))
    (out / 'same_fly_checks.json').write_text(json.dumps(chk, indent=1, default=str) + '\n')
    print(json.dumps(chk, indent=1, default=str), flush=True)
    if not chk['passed']:
        return chk
    t = np.round(np.arange(-5.0, 19.81, 0.2), 6)
    res = dict(odors=np.asarray(odors), t_window=np.asarray(WINDOWS['e90']), pre=np.asarray(PRE_NOISE))
    for k in ('claw', 'soma'):
        tr_fc = cut_trials(S[f'{k}_Fc'], S['t'], S['onsets'], t)
        tr_z = cut_trials(S[f'{k}_Fc_zscore'], S['t'], S['onsets'], t)
        sd = noise_sd(tr_fc, t)
        keep = sd > 0
        for name, X, scale in (('noise', tr_fc, sd), ('z', tr_z, None)):
            w = window_means(X, t, {'w': WINDOWS['e90'], 'p': WINDOWS['pre']})
            ch = w['w'] - w['p']                                        # trials x ROIs
            if scale is not None:
                ch = ch / np.where(scale > 0, scale, np.nan)[None]
            ch = ch[:, keep]
            arr = np.stack([np.stack([ch[(np.asarray(o_t) == o) & (rep == r)][0] for o in odors], axis=-1)
                            for r in range(rep.max() + 1)])            # repeats x ROIs x odors
            res[f'{k}_{name}'] = arr
        res[f'{k}_cells'] = S[f'{k}_cells'][keep]
    np.savez_compressed(out / 'same_fly_trials.npz', **res)
    return chk


# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None)
    ap.add_argument('--no-sha', action='store_true', help='skip the sha256 of the source files (tests)')
    ap.add_argument('--copy', action='store_true', help='copy the sources into source/ and stop')
    ap.add_argument('--part0b', action='store_true', help="the same fly's claw and soma trials (part 0b) and stop")
    a = ap.parse_args(argv)
    out = a.out or out_dir()
    out.mkdir(parents=True, exist_ok=True)
    if a.copy:
        copy_sources(out)
        return
    if a.part0b:
        chk = part_0b(out)
        sys.exit(0 if chk['passed'] else 1)
    when = datetime.now().astimezone().isoformat(timespec='seconds')
    root = source_copy_dir(out) if (source_copy_dir(out) / 'data_by_acq').is_dir() else REMY_ROOT
    acqs = acquisitions(root)
    acqs['excluded'] = [(r.stage, r.acq) in EXCLUDED for r in acqs.itertuples()]
    wms, nss, locks, checks = [], [], [], []
    frames_npz: Dict[str, np.ndarray] = {}
    per_stage: Dict[str, dict] = {s: dict(rdms={}, frames={}, t=None, odors=None) for s in STAGE_DIRS}
    for r in acqs.itertuples():
        print(f'{r.stage} {r.acq}', flush=True)
        red = reduce_acq(r.stage, r.acq, r.path, sha=not a.no_sha)
        red['lock']['excluded'] = EXCLUDED.get((r.stage, r.acq), False)
        locks.append(red['lock'])
        if not r.excluded:
            wms.append(red['window_means'])
            nss.append(red['noise_sd'])
            frames_npz[f'{r.stage}__{r.acq}'] = red['frames'].astype(np.float32)
        ps = per_stage[r.stage]
        ps['rdms'][r.acq] = red['rdms']
        ps['frames'][r.acq] = red['frames']
        if ps['t'] is None:
            ps['t'], ps['odors'] = red['t'], red['odors']
        elif not (np.array_equal(ps['t'], red['t']) and ps['odors'] == red['odors']):
            raise RuntimeError(f'{r.acq}: frame times or odors differ from the stage\'s first fly')
        for v, d in red['check0'].items():
            checks.append(dict(check='0_window_means', stage=r.stage, acq=r.acq, variable=v, value=d,
                               passed=bool(d <= EXACT_TOL), note='max abs diff from respvec/xrds_mean_peak.nc'))
    n_dash = {}
    for stage, ps in per_stage.items():
        c1 = check_remy_rdms(stage, ps['rdms'], ps['odors'], root)
        for x in c1.itertuples():
            checks.append(dict(check='1_remy_stim_rdm', stage=stage, acq=x.acq, variable=x.variable, value=x.max_abs_diff,
                               passed=bool(x.passed), note=x.note or 'max abs diff from aggregate/xrds_stim_rdm_concat.nc'))
        c2 = check_dashboard(stage, ps['frames'], ps['t'], ps['odors'],
                             excluded=[a for (s_, a) in EXCLUDED if s_ == stage])
        n_dash[stage] = int(c2['matched'].sum())
        for x in c2.itertuples():
            checks.append(dict(check='2_dashboard', stage=stage, acq=x.acq, variable=DIST_VAR[stage],
                               value=x.r if x.matched else x.r_best, passed=bool(x.passed),
                               note=(f'dashboard fly {x.dash_fly}; next best r {x.r_second:.3f}' if x.matched else
                                     f'excluded; matches no dashboard fly (best r {x.r_best:.3f})')))
    ck = pd.DataFrame(checks)
    ck.to_csv(out / 'checks.csv', index=False)
    ok = bool(ck['passed'].all())
    print(ck.groupby(['check', 'stage'])['passed'].agg(['sum', 'size']).to_string(), flush=True)
    if not ok:
        print('checks failed:\n' + ck[~ck['passed']].to_string(), flush=True)
        sys.exit(1)
    pd.concat(wms, ignore_index=True).to_parquet(out / 'window_means.parquet')
    pd.concat(nss, ignore_index=True).to_parquet(out / 'noise_sd.parquet')
    np.savez_compressed(out / 'frame_rdms.npz', t=per_stage['claw']['t'], odors=np.asarray(per_stage['claw']['odors']),
                        **frames_npz)
    lock = dict(generated=when, root=str(root), original_root=str(REMY_ROOT),
                excluded={f'{k[0]}__{k[1]}': v for k, v in EXCLUDED.items()}, variables=list(VARIABLES), windows=WINDOWS, pre_noise=PRE_NOISE,
                dist_var=DIST_VAR, noise_var=NOISE_VAR, acquisitions=locks,
                checks=dict(passed=ok, n=int(len(ck)), by_check=ck.groupby('check')['passed'].all().to_dict()))
    (out / 'lock.json').write_text(json.dumps(lock, indent=1, default=str) + '\n')
    n_used = {s: int((~acqs['excluded'] & (acqs['stage'] == s)).sum()) for s in per_stage}
    write_readme(out, n_used, n_dash.get('bouton', 0))
    prov = dict(script=H.script_path(__file__), generated=when, **H.code_commits(),
                settings=dict(windows=WINDOWS, pre_noise=PRE_NOISE, exact_tol=EXACT_TOL, dash_window=DASH_WINDOW,
                              dash_r_min=DASH_R_MIN))
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print(f'wrote {out}', flush=True)


if __name__ == '__main__':
    main()
