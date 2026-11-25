#!/usr/bin/env python3
"""E137, the held-out panels read like for like: the mixture recordings by position in the session, and validation2 on
all its flies and odors (no simulation).

The entry is item E137 of experiment_priorities.md (proposed by viki 2026-11-24 22:54 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-25 02:08 under the idle-host rule). "The family", the views,
one-fly units and the floor are E134's (kcq/population_fit.py, kcq/stage_readings.py, experiments/
E134_population_baseline.py). No model is run and no model parameter changes: E134's saved runs are re-read and
re-scored.

Parts (CLI `--part`, one per design item):
- lock (design items 1 and 3): for the eleven mixture flies of olfdb dataset 8, Remy's per-trial processed files (her
  suite2p 'combined' respvec files on matrix, xrds_suite2p_respvec_mean_peak.nc and xrds_suite2p_respvec_baseline_std.nc;
  the copy in her kc_soma_nls folders on tensor, respvec/xrds_mean_peak.nc, where there is one) with their sha256, the
  session's stimulus file (the file's own trial list; Remy's olfactometer config of the session and the generator config
  of its movie type, with randomize_presentation_order), and the reproduction of the uploaded table from them (the
  cells, trials and blocks it averages). Each fly's per-trial values of the cells the uploaded table holds are reduced
  to one parquet file in <shared>/data/remy_mixture_kc_trials/ (the lock records its sha256), so that a host without
  matrix reads the same values. A fly whose files are not found is listed and stays read as uploaded. The five readings
  of the flies (`READINGS`), their units and criteria and the footing are fixed here, before any score. Written to
  results/E137_heldout_readings/lock.json and lock_flies.csv; the committed copy is experiments/E137_lock.json.
- position (design item 2): per fly, stimulus, block and presentation, the response rate and mean response; the third
  over the first within a block, the second block over the first, the probes and the solvent at the session's start
  and end; whether the decline is general or specific to the mixture's components.
- readings (design items 3 and 4): each reading's flies, per panel and units: the floor, one recording against the
  others and the response rates.
- rescore (design items 4 and 5): E134's saved mixture runs (the family as fitted and held on its eight wirings, the five
  beside models) read as E134 reads the KC somata (`kcq.stage_readings`' default reading: the curve set on the reading's
  flies, the flies' criterion), with the flies' trial noise added (`TrialNoise`), scored under each reading, beside the
  flat profile and the ORN recordings as the prediction; every strength error split into level, spread and order
  (E136's `strength_split`). A check reads the lock's as-uploaded flies with no noise and compares with E134's
  score_scores.csv.
- rise (design item 6): the mixture's rise over its three dilutions (response rate and mean response at 0 and -1 over
  -2) and the mixture over its components, the flies' under each reading and the models'.
- validation2 (design item 7): validation2 on all its 20 flies and 22 odors under the two partial-panel rules of hal's
  `kcq.population_fit.partial_reference` ('anchor' on paa, PEA and guai; 'fly_scale'), the six-fly reading of E134
  and the two noise groups beside.
- order (design item 8): for megamat and validation2, the share of the strength unit an odor's place in its block
  carries; for the APL16 recordings, whether Yang's processed tables hold the order (otherwise a question for Yang).
- classes: the classes the entry fixes beforehand (`classify`), the outcome they point to, summary_strength.csv.
- figures: position_rates.png and strength_by_reading.png.

The public interface other experiments read (E138's mixture rows; signatures stable from 13bfca9 on):
- `load_lock()`: the committed lock, per fly its files with their sha256, its stimulus order, its reproduction and
  the readings it can be read under;
- `fly_mixture_view(panel, reading, trials=True, lock=None, units=None)`: the flies' populations of one mixture panel's
  KC somata under a named reading (one `kcq.population_fit.Population` per fly, with the single trials where on hand).

The readings (`READINGS`; design item 3). Remy's mixture sessions are not randomized: solvent, the five components (each
three consecutive presentations, a "block"), the mixture at -2, -1 and 0, then (July sessions) the partial mixtures and
the components again, and in the "probes" sessions two probe odors at the start and the end. The mixtures have one
block each. A reading takes, of the panel's eight stimuli:
- 'as_uploaded': olfdb dataset 8 as uploaded (each stimulus's first block, session z), with the single trials of that
  block for the noise;
- 'first_presentations': every stimulus's first presentation of the session (one trial);
- 'components_first_block': every stimulus's first block (the uploaded table's trials, in the reading's units);
- 'components_second_block': the components' second block where there is one, else their first (the two March / April
  control flies have one block); the mixtures' block;
- 'components_both_blocks': the components' two blocks together where there are two, else the first; the mixtures'
  block.
Each reading has its own flies' reference, so its own one-fly units and floor.

Units (design item 4). 'noise': baseline-noise units on megamat's footing (data_prep/stage_population_data.py's rule for
validation2): a cell's response (Remy's mean_peak Fc: the 2 to 8 s mean minus the -5 to 0 s mean) over the root mean
square of its per-trial baseline SD (her baseline_std file) over the session's trials, times megamat's median noise
bridge over the eleven mixture flies' median (`FOOTING_GROUPS`), at megamat's soma level (the E134 lock's 5.9486).
'z': Remy's session z (Fc_zscore), at the E134 lock's mixture criterion (1.9086). The position readings are read in
'noise' (`PRIMARY_UNITS`), 'as_uploaded' in 'z'; every reading is also scored in the other units beside.

Trial noise (`TrialNoise`, design item 4): the flies' single-trial residual SD about each cell's block mean, pooled over
the panel's flies, cells, stimuli and the reading's blocks (the n / (n - 1) correction), in the reading's units; a model
value is its value after the curve plus Gaussian noise of that SD over sqrt(n), n the number of trials the reading
averages for that stimulus (the flies' median over flies); seed `stage_readings.NOISE_SEED` + the wiring's index.

Holds (CLAUDE.md): Remy's processed files only (her per-trial respvec files and stimulus configs; no movie is read),
Yang's processed products only; no connectome table; no olfsysm change; nothing fitted to a held-out panel (the one
fitted quantity is a scale per fly among validation2's own flies, one of two candidate readings of the flies).

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part lock     # then commit the lock copy
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part position
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part readings
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part rescore
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part rise
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part validation2
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part order
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part classes
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part figures
"""
from __future__ import annotations

import argparse
import hashlib
import json
import socket
import subprocess
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E137_heldout_readings'
SCRIPT = 'experiments/E137_heldout_readings.py'
HERE = Path(__file__).resolve().parents[1]
LOCK_COPY = HERE / 'experiments' / 'E137_lock.json'

MIXTURE_PANELS: Tuple[str, ...] = ('kmix', 'cmix')
# the five readings of the flies (design item 3), fixed in the lock before any score
READINGS: Tuple[str, ...] = ('as_uploaded', 'first_presentations', 'components_first_block',
                             'components_second_block', 'components_both_blocks')
POSITION_READINGS: Tuple[str, ...] = READINGS[1:]
UNITS: Tuple[str, ...] = ('noise', 'z')
# each reading's units (the class reads these); the other units are scored beside
PRIMARY_UNITS: Dict[str, str] = {'as_uploaded': 'z', **{r: 'noise' for r in POSITION_READINGS}}
# validation2's two candidate rules for flies with partial panels (design item 7; hal's `partial_reference` rules)
V2_RULES: Tuple[str, ...] = ('common_mean', 'fly_scale')
V2_RULE_TO_P: Dict[str, str] = {'common_mean': 'anchor', 'fly_scale': 'fly_scale'}
V2_COMMON_ODORS: Tuple[str, ...] = ('paa', 'PEA', 'guai')

# Remy's files
REMY_MATRIX_MIX = Path('/local/matrix/Remy-Data/projects/natural_mixtures/processed_data')
MEAN_PEAK = 'xrds_suite2p_respvec_mean_peak.nc'
BASELINE_STD = 'xrds_suite2p_respvec_baseline_std.nc'
# her names of the panel's stimuli (as the respvec files give them) in the panel's (olfdb's, through
# stage_population_data._olfdb_mix_name)
REMY_NAMES: Dict[str, str] = {'IaOH': 'IAol', 'kiwi': 'kmix'}
SOLVENT = 'pfo'
PROBES: Tuple[str, ...] = ('t2h', '3mt1p')
CELL_RULE = "the cells whose xid0 is one of the mean_peak file's good_xid (Remy's good clusters)"
UPLOADED_RULE = "each stimulus's first block (its first three presentations of the session), the mean of Fc_zscore"
REPRO_TOL = 1e-5                 # the uploaded table has float32 values (about 7 significant digits)
# the footing of the 'noise' units: one group of all eleven mixture flies (one rig and pipeline)
FOOTING_GROUPS = 'all'
N_FLOOR = 1000
# the classes (E134's stage rule, on strength by odor: the adopted measure)
INSIDE_MIN_WIRINGS = 6
FAR_ERROR = 2.0
FAMILY_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
BESIDE_MODELS: Tuple[str, ...] = ('published_hemibrain', 'uniform_draw', 'bouton_baseline', 'reference',
                                  'firing_ceiling')
FAMILY_GROUPS: Tuple[str, ...] = ('fitted', 'held')


def shared() -> Path:
    from kcq import harness as H
    return H.shared_dir()


def out_dir() -> Path:
    """<shared>/results/E137_heldout_readings/ (created)."""
    from kcq import harness as H
    return H.results_dir(EXPERIMENT)


def data_dir() -> Path:
    """<shared>/data/remy_mixture_kc_trials/ (the reduced per-trial files of the eleven flies)."""
    return shared() / 'data' / 'remy_mixture_kc_trials'


def remy_tensor_kc() -> Path:
    """Remy's kc_soma_nls folders on tensor (natural_mixtures/analysis_outputs/data_by_acq/kc_soma_nls)."""
    return shared().parents[1] / 'Remy' / 'natural_mixtures' / 'analysis_outputs' / 'data_by_acq' / 'kc_soma_nls'


def olf_configs_dir() -> Path:
    """Remy's olfactometer configs on tensor (natural_mixtures/olfactometer_configs)."""
    return shared().parents[1] / 'Remy' / 'natural_mixtures' / 'olfactometer_configs'


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def sha256_file(path: Path, block: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(block), b''):
            h.update(b)
    return h.hexdigest()


def _rel(path: Path) -> str:
    from data_prep import stage_population_data as D
    return D.rel_path(Path(path))


def _abs(rel: str) -> Path:
    """A lock path back to a path on this host ('<shared>/...', '<queue>/...' or absolute)."""
    if rel.startswith('<shared>/'):
        return shared() / rel[len('<shared>/'):]
    if rel.startswith('<queue>/'):
        return HERE / rel[len('<queue>/'):]
    return Path(rel)


def file_entry(path: Path) -> dict:
    st = Path(path).stat()
    return dict(path=_rel(path), bytes=int(st.st_size),
                mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime)), sha256=sha256_file(path))


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    tmp.replace(path)


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git and has no uncommitted change."""
    try:
        r = subprocess.run(['git', '-C', str(HERE), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(HERE), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# stimuli

def _split_stim(s: str) -> Tuple[str, float]:
    name, conc = str(s).split(' @ ')
    return name.strip(), float(conc)


def panel_concentrations() -> Dict[str, Dict[str, float]]:
    """{panel: {stimulus label: log10 concentration}} of the uploaded table's stimuli (olfdb dataset 8's names)."""
    from data_prep import stage_population_data as D
    d = pd.read_csv(D.olfdb_path(D.OLFDB_DATASETS['kc_mixtures']), usecols=['panel', 'odor1']).drop_duplicates()
    out: Dict[str, Dict[str, float]] = {}
    for p, o in zip(d.panel, d.odor1):
        name, conc = _split_stim(o)
        out.setdefault(p, {})[D._olfdb_mix_name(o)] = conc
    return out


def stim_label(stim: str, panel: str, concs: Optional[Mapping[str, Mapping[str, float]]] = None) -> Optional[str]:
    """A respvec file's stimulus ('kiwi @ -2.0', 'IaOH @ -3.6') as the panel names it ('kmix-2', 'IAol'), None if it is
    not one of the panel's stimuli at the panel's concentration."""
    from data_prep import stage_population_data as D
    concs = panel_concentrations() if concs is None else concs
    name, conc = _split_stim(stim)
    name = REMY_NAMES.get(name, name)
    lab = D._olfdb_mix_name(f'{name} @ {conc:g}')
    c = concs.get(panel, {}).get(lab)
    return lab if c is not None and abs(c - conc) < 1e-6 else None


def session_blocks(stims: Sequence[str]) -> pd.DataFrame:
    """Per trial of a session (in order): the stimulus, its block (the index of its run of consecutive presentations
    among the stimulus's runs), its presentation within the block and within the session."""
    rows, runs, last, rep = [], {}, None, 0
    seen: Dict[str, int] = {}
    for i, s in enumerate(stims):
        if s != last:
            runs[s] = runs.get(s, -1) + 1
            rep = 0
        else:
            rep += 1
        rows.append(dict(trial=i, stim=str(s), block=runs[s], rep=rep, presentation=seen.get(s, 0)))
        seen[s] = seen.get(s, 0) + 1
        last = s
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# design item 1: the flies and their files

def mixture_flies() -> pd.DataFrame:
    """The eleven flies of olfdb dataset 8 (fly '<date>/<fly_num>', date, fly_num, panel, n_cells), from the E134
    lock's copy of the download."""
    from data_prep import stage_population_data as D
    d = pd.read_csv(D.olfdb_path(D.OLFDB_DATASETS['kc_mixtures']))
    g = d.groupby(['date', 'fly_num', 'panel']).roi.nunique().reset_index(name='n_cells')
    g['fly'] = [D._fly_label(a, b) for a, b in zip(g.date, g.fly_num)]
    return g[['fly', 'date', 'fly_num', 'panel', 'n_cells']].sort_values(['panel', 'fly']).reset_index(drop=True)


def uploaded_table(fly: str) -> pd.DataFrame:
    """One fly's uploaded table (olfdb dataset 8): cells (olfdb's roi, as int) x the panel's stimuli (session z)."""
    from data_prep import stage_population_data as D
    d = pd.read_csv(D.olfdb_path(D.OLFDB_DATASETS['kc_mixtures']))
    d = d[[D._fly_label(a, b) == fly for a, b in zip(d.date, d.fly_num)]]
    panel = str(d.panel.iloc[0])
    t = d.assign(odor=[D._olfdb_mix_name(o) for o in d.odor1]).pivot_table(index='roi', columns='odor', values='value')
    return t.reindex(columns=list(D.panel_odors(panel)))


def find_fly_files(fly: str, panel: str, root: Optional[Path] = None) -> dict:
    """Remy's per-trial processed files of one fly's mixture session: the matrix combined folder whose mean_peak file
    holds every stimulus of the panel ('mean_peak', 'baseline_std', 'folder', 'session'), the tensor kc_soma_nls copy
    of the same session where there is one ('tensor_mean_peak'), the session's olfactometer config ('config') and its
    movie type's generator config ('generator'); a missing file's entry is None."""
    import xarray as xr
    root = REMY_MATRIX_MIX if root is None else Path(root)
    date, num = fly.split('/')
    concs = panel_concentrations()
    want = set(concs[panel])
    out = dict(mean_peak=None, baseline_std=None, folder=None, session=None, tensor_mean_peak=None, config=None,
               generator=None, olf_config_name=None, movie_type=None)
    for c in sorted((root / date / num).glob(f'*/*/suite2p/combined/{MEAN_PEAK}')):
        with xr.open_dataset(c) as d:
            labs = {stim_label(s, panel, concs) for s in d['stim'].values}
            if not want <= labs:
                continue
            out.update(mean_peak=c, folder=c.parent, session=c.parents[3].name,
                       olf_config_name=str(d.attrs.get('olf_config', '')) or None,
                       movie_type=str(d.attrs.get('movie_type', '')) or None)
        b = c.parent / BASELINE_STD
        out['baseline_std'] = b if b.exists() else None
        break
    if out['session'] is not None:
        t = remy_tensor_kc() / f"{date}__fly{int(num):02d}__{out['session']}" / 'respvec' / 'xrds_mean_peak.nc'
        out['tensor_mean_peak'] = t if t.exists() else None
        if out['olf_config_name']:
            p = olf_configs_dir() / out['olf_config_name']
            out['config'] = p if p.exists() else None
        if out['movie_type']:
            p = olf_configs_dir() / f"{out['movie_type']}.yaml"
            out['generator'] = p if p.exists() else None
    return out


def config_order(path: Path) -> List[str]:
    """The stimulus of every trial of a session's olfactometer config (its pin_sequence through pins2odors), as
    '<name> @ <log10 conc>' with Remy's names (abbrev, else name)."""
    import yaml
    c = yaml.safe_load(Path(path).read_text())
    p2o = c.get('pins2odors') or {}
    out = []
    for g in c['pin_sequence']['pin_groups']:
        names = []
        for pin in g['pins']:
            o = p2o.get(pin) or p2o.get(str(pin))
            if not o:                       # a pin with no odor (the session's second, blank pin)
                continue
            names.append(f"{o.get('abbrev', o.get('name'))} @ {float(o.get('log10_conc', np.nan)):.1f}")
        out.append(', '.join(names))
    return out


def orders_match(config: Sequence[str], stims: Sequence[str]) -> bool:
    """Whether a config's trial order is the respvec file's: as many trials, the same concentration at every trial, and
    a one-to-one map between the config's names and the file's (Remy's files rename some: 'etoh' / 'EtOH', 'IAol' /
    'IaOH', 'kiwi approx.' / 'kiwi')."""
    if len(config) != len(stims):
        return False
    fwd: Dict[str, str] = {}
    back: Dict[str, str] = {}
    for a, b in zip(config, stims):
        try:
            ca, cb = [_split_stim(x)[1] for x in (a, b)]
        except ValueError:
            return False
        if abs(ca - cb) > 1e-6 or fwd.setdefault(a, b) != b or back.setdefault(b, a) != a:
            return False
    return True


def generator_settings(path: Optional[Path]) -> dict:
    """randomize_presentation_order and n_repeats of a movie type's generator config ({} when missing)."""
    if path is None:
        return {}
    import yaml
    c = yaml.safe_load(Path(path).read_text())
    return {k: c.get(k) for k in ('randomize_presentation_order', 'n_repeats') if k in c}


def reduce_fly(fly: str, panel: str, files: Mapping) -> pd.DataFrame:
    """One fly's per-trial values of the cells the uploaded table holds (`CELL_RULE`), all trials of the session: long
    rows cell, trial, stim (Remy's), label (the panel's name, else '<name> @ <conc>'), in_panel, block, rep,
    presentation, Fc (her mean_peak response), z (Fc_zscore), baseline_sd (her per-trial baseline SD of Fc)."""
    import xarray as xr
    concs = panel_concentrations()
    with xr.open_dataset(files['mean_peak']) as p, xr.open_dataset(files['baseline_std']) as b:
        good = set(int(x) for x in np.atleast_1d(p.attrs['good_xid']))
        sel = np.isin(p['xid0'].values, list(good))
        if not np.array_equal(p['cells'].values, b['cells'].values) or not np.array_equal(p['stim'].values,
                                                                                          b['stim'].values):
            raise ValueError(f'{fly}: mean_peak and baseline_std files differ in cells or trials')
        cells = p['cells'].values[sel].astype(int)
        stims = [str(s) for s in p['stim'].values]
        Fc = p['Fc'].values[:, sel]
        Z = p['Fc_zscore'].values[:, sel]
        B = b['Fc'].values[:, sel]
    blk = session_blocks(stims)
    labs = [stim_label(s, panel, concs) for s in stims]
    blk['label'] = [lab if lab is not None else (lambda n, c: f'{n} @ {c:g}')(*_split_stim(s))
                    for lab, s in zip(labs, stims)]
    blk['in_panel'] = [lab is not None for lab in labs]
    nT, nC = Fc.shape
    df = pd.DataFrame(dict(cell=np.tile(cells, nT), trial=np.repeat(np.arange(nT), nC)))
    for c in ('stim', 'label', 'in_panel', 'block', 'rep', 'presentation'):
        df[c] = np.repeat(blk[c].to_numpy(), nC)
    df['Fc'] = Fc.reshape(-1)
    df['z'] = Z.reshape(-1)
    df['baseline_sd'] = B.reshape(-1)
    return df


def reproduce_uploaded(fly: str, red: pd.DataFrame) -> dict:
    """Whether the uploaded table of one fly is the reduced file's `UPLOADED_RULE` on its cells: the cell sets, the
    largest absolute difference and the trials averaged per stimulus."""
    up = uploaded_table(fly)
    x = red[red.in_panel & (red.block == 0)]
    ours = x.groupby(['cell', 'label']).z.mean().unstack('label').reindex(columns=list(up.columns))
    same_cells = set(up.index.astype(int)) == set(ours.index.astype(int))
    ours = ours.reindex(up.index.astype(int))
    diff = float(np.nanmax(np.abs(ours.to_numpy() - up.to_numpy()))) if same_cells else float('nan')
    nper = x.drop_duplicates(['trial']).groupby('label').size().to_dict()
    return dict(n_cells_uploaded=int(len(up)), n_cells_ours=int(red.cell.nunique()), same_cells=bool(same_cells),
                max_abs_diff=diff, n_nan_uploaded=int(up.isna().sum().sum()), trials_per_stimulus=nper,
                cells=CELL_RULE, trials=UPLOADED_RULE, passed=bool(same_cells and np.isfinite(diff) and diff <= REPRO_TOL))


def noise_values(red: pd.DataFrame, factor: float = 1.0) -> np.ndarray:
    """The reduced rows' responses in baseline-noise units times `factor`: Fc over the root mean square of the cell's
    per-trial baseline SD over the session's trials."""
    rms = red.groupby('cell').baseline_sd.apply(lambda s: float(np.sqrt(np.nanmean(np.square(s)))))
    return red.Fc.to_numpy(dtype=float) / red.cell.map(rms).to_numpy(dtype=float) * float(factor)


def fly_bridge(red: pd.DataFrame) -> float:
    """A fly's noise bridge: the median over its cells of the trial-noise SD (`robust_trial_sd`) of the panel's
    stimuli's first blocks, in baseline-noise units (before the footing)."""
    from data_prep import stage_population_data as D
    sel = (red.in_panel & (red.block == 0)).to_numpy()
    x = red[sel].assign(v=noise_values(red)[sel])      # the baseline SD's RMS over the whole session
    labs = sorted(x.label.unique())
    cells = sorted(x.cell.unique())
    X = np.full((int(x.rep.max()) + 1, len(cells), len(labs)), np.nan)
    ci = {c: i for i, c in enumerate(cells)}
    li = {o: i for i, o in enumerate(labs)}
    X[x.rep.to_numpy(), x.cell.map(ci).to_numpy(), x.label.map(li).to_numpy()] = x.v.to_numpy()
    return float(np.nanmedian(D.robust_trial_sd(X)))


def e134_lock() -> dict:
    """The committed E134 lock (experiments/E134_lock.json)."""
    from data_prep import stage_population_data as D
    L = D.read_lock()
    if L is None:
        raise FileNotFoundError(D.LOCK_COPY)
    return L


def criteria_from_e134(L: Optional[Mapping] = None) -> dict:
    """The units' criteria and the footing's reference, from the E134 lock: 'noise' megamat's soma level, 'z' the
    mixtures' session-z criterion, and megamat's median noise bridge."""
    L = e134_lock() if L is None else L
    sc = L['soma_criteria']
    return dict(noise=float(sc['megamat']['criterion']), z=float(sc['mixtures']['z_criterion']),
                megamat_bridge=float(sc['noise_bridge_reference']['median']))


def reading_definitions() -> dict:
    """The readings as the lock records them (the module docstring)."""
    return {
        'as_uploaded': dict(components='first block', mixtures='first block', units='z',
                            rule="olfdb dataset 8 as uploaded; the single trials of each stimulus's first block for the "
                                 'noise'),
        'first_presentations': dict(components='first presentation', mixtures='first presentation', units='noise',
                                    rule="every stimulus's first presentation of the session; the noise SD from the "
                                         'first blocks, n = 1'),
        'components_first_block': dict(components='first block', mixtures='first block', units='noise',
                                       rule="every stimulus's first block (the uploaded table's trials)"),
        'components_second_block': dict(components='second block (first where none)', mixtures='first block',
                                        units='noise', rule="the components' second block where there is one, else "
                                                            'their first; the mixtures\' block'),
        'components_both_blocks': dict(components='both blocks (first where one)', mixtures='first block',
                                       units='noise', rule="the components' two blocks together where there are two, "
                                                           "else the first; the mixtures' block"),
    }


def part_lock() -> dict:
    """Design items 1 and 3: the eleven flies' files, stimulus orders and reproductions, the reduced per-trial files,
    the footing, the criteria and the readings; lock.json and lock_flies.csv (results/), and the committed copy
    experiments/E137_lock.json (commit it before any score)."""
    from kcq import harness as H
    t0 = time.time()
    crit = criteria_from_e134()
    data_dir().mkdir(parents=True, exist_ok=True)
    flies, rows, missing, bridges = {}, [], [], {}
    for r in mixture_flies().itertuples():
        f = find_fly_files(r.fly, r.panel)
        rec = dict(panel=r.panel, n_cells_uploaded=int(r.n_cells), session=f['session'], movie_type=f['movie_type'])
        if f['mean_peak'] is None or f['baseline_std'] is None:
            rec.update(found=False, readings=['as_uploaded'], files={})
            missing.append(dict(fly=r.fly, panel=r.panel, reason='no mean_peak / baseline_std file found'))
            flies[r.fly] = rec
            rows.append(dict(fly=r.fly, panel=r.panel, found=False))
            print(f'{r.fly}: files not found', flush=True)
            continue
        red = reduce_fly(r.fly, r.panel, f)
        date, num = r.fly.split('/')
        rp = data_dir() / f"{date}__fly{int(num):02d}__{f['session']}.parquet"
        red.to_parquet(rp, index=False)
        rep = reproduce_uploaded(r.fly, red)
        trials = red.drop_duplicates('trial').sort_values('trial')
        order = [str(s) for s in trials.stim]
        cfg = config_order(f['config']) if f['config'] is not None else None
        gen = generator_settings(f['generator'])
        blocks = trials[trials.in_panel].groupby('label').block.max().add(1).to_dict()
        files = {k: file_entry(f[k]) for k in ('mean_peak', 'baseline_std', 'tensor_mean_peak', 'config', 'generator')
                 if f[k] is not None}
        files['reduced'] = file_entry(rp)
        bridges[r.fly] = fly_bridge(red)
        rec.update(found=True, files=files, folder=_rel(f['folder']), olf_config_name=f['olf_config_name'],
                   order=order, labels=[str(x) for x in trials.label], blocks_per_stimulus={k: int(v) for k, v in
                                                                                             blocks.items()},
                   n_component_blocks=int(min(v for k, v in blocks.items() if 'mix' not in k)),
                   config_order_matches=(None if cfg is None else orders_match(cfg, order)),
                   config_order=cfg,
                   generator=gen, reproduction=rep, readings=list(READINGS), noise_bridge=bridges[r.fly],
                   n_trials=int(len(trials)))
        flies[r.fly] = rec
        rows.append(dict(fly=r.fly, panel=r.panel, found=True, session=f['session'], n_trials=len(trials),
                         n_component_blocks=rec['n_component_blocks'], config_order_matches=rec['config_order_matches'],
                         randomize_presentation_order=gen.get('randomize_presentation_order'),
                         reproduced=rep['passed'], max_abs_diff=rep['max_abs_diff'], n_cells=rep['n_cells_ours'],
                         noise_bridge=bridges[r.fly], reduced_sha256=files['reduced']['sha256'],
                         mean_peak=files['mean_peak']['path'], mean_peak_sha256=files['mean_peak']['sha256'],
                         baseline_std_sha256=files['baseline_std']['sha256']))
        print(f"{r.fly} {r.panel}: {f['session']}, reproduced {rep['passed']} ({rep['max_abs_diff']:.2e}), "
              f"{rec['n_component_blocks']} component blocks, {time.time() - t0:.0f} s", flush=True)
    group = float(np.median(list(bridges.values()))) if bridges else float('nan')
    footing = dict(rule="the 'noise' units times megamat's median noise bridge over the median of the mixture flies' "
                        'bridges (one group: FOOTING_GROUPS)', groups=FOOTING_GROUPS, megamat_bridge=crit['megamat_bridge'],
                   group_bridge=group, factor=crit['megamat_bridge'] / group if bridges else None,
                   per_fly_bridge=bridges,
                   per_panel_bridge={p: float(np.median([bridges[k] for k, v in flies.items() if v['panel'] == p
                                                         and k in bridges])) for p in MIXTURE_PANELS})
    lock = dict(experiment=EXPERIMENT, part='lock', host=_host(), time=_now(), commits=H.code_commits(),
                e134_lock=file_entry(LOCK_COPY.parent / 'E134_lock.json'),
                olfdb_dataset_8=file_entry(_e134_ds8_path()),
                flies=flies, missing=missing, readings=reading_definitions(), primary_units=PRIMARY_UNITS,
                units=dict(noise="baseline-noise units on megamat's footing: Remy's mean_peak Fc (2 to 8 s mean minus "
                                 "-5 to 0 s mean) over the RMS over the session's trials of her per-trial baseline SD "
                                 '(baseline_std file), times the footing factor',
                           z="Remy's session z (Fc_zscore), as uploaded"),
                criteria=dict(noise=crit['noise'], z=crit['z'],
                              source="the E134 lock's soma_criteria (megamat's level; the mixtures' z bridge)"),
                footing=footing,
                noise_rule="the flies' single-trial residual SD about each cell's block mean, pooled over the panel's "
                           "flies, cells, stimuli and the reading's blocks (n / (n - 1)); added to a model value over "
                           'sqrt(n), n the trials the reading averages for the stimulus',
                classes=dict(stage_rule="inside if the family's strength error is at or below the flies' most atypical "
                                        'fly on strength on 6 of 8 wirings; far if above 2 on 6 of 8; else between',
                             readable="a mixture panel is readable on strength if the family's class is the same under "
                                      'the four position readings',
                             replaces='a reading replaces the lock\'s where the panel is readable under it and the '
                                      'reproduction holds',
                             v2_readable='validation2 is readable if the class is the same on the six flies and on all '
                                         'twenty under both rules'),
                v2_rules={k: V2_RULE_TO_P[k] for k in V2_RULES}, v2_anchor=list(V2_COMMON_ODORS))
    _write_json(out_dir() / 'lock.json', lock)
    _write_json(LOCK_COPY, lock)
    pd.DataFrame(rows).to_csv(out_dir() / 'lock_flies.csv', index=False)
    (data_dir() / 'README.md').write_text(_data_readme(flies))
    print(f'lock: {len(flies)} flies, {len(missing)} missing, {time.time() - t0:.0f} s; commit {LOCK_COPY}', flush=True)
    return lock


def _e134_ds8_path() -> Path:
    from data_prep import stage_population_data as D
    return D.olfdb_path(D.OLFDB_DATASETS['kc_mixtures'])


def _data_readme(flies: Mapping) -> str:
    L = ['# Remy\'s kiwi / control mixture KC somata, per trial (E137)', '',
         'One parquet file per fly of olfdb dataset 8 (the eleven mixture flies), written by '
         '`experiments/E137_heldout_readings.py --part lock` (ejhonglab/model_experiment_queue) from Remy\'s processed '
         'suite2p respvec files on matrix (natural_mixtures/processed_data/<date>/<fly>/<session>/.../suite2p/combined/'
         f'{MEAN_PEAK} and {BASELINE_STD}), read only (her processed data may be read and copied, 2026-11-20).', '',
         f'Cells: {CELL_RULE}, the cells the uploaded table holds. Every trial of the session.', '',
         'Columns: cell (suite2p cell index, olfdb\'s roi), trial (index in the session), stim (Remy\'s name), label '
         "(the panel's name: EtOH, IAol, ..., kmix-2, kmix-1, kmix0; else '<name> @ <conc>'), in_panel, block (the "
         "stimulus's run of consecutive presentations, 0 first), rep (presentation within the block), presentation "
         '(within the session), Fc (her mean_peak response, the 2 to 8 s mean minus the -5 to 0 s mean of Fc), z '
         '(Fc_zscore, her session z), baseline_sd (her per-trial baseline SD of Fc).', '',
         'The sha256 of every file, and of its sources, is in results/E137_heldout_readings/lock.json.', '',
         '| fly | panel | session | file |', '|---|---|---|---|']
    for k, v in flies.items():
        if v.get('found'):
            L.append(f"| {k} | {v['panel']} | {v['session']} | {Path(v['files']['reduced']['path']).name} |")
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------------------
# the public interface (E138 reads it)

def load_lock(path: Optional[Path] = None) -> dict:
    """The committed lock (experiments/E137_lock.json, or `path`): 'flies' maps each mixture fly ('<date>/<fly_num>')
    to its panel, its files (path, bytes, sha256; 'reduced' the per-trial file the readings read), its session's
    stimulus order, whether the uploaded table is reproduced from them (and which cells, trials and blocks it averages),
    and the readings it is read under ('as_uploaded' only for a fly whose files are missing); 'readings' defines each
    of `READINGS`; 'criteria' the response criterion of each unit; 'footing' the noise units' factor. Raises
    FileNotFoundError before the lock is written, and PermissionError if the copy (the default path) has uncommitted
    changes."""
    p = LOCK_COPY if path is None else Path(path)
    if not p.exists():
        raise FileNotFoundError(f'{p}: run --part lock first')
    if path is None and not _committed(p):
        raise PermissionError(f'{p}: commit the lock before any score')
    return json.loads(p.read_text())


@lru_cache(maxsize=16)
def _reduced(path: str, sha: str) -> pd.DataFrame:
    p = _abs(path)
    got = sha256_file(p)
    if got != sha:
        raise ValueError(f"{p.name}: sha256 {got[:12]} differs from the lock's {sha[:12]}")
    return pd.read_parquet(p)


def fly_reduced(fly: str, lock: Optional[Mapping] = None) -> Optional[pd.DataFrame]:
    """One fly's reduced per-trial rows (the lock's file, sha256 checked); None for a fly whose files are missing."""
    lock = load_lock() if lock is None else lock
    rec = lock['flies'][fly]
    if not rec.get('found'):
        return None
    e = rec['files']['reduced']
    return _reduced(e['path'], e['sha256'])


def select_trials(red: pd.DataFrame, reading: str) -> pd.DataFrame:
    """The panel's trials a reading averages (module docstring), with a column 'noise_block' marking the blocks its
    trial noise is estimated on."""
    x = red[red.in_panel].copy()
    comp = ~x.label.str.contains('mix')
    nb = x.groupby('label').block.transform('max')
    if reading in ('as_uploaded', 'components_first_block'):
        keep = x.block == 0
        nblock = keep
    elif reading == 'first_presentations':
        keep = (x.block == 0) & (x.rep == 0)
        nblock = x.block == 0
    elif reading == 'components_second_block':
        keep = np.where(comp & (nb >= 1), x.block == 1, x.block == 0)
        nblock = keep
    elif reading == 'components_both_blocks':
        keep = np.where(comp & (nb >= 1), x.block <= 1, x.block == 0)
        nblock = keep
    else:
        raise ValueError(f'{reading=}')
    x['use'] = np.asarray(keep, dtype=bool)
    x['noise_block'] = np.asarray(nblock, dtype=bool)
    return x


def fly_values(red: pd.DataFrame, units: str, factor: float) -> np.ndarray:
    """The reduced rows' values in `units` ('noise': `noise_values` times the footing factor; 'z': session z)."""
    if units == 'noise':
        return noise_values(red, factor)
    if units == 'z':
        return red.z.to_numpy(dtype=float)
    raise ValueError(f'{units=}')


def reading_fly(fly: str, reading: str, units: str, lock: Mapping) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """One fly under a reading: (cells x the panel's stimuli trial means, cells x (stimulus, repeat) single trials of
    the means, {'ss', 'dof', 'n': {stimulus: trials averaged}} of the trial noise on the reading's noise blocks). A fly
    read as uploaded (files missing, or 'as_uploaded' in 'z') takes its trial means from olfdb dataset 8."""
    from data_prep import stage_population_data as D
    panel = lock['flies'][fly]['panel']
    odors = list(D.panel_odors(panel))
    red = fly_reduced(fly, lock)
    if red is None:
        up = uploaded_table(fly)
        up.index = up.index.astype(str)
        return up, None, dict(ss=0.0, dof=0, n={})
    x = select_trials(red, reading)
    x['v'] = fly_values(red, units, lock['footing']['factor'])[red.in_panel.to_numpy()]
    nb = x[x.noise_block]
    mu = nb.groupby(['cell', 'label', 'block']).v.transform('mean')
    cnt = nb.groupby(['cell', 'label', 'block']).v.transform('count')
    r = (nb.v - mu).to_numpy(dtype=float)
    ok = cnt.to_numpy() >= 2
    ss = float(np.nansum(r[ok] ** 2))
    dof = int(nb[ok].groupby(['cell', 'label', 'block']).size().sub(1).sum())
    u = x[x.use]
    means = u.groupby(['cell', 'label']).v.mean().unstack('label').reindex(columns=odors)
    means.index = means.index.astype(str)
    if reading == 'as_uploaded' and units == 'z':
        up = uploaded_table(fly)
        up.index = up.index.astype(str)
        means = up.reindex(index=means.index)
    u = u.assign(r=u.groupby(['cell', 'label']).cumcount())
    tr = u.pivot_table(index='cell', columns=['label', 'r'], values='v', aggfunc='first')
    tr.index = tr.index.astype(str)
    tr.columns = tr.columns.set_names(['odor', 'repeat'])
    n = u.drop_duplicates(['label', 'trial']).groupby('label').size().to_dict()
    return means, tr.reindex(index=means.index), dict(ss=ss, dof=dof, n={k: int(v) for k, v in n.items()})


@dataclass
class TrialNoise:
    """The flies' trial noise of one panel under a reading (module docstring): the single-trial SD and, per stimulus,
    the number of trials a mean averages."""
    sd: float
    n: Dict[str, int]
    dof: int = 0
    flies: List[str] = field(default_factory=list)

    def add(self, values: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
        """`values` (units x stimuli) plus independent Gaussian noise of SD sd / sqrt(n[stimulus])."""
        X = values.to_numpy(dtype=float)
        k = np.array([max(int(self.n.get(str(c), 1)), 1) for c in values.columns], dtype=float)
        return pd.DataFrame(X + rng.standard_normal(X.shape) * (self.sd / np.sqrt(k))[None, :], index=values.index,
                            columns=values.columns)

    def as_dict(self) -> dict:
        return dict(sd=float(self.sd), n=dict(self.n), dof=int(self.dof), flies=list(self.flies))


def reading_population(panel: str, reading: str, units: Optional[str] = None, lock: Optional[Mapping] = None):
    """The flies of one panel under a reading as data_prep/stage_population_data.py's `Population` (tables, trials,
    criterion, units), and their `TrialNoise` (None where no fly has trials)."""
    from data_prep import stage_population_data as D
    lock = load_lock() if lock is None else lock
    units = PRIMARY_UNITS[reading] if units is None else units
    flies = sorted(k for k, v in lock['flies'].items() if v['panel'] == panel)
    tables, trials, ss, dof, ns = {}, {}, 0.0, 0, {}
    for f in flies:
        m, t, nz = reading_fly(f, reading, units, lock)
        tables[f] = m
        if t is not None:
            trials[f] = t
        ss += nz['ss']
        dof += nz['dof']
        for k, v in nz['n'].items():
            ns.setdefault(k, []).append(v)
    crit = float(lock['criteria'][units])
    pop = D.Population(stage='kc_soma', panel=panel, cohort='wild_type', odors=tuple(D.panel_odors(panel)),
                       flies=tuple(flies), tables=tables, criterion=crit,
                       units=D.UNITS_BASELINE if units == 'noise' else D.UNITS_SESSION_Z, graded=True,
                       trials=trials or None, has_trials=bool(trials), unit_kind='kc', named_units=False,
                       recording_site=D.RECORDING_SITE['kc_soma'],
                       notes=[f'E137 reading {reading} in {units} units'])
    noise = (TrialNoise(sd=float(np.sqrt(ss / dof)), n={k: int(np.median(v)) for k, v in ns.items()}, dof=dof,
                        flies=list(trials)) if dof > 0 else None)
    return pop, noise


def fly_mixture_view(panel: str, reading: str, trials: bool = True, lock: Optional[Mapping] = None,
                     units: Optional[str] = None):
    """The flies' KC somata of one mixture panel ('kmix' or 'cmix') under a named reading (one of `READINGS`), as the
    re-score reads them: a list of `kcq.population_fit.Population`, one per fly (named '<date>/<fly_num>'), cells x the
    panel's stimuli (data_prep.stage_population_data's KMIX_ODORS / CMIX_ODORS), in the reading's units
    (`PRIMARY_UNITS`, or `units`: 'noise' or 'z') and at the lock's criterion; with `trials`, a tuple (that list,
    {fly: DataFrame of single trials with columns (odor, repeat), None for a fly read as uploaded}). Reads only the
    lock's files (sha256 checked)."""
    from kcq import population_fit as P
    pop, _ = reading_population(panel, reading, units, lock)
    fps = P.fly_populations(pop)
    if not trials:
        return fps
    return fps, {f: (pop.trials or {}).get(f) for f in pop.flies}


# ---------------------------------------------------------------------------
# design item 2: the flies by position

def part_position() -> pd.DataFrame:
    """Per fly and stimulus (the panel's, the solvent and the probes), each presentation's response rate (single trial
    above the criterion) and mean response, and each block's (its mean above the criterion), in both units;
    position_trials.csv, position_blocks.csv, position_summary.csv (the third over the first within the first block,
    the second block over the first, the session's end over its start for the solvent and the probes) and
    position_decline.csv (per fly, the median of each ratio over the components, the mixtures, the solvent and the
    probes: whether the decline is general or specific to the components)."""
    lock = load_lock()
    rows, brows = [], []
    for fly, rec in lock['flies'].items():
        red = fly_reduced(fly, lock)
        if red is None:
            continue
        for units in UNITS:
            v = fly_values(red, units, lock['footing']['factor'])
            crit = float(lock['criteria'][units])
            x = red.assign(v=v)
            g = x.groupby(['trial', 'stim', 'label', 'in_panel', 'block', 'rep', 'presentation'])
            t = g.v.agg(rate=lambda s: float((s > crit).mean()), mean='mean').reset_index()
            t.insert(0, 'units', units)
            t.insert(0, 'panel', rec['panel'])
            t.insert(0, 'fly', fly)
            rows.append(t)
            bm = x.groupby(['label', 'block', 'cell']).v.mean()
            b = bm.groupby(['label', 'block']).agg(rate=lambda s: float((s > crit).mean()), mean='mean').reset_index()
            b['n_trials'] = x.drop_duplicates('trial').groupby(['label', 'block']).size().reindex(
                pd.MultiIndex.from_frame(b[['label', 'block']])).to_numpy()
            b.insert(0, 'units', units)
            b.insert(0, 'panel', rec['panel'])
            b.insert(0, 'fly', fly)
            brows.append(b)
    T = pd.concat(rows, ignore_index=True)
    B = pd.concat(brows, ignore_index=True)
    T.to_csv(out_dir() / 'position_trials.csv', index=False)
    B.to_csv(out_dir() / 'position_blocks.csv', index=False)
    S = position_summary(T, B)
    S.to_csv(out_dir() / 'position_summary.csv', index=False)
    D_ = position_decline(S)
    D_.to_csv(out_dir() / 'position_decline.csv', index=False)
    return S


def stim_kind(label: str) -> str:
    """'mixture', 'component', 'solvent', 'probe' or 'other' (the partial mixtures)."""
    name = str(label).split(' @ ')[0]
    if name.startswith(('kmix', 'cmix')) and ' @ ' not in str(label):
        return 'mixture'
    if ' @ ' not in str(label):
        return 'component'
    if name == SOLVENT:
        return 'solvent'
    if name in PROBES:
        return 'probe'
    return 'other'


def position_summary(T: pd.DataFrame, B: pd.DataFrame) -> pd.DataFrame:
    """Per fly, units and stimulus: the first block's third over first presentation (rate and mean), the second block
    over the first (block means), and the last block over the first (the solvent and probes at the end over the
    start)."""
    rows = []
    for (fly, units, lab), t in T.groupby(['fly', 'units', 'label'], sort=False):
        b = B[(B.fly == fly) & (B.units == units) & (B.label == lab)].set_index('block')
        t0 = t[t.block == 0].set_index('rep')
        rec = dict(fly=fly, panel=t.panel.iloc[0], units=units, label=lab, kind=stim_kind(lab),
                   n_blocks=int(b.index.max()) + 1, first_trial=int(t.trial.min()), last_trial=int(t.trial.max()))
        for q in ('rate', 'mean'):
            a, c = t0[q].get(0, np.nan), t0[q].get(int(t0.index.max()), np.nan)
            rec[f'{q}_first'] = float(a)
            rec[f'{q}_third_over_first'] = float(c / a) if a and np.isfinite(a) and a != 0 else np.nan
            rec[f'{q}_block0'] = float(b[q].get(0, np.nan))
            rec[f'{q}_block1'] = float(b[q].get(1, np.nan))
            rec[f'{q}_block1_over_block0'] = (float(b[q][1] / b[q][0]) if 1 in b.index and b[q][0] != 0 else np.nan)
            last = int(b.index.max())
            rec[f'{q}_last_over_first_block'] = (float(b[q][last] / b[q][0]) if last > 0 and b[q][0] != 0 else np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def position_decline(S: pd.DataFrame) -> pd.DataFrame:
    """Per fly and units, the median over each kind of stimulus of the within-block ratio and of the later block over
    the first (solvent and probes: end over start), beside the components'."""
    rows = []
    for (fly, units), s in S.groupby(['fly', 'units'], sort=False):
        rec = dict(fly=fly, panel=s.panel.iloc[0], units=units)
        for kind in ('component', 'mixture', 'solvent', 'probe'):
            k = s[s.kind == kind]
            for q in ('rate', 'mean'):
                rec[f'{kind}_{q}_third_over_first'] = float(k[f'{q}_third_over_first'].median()) if len(k) else np.nan
                rec[f'{kind}_{q}_later_over_first_block'] = (float(k[f'{q}_last_over_first_block'].median())
                                                             if len(k) else np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# design items 3 and 4: the readings of the flies

def fly_reference(panel: str, reading: str, units: str, lock: Mapping, odors: Optional[Sequence[str]] = None):
    """(the flies' `Population` and `TrialNoise`, their hal's populations, `FlyReference` on `odors` (default the
    panel's), its breadth bins)."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    pop, noise = reading_population(panel, reading, units, lock)
    fps = P.fly_populations(pop)
    od = list(D.panel_odors(panel)) if odors is None else [o for o in odors if o in D.panel_odors(panel)]
    ref = P.fly_reference(fps, od, D.breadth_bins(len(od)))
    return pop, noise, fps, ref


def flies_summary(fps, ref, n_floor: int = N_FLOOR) -> dict:
    """The flies' floor (mean of `P.floor`'s draws) and one recording against the others (mean and largest), on every
    view and the population error."""
    from kcq import population_fit as P
    fl = P.floor_summary(P.floor(ref, n_boot=n_floor))
    loo = P.leave_one_out(fps, ref=ref)
    out = {}
    for c in ('error', 'view_strength', 'view_sharing', 'view_breadth'):
        out[f'floor_{c}'] = float(fl[c])
        out[f'loo_mean_{c}'] = float(loo[c].mean())
        out[f'loo_max_{c}'] = float(loo[c].max())
    out['loo_max_fly'] = str(loo['view_strength'].idxmax())
    return out


def part_readings() -> pd.DataFrame:
    """Each reading's flies per panel and units: the floor and one recording against the others on every view, the
    response rates (mean, least and largest over flies), the trial noise; readings_flies.csv, readings_rates.csv."""
    from data_prep import stage_population_data as D
    lock = load_lock()
    rows, rrows = [], []
    for panel in MIXTURE_PANELS:
        for reading in READINGS:
            for units in UNITS:
                pop, noise, fps, ref = fly_reference(panel, reading, units, lock)
                rates = {f: D.fly_rate(pop.tables[f], pop.criterion) for f in pop.flies}
                rec = dict(panel=panel, reading=reading, units=units, primary=PRIMARY_UNITS[reading] == units,
                           n_flies=len(fps), criterion=pop.criterion, rate_mean=float(np.mean(list(rates.values()))),
                           rate_min=float(np.min(list(rates.values()))), rate_max=float(np.max(list(rates.values()))),
                           noise_sd=np.nan if noise is None else noise.sd,
                           noise_n=None if noise is None else json.dumps(noise.n), **flies_summary(fps, ref))
                rows.append(rec)
                for f, r in rates.items():
                    per = (pop.tables[f] > pop.criterion).mean(axis=0)
                    rrows.append(dict(panel=panel, reading=reading, units=units, fly=f, rate=r,
                                      **{f'rate_{o}': float(per[o]) for o in per.index}))
                print(f'readings {panel} {reading} {units}: rate {rec["rate_mean"]:.4f}, floor strength '
                      f'{rec["floor_view_strength"]:.3f}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'readings_flies.csv', index=False)
    pd.DataFrame(rrows).to_csv(out_dir() / 'readings_rates.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# design items 4 and 5: the re-score

def mixture_groups() -> Dict[str, list]:
    """E134's saved mixture runs by group: 'fitted' and 'held' (the family's eight wirings, in E134's order) and each
    beside model ('<model>', one run)."""
    from kcq import stage_readings as R
    out = {}
    for g in FAMILY_GROUPS:
        out[g] = [R.load_step3_run(f'mixtures_{g}_{w}') for w in FAMILY_WIRINGS]
    for m in BESIDE_MODELS:
        try:
            out[m] = [R.load_step3_run(f'mixtures_fitted_{m}')]
        except FileNotFoundError:
            print(f'rescore: mixtures_fitted_{m} not found', flush=True)
    return out


def read_model_panel(runs: Sequence, panel: str, flies, noise: Optional[TrialNoise]) -> list:
    """Each run's KC somata on a mixture panel read as E134 reads them (`stage_readings.DEFAULT_READINGS['kc_soma']`:
    the APL-on spike counts, the quantile curve set per wiring on these flies' pooled trial means, the flies'
    criterion), with `noise` added after the curve (seed NOISE_SEED + the run's index; none where None): one hal's
    `Population` per run."""
    from kcq import population_fit as P
    from kcq import stage_readings as R
    rd = R.DEFAULT_READINGS['kc_soma']
    pops = []
    for i, run in enumerate(runs):
        odors = dict(R.run_panels(run, 'kc_soma'))[panel]
        U = R.sample_units(run, 'kc_soma', 'on', odors, rd, i)
        curve = R.set_curve(R.pooled_values([U]), flies, odors, rd, i)
        V = R.apply_curve(U, curve)
        if noise is not None:
            V = noise.add(V, np.random.default_rng(R.NOISE_SEED + i))
        pops.append(P.Population.from_criterion(V, float(flies.criterion), at_or_above=False,
                                                name=f'{run.model}/{run.wiring}', named_units=False, graded=True))
    return pops


def flat_views(odors: Sequence[str]) -> Dict[str, pd.Series]:
    """The flat profile: every odor at the panel's mean (both strength parts 1 at every odor)."""
    s = pd.Series(1.0, index=list(odors))
    return {'strength/rate': s, 'strength/mean': s.copy()}


def flat_score(ref) -> dict:
    """The flat profile scored like a model's mean view (strength only)."""
    from kcq import population_fit as P
    sc = P.score_views(flat_views(ref.odors), ref)
    return dict(view_strength=float(sc.views.get('strength', np.nan)),
                part_strength_rate=float(sc.parts.get('strength/rate', np.nan)),
                part_strength_mean=float(sc.parts.get('strength/mean', np.nan)))


def orn_prediction(panel: str, ref) -> dict:
    """The ORN recordings of the panel as the prediction (hal's `no_step_score`): their mean views scored against the
    KC flies."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    orn = P.fly_populations(D.load_population('orn', panel))
    sc = P.no_step_score(orn, ref)
    return dict(view_strength=float(sc.views.get('strength', np.nan)), view_sharing=float(sc.views.get('sharing',
                np.nan)), view_breadth=float(sc.views.get('breadth', np.nan)), error=float(sc.error))


def _split(pops, ref) -> pd.DataFrame:
    from experiments.E136_stage_population_fit import strength_split
    return strength_split(pops, ref)


def e134_check(groups: Mapping[str, list]) -> dict:
    """The lock's as-uploaded flies (E134's population file) with no noise reproduce E134's mixture scores
    (score_scores.csv, groups mixtures_fitted / mixtures_held, the default reading): the largest absolute difference."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    from kcq import stage_readings as R
    p = shared() / 'results' / 'E134_population_baseline' / 'score_scores.csv'
    if not p.exists():
        return dict(checked=False)
    s = pd.read_csv(p)
    rd = R.DEFAULT_READINGS['kc_soma'].name
    out = {}
    for g in FAMILY_GROUPS:
        for panel in MIXTURE_PANELS:
            fl = D.load_population('kc_soma', panel)
            pops = read_model_panel(groups[g], panel, fl, None)
            fps = P.fly_populations(fl)
            ref = P.fly_reference(fps, list(D.panel_odors(panel)), D.breadth_bins(len(D.panel_odors(panel))))
            tab = P.wirings_as_flies(pops, ref)
            x = s[(s.group == f'mixtures_{g}') & (s.stage == 'kc_soma') & (s.panel == panel) & (s.condition == 'on')
                  & (s.reading == rd)].set_index('wiring')
            names = [str(n).partition('/')[2] or str(n) for n in tab.index]
            t = tab.set_axis(names)
            d = [abs(float(t.loc[w, c]) - float(x.loc[w, c])) for w in x.index if w in t.index
                 for c in ('error', 'view_strength', 'view_sharing', 'view_breadth')]
            out[f'{g}/{panel}'] = dict(max_abs_diff=float(np.max(d)) if d else np.nan, n=len(d),
                                       mean_strength=float(t.loc['mean', 'view_strength']))
    out['checked'] = True
    out['passed'] = bool(all(v['max_abs_diff'] <= 1e-9 for k, v in out.items() if isinstance(v, dict)))
    return out


def part_rescore() -> pd.DataFrame:
    """E134's saved mixture runs re-read with the flies' trial noise and scored under each reading and units, beside the
    flat profile and the ORN recordings; every strength error split (E136's `strength_split`); rescore_scores.csv,
    rescore_split.csv, rescore_beside.csv, rescore_e134_check.json."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    lock = load_lock()
    t0 = time.time()
    groups = mixture_groups()
    chk = e134_check(groups)
    _write_json(out_dir() / 'rescore_e134_check.json', chk)
    print(f'rescore: E134 check {chk}', flush=True)
    rows, srows, brows = [], [], []
    for panel in MIXTURE_PANELS:
        for reading in READINGS:
            for units in UNITS:
                for with_noise in (True, False):
                    pop, noise, fps, ref = fly_reference(panel, reading, units, lock)
                    fs = flies_summary(fps, ref)
                    base = dict(panel=panel, reading=reading, units=units, primary=PRIMARY_UNITS[reading] == units,
                                noise=bool(with_noise and noise is not None), n_flies=len(fps), **fs)
                    if with_noise:
                        brows.append(dict(base, kind='flat', **flat_score(ref)))
                        brows.append(dict(base, kind='orn_recordings', **orn_prediction(panel, ref)))
                    for g, runs in groups.items():
                        pops = read_model_panel(runs, panel, pop, noise if with_noise else None)
                        tab = P.wirings_as_flies(pops, ref) if len(pops) > 1 else P.score_many(pops, ref)
                        for name, r in tab.iterrows():
                            w = str(name).partition('/')[2] or str(name)
                            rates = (float(np.nanmean([p.responds.to_numpy().mean() for p in pops]))
                                     if w in ('mean', 'common') else
                                     float(next(p.responds.to_numpy().mean() for p in pops if p.name == name)))
                            rows.append(dict(base, group=g, wiring=w, rate_model=rates, **r.to_dict()))
                        sp = _split(pops, ref)
                        for name, r in sp.iterrows():
                            srows.append(dict(base, group=g, wiring=str(name).partition('/')[2] or str(name),
                                              **r.to_dict()))
                print(f'rescore {panel} {reading} {units}: {time.time() - t0:.0f} s', flush=True)
    sc = pd.DataFrame(rows)
    sc.to_csv(out_dir() / 'rescore_scores.csv', index=False)
    pd.DataFrame(srows).to_csv(out_dir() / 'rescore_split.csv', index=False)
    pd.DataFrame(brows).to_csv(out_dir() / 'rescore_beside.csv', index=False)
    return sc


# ---------------------------------------------------------------------------
# design item 6: the rise

def _rise(table: pd.DataFrame, responds: pd.DataFrame, mix: str) -> dict:
    """Response rate and mean response at the mixture's 0 and -1 over -2, and the mixture at 0 over the components'
    mean, of one units x stimuli table."""
    rate, mean = responds.mean(axis=0), table.mean(axis=0)
    comps = [c for c in table.columns if 'mix' not in str(c)]
    out = {}
    for q, s in (('rate', rate), ('mean', mean)):
        lo = float(s.get(f'{mix}-2', np.nan))
        out[f'{q}_0_over_m2'] = float(s.get(f'{mix}0', np.nan)) / lo if lo else np.nan
        out[f'{q}_m1_over_m2'] = float(s.get(f'{mix}-1', np.nan)) / lo if lo else np.nan
        cm = float(s[comps].mean())
        out[f'{q}_mix0_over_components'] = float(s.get(f'{mix}0', np.nan)) / cm if cm else np.nan
    return out


def part_rise() -> pd.DataFrame:
    """The mixture's rise over its dilutions and the mixture over its components: per fly under each reading and units,
    the family's wirings and the beside models after the reading's curve and noise, and on raw spike counts (two or more
    spikes, APL on); rise.csv."""
    from kcq import stage_readings as R
    lock = load_lock()
    groups = mixture_groups()
    rows = []
    for panel in MIXTURE_PANELS:
        for reading in READINGS:
            for units in UNITS:
                pop, noise, fps, ref = fly_reference(panel, reading, units, lock)
                for f in fps:
                    rows.append(dict(panel=panel, reading=reading, units=units, side='flies', group='flies',
                                     wiring=f.name, **_rise(f.values, f.responds, panel)))
                for g, runs in groups.items():
                    for p in read_model_panel(runs, panel, pop, noise):
                        rows.append(dict(panel=panel, reading=reading, units=units, side='model', group=g,
                                         wiring=str(p.name).partition('/')[2], **_rise(p.values, p.responds, panel)))
        for g, runs in groups.items():
            for run in runs:
                odors = dict(R.run_panels(run, 'kc_soma'))[panel]
                c = R._abbrev(run.conditions['on']['spike_counts'])[list(odors)]
                rows.append(dict(panel=panel, reading='raw_counts', units='spikes', side='model', group=g,
                                 wiring=run.wiring, **_rise(c, (c >= R.MODEL_SPIKES).astype(float), panel)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'rise.csv', index=False)
    s = df.groupby(['panel', 'reading', 'units', 'side', 'group'])[[c for c in df if c.endswith(('m2', 'components'))]]
    s.median().reset_index().to_csv(out_dir() / 'rise_summary.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# design item 7: validation2 on all flies and odors

def v2_reference(flies, odors: Sequence[str], rule: str):
    """The flies' reference of validation2's KC somata on all flies, each fly on its own odors, under a partial-panel
    rule (`V2_RULES`, through hal's `kcq.population_fit.partial_reference`: 'common_mean' its 'anchor' rule on
    V2_COMMON_ODORS, 'fly_scale' its 'fly_scale')."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    kw = dict(anchor=list(V2_COMMON_ODORS)) if rule == 'common_mean' else {}
    return P.partial_reference(flies, rule=V2_RULE_TO_P[rule], odors=list(odors), breadth_edges=D.breadth_bins, **kw)


def v2_groups_of() -> Dict[str, str]:
    from data_prep import stage_population_data as D
    return D.v2_groups()


def part_validation2() -> pd.DataFrame:
    """Validation2's KC somata: the family (E107's saved runs, as fitted and held) and the beside models (E134's
    step-3 validation2 runs), read as E134 reads them (the curve and noise on all 20 flies), scored against the six
    flies on 15 odors (E134's reading), all twenty flies on all 22 odors under both rules, and each noise group's
    flies under both rules; the flat profile beside; validation2_scores.csv, validation2_flies.csv,
    validation2_classes.json."""
    from data_prep import stage_population_data as D
    from experiments import E134_population_baseline as E134
    from kcq import population_fit as P
    from kcq import stage_readings as R
    lock = R.require_lock()
    flies = R.FlyCache()
    fl = flies.get(('kc_soma', 'validation2', 'on'))
    groups = {g: [R.load_e107_run('validation2', g, w) for w in FAMILY_WIRINGS] for g in FAMILY_GROUPS}
    for m in BESIDE_MODELS:
        try:
            groups[m] = [R.load_step3_run(f'validation2_fitted_{m}')]
        except FileNotFoundError:
            pass
    pops = {}
    for g, runs in groups.items():
        fam = R.read_family(runs, stages=['kc_soma'], lock=lock, flies=flies)
        pops[g] = [r.populations[('kc_soma', 'validation2', 'on')] for r in fam.values()]
    fps_all = P.fly_populations(fl)
    six, od6 = E134.fly_side_read(fl)
    odors22 = [o for o in D.panel_odors('validation2') if o in pops['fitted'][0].values.columns]
    vg = v2_groups_of()
    refs = {('six', 'e134'): (six, P.fly_reference(six, od6, D.breadth_bins(len(od6))))}
    for rule in V2_RULES:
        refs[('all', rule)] = (fps_all, v2_reference(fps_all, odors22, rule))
        for grp in ('tsub03', 'other'):
            sub = [f for f in fps_all if vg.get(f.name) == grp]
            refs[(grp, rule)] = (sub, v2_reference(sub, odors22, rule))
    rows, frows = [], []
    for (fset, rule), (fps, ref) in refs.items():
        fs = flies_summary(fps, ref)
        base = dict(flies=fset, rule=rule, n_flies=len(fps), n_odors=len(ref.odors), **fs)
        frows.append(dict(base, **{f'flat_{k}': v for k, v in flat_score(ref).items()}))
        for g, ps in pops.items():
            tab = P.wirings_as_flies(ps, ref) if len(ps) > 1 else P.score_many(ps, ref)
            for name, r in tab.iterrows():
                rows.append(dict(base, group=g, wiring=str(name).partition('/')[2] or str(name), **r.to_dict()))
        print(f'validation2 {fset} {rule}: {len(fps)} flies, {len(ref.odors)} odors', flush=True)
    sc = pd.DataFrame(rows)
    sc.to_csv(out_dir() / 'validation2_scores.csv', index=False)
    pd.DataFrame(frows).to_csv(out_dir() / 'validation2_flies.csv', index=False)
    cls = classify_v2(sc)
    _write_json(out_dir() / 'validation2_classes.json', cls)
    return sc


def stage_class(x: pd.DataFrame, loo_max: float) -> str:
    """E134's stage rule on strength by odor over the family's wirings: 'inside' if at or below the flies' most atypical
    fly on 6 of 8, 'far' if above 2 on 6 of 8, else 'between'."""
    v = x.view_strength.to_numpy(dtype=float)
    if int((v <= loo_max).sum()) >= INSIDE_MIN_WIRINGS:
        return 'inside'
    if int((v > FAR_ERROR).sum()) >= INSIDE_MIN_WIRINGS:
        return 'far'
    return 'between'


def classify_v2(sc: pd.DataFrame) -> dict:
    """Validation2's class (fitted and held) on six flies and on all twenty under both rules; 'readable on strength' if
    the same in all three."""
    out = {}
    for g in FAMILY_GROUPS:
        cls, vals = {}, {}
        for (fset, rule), x in sc[(sc.group == g)].groupby(['flies', 'rule'], sort=False):
            w = x[x.wiring.isin(FAMILY_WIRINGS)]
            cls[f'{fset}/{rule}'] = stage_class(w, float(x.loo_max_view_strength.iloc[0]))
            m = x[x.wiring == 'mean'].iloc[0]
            vals[f'{fset}/{rule}'] = dict(strength=float(m.view_strength), sharing=float(m.view_sharing),
                                          breadth=float(m.view_breadth), error=float(m.error),
                                          loo_max_strength=float(m.loo_max_view_strength),
                                          floor_strength=float(m.floor_view_strength))
        keys = ['six/e134'] + [f'all/{r}' for r in V2_RULES]
        got = [cls.get(k) for k in keys]
        out[g] = dict(classes=cls, values=vals, readable=bool(None not in got and len(set(got)) == 1))
    return out


# ---------------------------------------------------------------------------
# design item 8: the order at the other panels

def order_share(trials: pd.DataFrame, crit: float) -> pd.DataFrame:
    """Per fly of long single-trial rows (fly, cell, odor, trial, value) with each trial's place in its block
    ('pos'), the share of strength by odor that the place carries: each trial's response rate (cells above `crit`)
    and mean response fitted per fly as odor effects plus a slope on the place (least squares); the odors' strength
    entries (each odor's mean over its trials, over the mean over odors) with and without the place's part; their RMS
    difference. Rows fly, part, slope, rms_shift (in strength entries), n_trials."""
    rows = []
    for fly, t in trials.groupby('fly', sort=False):
        per = t.groupby(['trial', 'odor', 'pos']).value.agg(rate=lambda s: float((s > crit).mean()), mean='mean')
        per = per.reset_index()
        for q in ('rate', 'mean'):
            y = per[q].to_numpy(dtype=float)
            X = pd.get_dummies(per.odor).to_numpy(dtype=float)
            pc = per.pos.to_numpy(dtype=float) - per.pos.mean()
            A = np.column_stack([X, pc])
            beta = np.linalg.lstsq(A, y, rcond=None)[0]
            slope = float(beta[-1])
            raw = per.assign(y=y).groupby('odor').y.mean()
            cor = per.assign(y=y - slope * pc).groupby('odor').y.mean()
            sr, scr = raw / raw.mean(), cor / cor.mean()
            rows.append(dict(fly=fly, part=f'strength/{q}', slope=slope, rms_shift=float(np.sqrt(np.mean(
                (sr - scr) ** 2))), n_trials=int(len(per)), n_odors=int(raw.size)))
    return pd.DataFrame(rows)


def megamat_order_trials() -> Tuple[pd.DataFrame, float]:
    """Remy's four megamat flies' single trials in baseline-noise units with each trial's place in its block (the
    session's blocks of the 17 odors, each once per block), and megamat's soma criterion."""
    from data_prep import stage_population_data as D
    r = pd.read_parquet(D.megamat_copy_dir() / 'frames_reduced.parquet')
    r = r[np.isfinite(r.sd) & (r.sd > 0)]
    r = r.assign(odor=r.stim.str.split(' @ ').str[0])
    r = r[r.odor.isin(D.MEGAMAT_ODORS)]
    tr = r.drop_duplicates(['fly', 'trial']).sort_values(['fly', 'trial'])
    tr = tr.assign(rep=tr.groupby(['fly', 'odor']).cumcount())
    tr['pos'] = tr.groupby(['fly', 'rep']).cumcount()
    r = r.merge(tr[['fly', 'trial', 'pos']], on=['fly', 'trial'])
    return pd.DataFrame(dict(fly=r.fly, cell=r.cell, odor=r.odor, trial=r.trial, pos=r.pos,
                             value=r.dfc / r.sd)), float(D.megamat_soma_criterion())


def v2_order_trials() -> Tuple[pd.DataFrame, float]:
    """Remy's validation2 flies' single trials (`v2_soma_long`'s units, before the footing, which is one factor per
    group and does not move strength by odor) with each trial's place in its repeat (the trials of one stim_occ in the
    acquisition's order), and megamat's criterion over the fly's footing factor."""
    import xarray as xr
    from data_prep import stage_population_data as D
    d = D.v2_copy_dir()
    fac = D.soma_factors(('kc_soma', 'validation2', 'wild_type'))
    drops = D.v2_qc_drops()
    rows = []
    for acq in D.v2_acqs():
        p = xr.open_dataset(d / acq / 'respvec' / 'xrds_peak_amp.nc')
        b = xr.open_dataset(d / acq / 'respvec' / 'xrds_baseline_std.nc')
        lab = f"{p.attrs['acq.date_imaged']}/{int(p.attrs['acq.fly_num'])}"
        mask = D.v2_mask(p[D.V2_CELL_MASK].values)
        amp = p['Fc'].values[:, mask]
        sd = np.sqrt(np.nanmean(b['Fc'].values[:, mask] ** 2, axis=0))
        names, stims = D._v2_stims(p['stim'].values)
        dropped = set(drops[drops.fly == lab].stim)
        occ = p['stim_occ'].values
        ok = np.isfinite(sd) & (sd > 0)
        pos_by = {}
        for k, (nm, st) in enumerate(zip(names, stims)):
            pos = pos_by.get(int(occ[k]), 0)          # the place among the acquisition's trials of this repeat
            pos_by[int(occ[k])] = pos + 1
            if nm is None or st in dropped:
                continue
            v = amp[k, ok] / sd[ok] * float(fac.get(lab, 1.0))
            rows.append(pd.DataFrame(dict(fly=lab, cell=np.arange(int(ok.sum())), odor=nm, trial=f'{acq}/{k}', pos=pos,
                                          value=v)))
        p.close()
        b.close()
    return pd.concat(rows, ignore_index=True), float(D.megamat_soma_criterion())


def part_order() -> dict:
    """The share of strength by odor an odor's place in its block carries at megamat and validation2 (`order_share`,
    beside the strength parts' one-fly units of E134's flies), and the APL16 recordings' order: Yang's processed tables
    (data/yang_apl16_kc_noise_units) hold each trial's odor and repeat but not its place in the session, which is a
    question for Yang; order_megamat.csv, order_validation2.csv, order.json."""
    from experiments import E134_flies as F
    from kcq import population_fit as P
    from data_prep import stage_population_data as D
    out = {}
    for panel, fn in (('megamat', megamat_order_trials), ('validation2', v2_order_trials)):
        t, crit = fn()
        sh = order_share(t, crit)
        sh.to_csv(out_dir() / f'order_{panel}.csv', index=False)
        fps, od = F.reading(('kc_soma', panel, 'wild_type'))     # E134's lock reading (validation2: six flies)
        ref = P.fly_reference(fps, od, D.breadth_bins(len(od)))
        units = {q: float(ref.one_fly[q]) for q in ('strength/rate', 'strength/mean')}
        out[panel] = dict(one_fly_units=units, rms_shift={q: float(sh[sh.part == q].rms_shift.mean())
                                                          for q in units},
                          share_of_unit={q: float(sh[sh.part == q].rms_shift.mean()) / units[q] for q in units},
                          slope_median={q: float(sh[sh.part == q].slope.median()) for q in units},
                          n_flies=int(sh.fly.nunique()))
    y = pd.read_parquet(shared() / 'data' / 'yang_apl16_kc_noise_units' / 'per_trial_long.parquet', columns=None)
    out['apl16'] = dict(columns=list(y.columns), holds_order=bool({'trial', 'trial_index', 'presentation',
                                                                   'onset'} & set(y.columns)),
                        question_for_yang="Yang's processed APL16 tables give each trial's odor and repeat but not its "
                                          'place in the session: in what order were the 16 odors (and the solvent) '
                                          'presented in each recording, and was it drawn anew per block and fly?')
    _write_json(out_dir() / 'order.json', out)
    return out


# ---------------------------------------------------------------------------
# the classes

def classify(scores: pd.DataFrame, lock: Optional[Mapping] = None) -> dict:
    """The entry's classes: per mixture panel and family group, the class (E134's stage rule on strength by odor) under
    each reading in its units (`PRIMARY_UNITS`) and in the other units beside; 'readable on strength' if the same
    under the four position readings, with the bounds of the family's strength error; a reading 'replaces' the lock's
    where the panel is readable under it and the panel's reproductions hold."""
    lock = load_lock() if lock is None else lock
    out = {}
    s = scores[scores.noise]
    for panel in MIXTURE_PANELS:
        rep_ok = all(v['reproduction']['passed'] for v in lock['flies'].values() if v['panel'] == panel and v['found'])
        found = all(v['found'] for v in lock['flies'].values() if v['panel'] == panel)
        for g in FAMILY_GROUPS:
            res = {}
            for units_mode in ('primary', 'other'):
                cls, vals = {}, {}
                for reading in READINGS:
                    u = PRIMARY_UNITS[reading] if units_mode == 'primary' else [x for x in UNITS
                                                                                if x != PRIMARY_UNITS[reading]][0]
                    x = s[(s.panel == panel) & (s.group == g) & (s.reading == reading) & (s.units == u)]
                    w = x[x.wiring.isin(FAMILY_WIRINGS)]
                    if len(w) < len(FAMILY_WIRINGS):
                        continue
                    lm = float(x.loo_max_view_strength.iloc[0])
                    cls[reading] = stage_class(w, lm)
                    m = x[x.wiring == 'mean'].iloc[0]
                    vals[reading] = dict(units=u, strength=float(m.view_strength), sharing=float(m.view_sharing),
                                         breadth=float(m.view_breadth), error=float(m.error), loo_max_strength=lm,
                                         floor_strength=float(m.floor_view_strength),
                                         common_strength=float(x[x.wiring == 'common'].view_strength.iloc[0]),
                                         wirings_inside=int((w.view_strength <= lm).sum()),
                                         wirings_far=int((w.view_strength > FAR_ERROR).sum()))
                pos = [cls.get(r) for r in POSITION_READINGS]
                readable = bool(None not in pos and len(set(pos)) == 1)
                st = [vals[r]['strength'] for r in POSITION_READINGS if r in vals]
                res[units_mode] = dict(classes=cls, values=vals, readable=readable,
                                       bounds=[float(min(st)), float(max(st))] if st else None,
                                       replaces={r: bool(readable and rep_ok and found) for r in POSITION_READINGS})
            out[f'{panel}/{g}'] = dict(res, reproduced=rep_ok, all_files_found=found)
    return out


def outcome(cls: Mapping, v2: Optional[Mapping] = None) -> dict:
    """The outcomes (a) to (e) the classes point to (several may hold): (e) files missing or a reproduction failing;
    (c) a mixture panel not readable; (a) readable and the kiwi panel no longer far; (b) readable and still far; (d)
    validation2's class changing with the flies read."""
    out = []
    k = cls.get('kmix/fitted', {}).get('primary', {})
    if not all(cls[p].get('reproduced') and cls[p].get('all_files_found') for p in cls):
        out.append('e')
    readable = all(cls[f'{p}/fitted']['primary']['readable'] for p in MIXTURE_PANELS)
    if not readable:
        out.append('c')
    else:
        kc = set(k.get('classes', {}).get(r) for r in POSITION_READINGS)
        out.append('b' if kc == {'far'} else 'a')
    if v2 is not None and 'fitted' in v2 and not v2['fitted'].get('readable', True):
        out.append('d')
    return dict(outcomes=out)


def summary_strength(sc: pd.DataFrame, beside: pd.DataFrame, cls: Mapping) -> pd.DataFrame:
    """One row per mixture panel and reading (its primary units, with the flies' trial noise): the family's strength
    error as fitted (the wirings' mean, least and largest, and the common error) and held, its sharing, breadth and
    population error, the flies' floor and most atypical fly on strength, the flat profile, the ORN recordings as the
    prediction, each beside model's strength error, and the class."""
    rows = []
    s = sc[sc.noise & sc.primary]
    for panel in MIXTURE_PANELS:
        for reading in READINGS:
            x = s[(s.panel == panel) & (s.reading == reading)]
            f = x[x.group == 'fitted']
            fw = f[f.wiring.isin(FAMILY_WIRINGS)].view_strength
            m = f[f.wiring == 'mean'].iloc[0]
            b = beside[(beside.panel == panel) & (beside.reading == reading) & beside.primary]
            rec = dict(panel=panel, reading=reading, units=PRIMARY_UNITS[reading],
                       family_strength=float(m.view_strength), family_strength_min=float(fw.min()),
                       family_strength_max=float(fw.max()),
                       family_strength_common=float(f[f.wiring == 'common'].view_strength.iloc[0]),
                       held_strength=float(x[(x.group == 'held') & (x.wiring == 'mean')].view_strength.iloc[0]),
                       family_sharing=float(m.view_sharing), family_breadth=float(m.view_breadth),
                       family_error=float(m.error), flies_floor_strength=float(m.floor_view_strength),
                       flies_atypical_strength=float(m.loo_max_view_strength),
                       flat_strength=float(b[b.kind == 'flat'].view_strength.iloc[0]),
                       orn_strength=float(b[b.kind == 'orn_recordings'].view_strength.iloc[0]))
            for mdl in BESIDE_MODELS:
                y = x[x.group == mdl]
                rec[f'{mdl}_strength'] = float(y.view_strength.iloc[0]) if len(y) else np.nan
            rec['class_fitted'] = cls[f'{panel}/fitted']['primary']['classes'].get(reading)
            rec['class_held'] = cls[f'{panel}/held']['primary']['classes'].get(reading)
            rows.append(rec)
    return pd.DataFrame(rows)


def part_classes() -> dict:
    """classes.json (the mixture panels' classes, validation2's where its scores exist, the outcome),
    summary_strength.csv and provenance.json (`H.code_commits()`, the lock's sha256)."""
    lock = load_lock()
    sc = pd.read_csv(out_dir() / 'rescore_scores.csv')
    cls = classify(sc, lock)
    v2p = out_dir() / 'validation2_classes.json'
    v2 = json.loads(v2p.read_text()) if v2p.exists() else None
    res = dict(mixtures=cls, validation2=v2, outcome=outcome(cls, v2), time=_now(), host=_host())
    _write_json(out_dir() / 'classes.json', res)
    summary_strength(sc, pd.read_csv(out_dir() / 'rescore_beside.csv'), cls).to_csv(
        out_dir() / 'summary_strength.csv', index=False)
    provenance('classes', lock=lock)
    return res


# the figures' colors (the dataviz skill's reference categorical slots, in fixed order) and ink
KIND_COLORS = {'component': '#2a78d6', 'mixture': '#eb6834', 'probe': '#1baf7a', 'other': '#6250d6',
               'solvent': '#8a8984'}
INK, INK_2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def part_figures() -> None:
    """position_rates.png (each fly's single-trial response rate along its session, by kind of stimulus, in the noise
    units) and strength_by_reading.png (the family's strength error at the KC somata of each mixture panel under each
    reading, its wirings, the flat profile, the ORN recordings, the flies' floor and most atypical fly)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK_2, 'axes.labelcolor': INK, 'xtick.color': INK_2,
                         'ytick.color': INK_2, 'axes.spines.top': False, 'axes.spines.right': False})
    T = pd.read_csv(out_dir() / 'position_trials.csv')
    T = T[T.units == 'noise'].assign(kind=lambda d: [stim_kind(x) for x in d.label])
    flies = list(dict.fromkeys(T.sort_values(['panel', 'fly'], ascending=[False, True]).fly))
    nc = 4
    nr = int(np.ceil(len(flies) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.3 * nr), sharey=True, squeeze=False)
    for ax, fly in zip(axs.flat, flies):
        t = T[T.fly == fly]
        for kind, col in KIND_COLORS.items():
            k = t[t.kind == kind]
            if len(k):
                ax.plot(k.trial, k.rate, 'o', ms=4, color=col, label=kind, mec='white', mew=0.6)
        ax.set_title(f"{fly} ({'kiwi' if t.panel.iloc[0] == 'kmix' else 'control'} mixtures)", fontsize=8.5,
                     color=INK)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_xlabel('trial in the session', color=INK_2)
    for ax in axs.flat[len(flies):]:
        ax.axis('off')
    for ax in axs[:, 0]:
        ax.set_ylabel('response rate,\nsingle trial')
    h, lab = axs.flat[0].get_legend_handles_labels()
    for ax in axs.flat[1:len(flies)]:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    fig.legend(h, lab, loc='lower right', ncol=1, frameon=False, bbox_to_anchor=(0.98, 0.04))
    fig.suptitle("The mixture flies' KC somata along the session (each cell's single-trial response above megamat's "
                 'soma level, baseline-noise units)', fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir() / 'position_rates.png', dpi=150)
    plt.close(fig)

    S = pd.read_csv(out_dir() / 'summary_strength.csv')
    sc = pd.read_csv(out_dir() / 'rescore_scores.csv')
    sc = sc[sc.noise & sc.primary & (sc.group == 'fitted') & sc.wiring.isin(FAMILY_WIRINGS)]
    names = {'as_uploaded': 'as uploaded\n(session z)', 'first_presentations': 'first\npresentations',
             'components_first_block': 'components,\nfirst block', 'components_second_block': 'components,\nsecond block',
             'components_both_blocks': 'components,\nboth blocks'}
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.8), sharey=True)
    for ax, panel in zip(axs, MIXTURE_PANELS):
        s = S[S.panel == panel].set_index('reading').reindex(list(READINGS))
        x = np.arange(len(READINGS))
        for i, r in enumerate(READINGS):
            w = sc[(sc.panel == panel) & (sc.reading == r)].view_strength
            ax.plot(np.full(len(w), i) - 0.12, w, 'o', ms=3.5, color='#9fc3ee', mec='white', mew=0.4,
                    label='the family, each wiring' if i == 0 else None)
        ax.plot(x - 0.12, s.family_strength, 'o', ms=8, color='#2a78d6', mec='white', mew=1.2,
                label="the family, wirings' mean")
        ax.plot(x + 0.12, s.flat_strength, 'D', ms=6, color='#8a8984', mec='white', mew=1, label='flat profile')
        ax.plot(x + 0.24, s.orn_strength, 's', ms=6, color='#1baf7a', mec='white', mew=1,
                label='ORN recordings as the prediction')
        ax.plot(x, s.flies_atypical_strength, '_', ms=22, mew=2, color='#eb6834', label='most atypical fly')
        ax.plot(x, s.flies_floor_strength, '_', ms=22, mew=2, color=INK_2, label="the flies' floor")
        ax.axhline(FAR_ERROR, color=INK_2, lw=0.8, ls=':')
        ax.text(-0.45, FAR_ERROR + 0.04, 'far', ha='left', va='bottom', color=INK_2, fontsize=8)
        ax.set_xticks(x, [names[r] for r in READINGS], fontsize=8)
        ax.set_title(f"KC somata, {'kiwi' if panel == 'kmix' else 'control'} mixtures", color=INK)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_ylim(0, max(3.0, float(np.nanmax(S.family_strength_max)) + 0.2))
    axs[0].set_ylabel('strength by odor,\none-fly units')
    axs[1].legend(loc='upper right', frameon=False, fontsize=7.5)
    fig.suptitle("The family's error on strength by odor under each reading of the flies (with the flies' trial "
                 'noise on the model)', fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir() / 'strength_by_reading.png', dpi=150)
    plt.close(fig)


def provenance(part: str, lock: Optional[Mapping] = None, **extra) -> None:
    """provenance.json: the script, `H.code_commits()`, the lock copy's sha256, host and time (one entry per part,
    the last run of each kept)."""
    from kcq import harness as H
    p = out_dir() / 'provenance.json'
    prev = json.loads(p.read_text()) if p.exists() else {}
    prev[part] = dict(script=SCRIPT, **H.code_commits(), lock_sha256=sha256_file(LOCK_COPY) if LOCK_COPY.exists()
                      else None, host=_host(), time=_now(), **extra)
    _write_json(p, prev)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'position', 'readings', 'rescore', 'rise',
                                                      'validation2', 'order', 'classes', 'figures'))
    a = ap.parse_args(argv)
    fn = dict(lock=part_lock, position=part_position, readings=part_readings, rescore=part_rescore, rise=part_rise,
              validation2=part_validation2, order=part_order, classes=part_classes, figures=part_figures)[a.part]
    fn()
    if a.part != 'classes':
        provenance(a.part)


if __name__ == '__main__':
    main()
