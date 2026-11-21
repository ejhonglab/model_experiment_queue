#!/usr/bin/env python3
"""Copies an experiment's model outputs, with their parameters, to the shared archive.

The archive is <shared>/data/model_outputs/E<n>_<short-description>/ (same name as the
experiment's results/ and outputs/ folders), organized by experiment regardless of the host
that ran it. Per experiment it holds: every KC x odor matrix the run saved (parquet;
<variant>_<condition>_spike_counts.parquet = total spikes per KC per odor,
<variant>_<condition>_responses.parquet = at least 2 spikes), any other data files the run
wrote next to them (tuned_params.json, wPNKC matrices, sweep tables, run logs),
params.csv (the stats table's identifying and parameter columns: tag, variant, base, seed,
thresholds, APL / KC-KC scales, sweep parameters), provenance.json from results/, and a
README.md listing the variants. The top-level README.md gets one row per experiment.

    ./venv/bin/python archive_tools/archive_model_outputs.py E16_claw_saturation [E19_... ...]
    ./venv/bin/python archive_tools/archive_model_outputs.py --all          # every outputs/E* folder
    ./venv/bin/python archive_tools/archive_model_outputs.py --inputs       # ORN input + base wirings
    ./venv/bin/python archive_tools/archive_model_outputs.py --all --params-only   # rebuild every params.csv
    ./venv/bin/python archive_tools/archive_model_outputs.py E98_apl_response_strength --exclude tmp,cache,B/old_46stim
        # leave out folders under outputs/<name>/ that hold no model outputs of the experiment (pytest temp
        # folders, caches of other experiments' runs, superseded runs); glob patterns allowed (B/*/*_held_*)

Experiments whose step folders repeat their files' names (FOLDER_PREFIXES; E102's B/<step>/<wiring>_...) have the
step put before each name in the archive, so that every variant's name is unique (olfdb finds a model's files by its
tag in their names, not by their folder); a copy archived earlier under the bare name is removed. DEFAULT_EXCLUDE
lists, per experiment, what is never copied (scratch checks, copies of runs archived elsewhere).

File and folder names are ASCII only (CLAUDE.md); names with other characters are rejected.
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import shutil
import socket
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

from kcq import harness as H

ARCHIVE = H.shared_dir() / 'data' / 'model_outputs'
# identifying columns, every model parameter a script has put in its table (fit_mb_model kwargs,
# olfsysm parameters, input / wiring modifications), fit settings, and a few output rates. Add
# to this list whenever a script records a new parameter column; params.csv is what the
# archive README points at for regenerating a variant.
PARAM_COLS = ['tag', 'variant', 'base', 'config', 'kind', 'point', 'part', 'seed', 'target_rate', 'target',
    'fixed_thr', 'fixed_thr_mean', 'fixed_thr_is_vector', 'thr_factor', 'thr_scale_homeo_full',
    'thr_scale_homeo_data', 'stage0_fixed_thr', 'stage0_thr_factor', 'wAPLKC', 'wKCAPL',
    'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'kc_apl_factor', 'wKCKC_scale',
    'kc_kc_inh', 'kc_kc_tau', 'wkckc', 'claw_thr', 'claw_thr_width', 'claw_thr_subgain', 'claw_sat', 'form',
    'claw_pow', 'claw_pow_ref', 'bouton_ves_p', 'bouton_tau_r', 'alpha', 'gain_min', 'gain_max', 'gain', 'glomeruli',
    'gain_glomeruli', 'factor', 'scaled_odors', 'k', 'm', 'breadth', 'breadth_k', 'n_claws_active_to_spike',
    'coinput_factor', 'apl_unit_coupling', 'coupling', 'kc_set', 'kc_scale', 'bouton_scale', 'ext_inh_scale',
    'pn_inhadd', 'pn_inhsc', 'pn_apl_multiplier', 'multibouton_claws', 'claw_weight', 'wiring', 'wiring_seed',
    'wiring_file', 'wiring_tag',
    'apl_drive', 'pn_claw_to_apl', 'apl_inh_divisive', 'bouton_inh_divisive', 'base_params', 'breadth_def', 'apl_kc_factor',
    'arm', 'orn_input', 'multiplier_min', 'multiplier_max', 'multiplier_geomean', 'scale_min', 'scale_max', 'cal_readout_param',
    'init_params', 'aplpn_max', 'aplpn_fixed', 'aplkc_min',
    'kc_refractory_period', 'kc_adapt_inc', 'kc_adapt_tau', 'adapt_strength', 'apl_factor', 'ref_experiment', 'ref_fit',
    'time_stim_end', 'time_end', 'apl_taum', 'count_kernel_rise', 'count_kernel_decay', 'option', 'value', 'fixed_thr_off',
    'kc_thr_input_exp', 'breadth_kept', 'prepass_thr_factor', 'prepass_fixed_thr',
    'rel', 'start_experiment', 'start_tag', 'breadth_from', 'kckc', 'kckc_factor',
    'max_count_on', 'max_count_off', 'claw_inh_pool', 'allow_net_inh_per_claw', 'ext_inh_kc_max',
    'aplpn_hold', 'aplkc_flat', 'aplkc_flat_size', 'aplkc_flat_unit', 'stimuli',
    'apl_residual', 'apl_residual_fitted', 'apl_residual_lo', 'apl_residual_hi', 'off_wAPLKC_scale', 'off_wAPLPN_scale',
    'change', 'search', 'size', 'panel',
    'init_from', 'init_thr', 'loss_variant', 'method', 'maxfev', 'kckc0', 'thr0_factor',
    'n_kcs', 'n_glomeruli', 'rate_on', 'rate_off', 'rate_acute_off', 'rate_homeo_full_off',
    'rate_homeo_data_off', 'error']
STATS_NAMES = ['stats_variants.csv', 'models.csv', 'wkckc_controls.csv', 'sweep.csv']
# archive folders whose params.csv (and file names / README) were rebuilt by hand, so that olfdb lists one
# row per saved model; re-archiving them would overwrite that with the stats tables' rows (and, for E6, copy
# back the old file names), so they are skipped unless --force-hand-built is passed. Each folder's
# params_migration_log.csv says what was done.
HAND_BUILT = {
    'E6_readout_window': 'files renamed <base>_<window>_<condition>_spike_counts.parquet, one row per base x window (2026-11-09)',
    'E36_spike_apl': 'rows spk_base_thr40 / spk_base_thr70 added from tuned_params.json (2026-11-09)',
    'E43_model_trial_noise': 'params.csv is the two models ref40 / base70; the noise statistics are in noise_statistics_rows.csv (2026-11-09)',
    'E49_low_conc_kc': 'params.csv written by hand: E49 has none of the tables params.csv is built from (2026-11-09)',
}

# experiments whose params.csv is built from these tables only (names in results/<name>/ or outputs/<name>/): their other
# summary tables are readings, not model variants (E120: part A's a1_ / a2_ / a3_summary.csv read other experiments'
# archives, b_stage_summary_*.csv reads part B's runs; variants.csv has one row per saved stage run)
PARAMS_TABLES = {
    'E120_bouton_stage_measured': ('variants.csv',),
}

# experiments whose step folders repeat their files' names: per folder pattern (fnmatch on the file's folder under
# outputs/<name>/), the prefix put before the file's name in the archive, '{0}', '{1}', ... the folder's path parts,
# so that each variant's name is unique (E102: B/<step>/<wiring>_held_on_... in every step's folder)
FOLDER_PREFIXES = {
    'E102_antennal_lobe_normalization': (('B/*', 'B_{1}_'), ('D/runs/B/*', 'D_{3}_'), ('E/B/*', 'E_{2}_'),
                                         ('H/B/*', 'H_{2}_')),
}
# per experiment, the paths under outputs/<name>/ never copied (as --exclude, which adds to them): E102's start checks
# (one evaluation of E10 each, in a scratch folder whose file names are the refits'), the refits' 16-odor copies under
# D/runs (E10's own outputs are archived at the top) and the references part E's score copies beside its points
DEFAULT_EXCLUDE = {
    'E102_antennal_lobe_normalization': ('D/start_check', 'D/runs/B/*/*_held_*', 'E/B/*/unchanged_none_*',
                                         'E/B/*/refit_D_*', 'E/B/*/E95_f1_*'),
    # E97: part B's per-point claims (lane bookkeeping) and the scoring caches
    'E97_glomerulus_claw_shares': ('claims', 'distances_cache.pkl', 'coresponse_cache.pkl'),
}
# folders of outputs withheld from queue work (e.g. E97's withheld_2026-11-15/, built from connectome tables that queue
# work does not use; CLAUDE.md, "Connectome data"): never copied, at any depth
WITHHELD_PREFIX = 'withheld_'


def excluded(rel: str, exclude) -> bool:
    """Whether the path `rel` (relative to outputs/<name>/) is excluded: under a folder of withheld outputs
    (WITHHELD_PREFIX), equal to or under one of `exclude`, or matching one as a glob pattern."""
    if any(part.startswith(WITHHELD_PREFIX) for part in rel.split('/')[:-1]):
        return True
    return any(rel == x or rel.startswith(x.rstrip('/') + '/') or fnmatch.fnmatchcase(rel, x) for x in exclude)


def archive_rel(name: str, rel: Path) -> Path:
    """Where a file of outputs/<name>/ goes in the archive: `rel`, with FOLDER_PREFIXES' prefix for its folder put
    before its name when one applies."""
    folder = rel.parent.as_posix()
    for pattern, prefix in FOLDER_PREFIXES.get(name, ()):
        if fnmatch.fnmatchcase(folder, pattern):
            return rel.parent / (prefix.format(*rel.parent.parts) + rel.name)
    return rel


# condition names used in the output file names (longest first, so 'acute_off' wins over 'off')
CONDITIONS = ('homeo_full_off', 'homeo_data_off', 'acute_off', 'off', 'on')


def split_variant_condition(filename: str):
    """'<variant>_<condition>_spike_counts.parquet' -> (variant, condition); condition may
    contain underscores (acute_off), so it is matched against CONDITIONS."""
    stem = filename
    for suffix in ('_spike_counts.parquet', '_responses.parquet'):
        if stem.endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    for c in CONDITIONS:
        if stem.endswith('_' + c):
            return stem[:-len(c) - 1], c
    return stem, ''


def ascii_only(name: str) -> None:
    if not name.isascii():
        raise ValueError(f'non-ASCII name: {name!r}')


def archive_experiment(name: str, host: str, params_only: bool = False, force_hand_built: bool = False,
                       exclude: tuple = ()) -> dict:
    """Copies outputs/<name>/ to the archive and (re)builds its params.csv and README.

    With `params_only`, only params.csv and the provenance files are refreshed (from the shared
    results/ tables and the tables already in the archive folder), so it works for experiments
    whose local outputs are on another host; the README row is left as it is. `exclude`: paths relative to
    outputs/<name>/ (folders or files) that are not copied.
    """
    src = H.HERE / 'outputs' / name
    res = H.results_dir(name) if (H.shared_dir() / 'results' / name).exists() else None
    ascii_only(name)
    if name in HAND_BUILT and not force_hand_built:
        print(f'!!! {name} SKIPPED: its params.csv was rebuilt by hand ({HAND_BUILT[name]}); re-archiving would overwrite it. '
            'Update it by hand (and its params_migration_log.csv), or pass --force-hand-built and redo the hand edits.', file=sys.stderr)
        return dict(name=name, host=host, n_files=0, n_variants=0, conditions='', skipped=True)
    dst = ARCHIVE / name
    if params_only:
        if not dst.exists():
            raise FileNotFoundError(f'{dst} (nothing archived yet; run without --params-only on the host that ran it)')
        src = dst
    elif not src.exists():
        raise FileNotFoundError(src)
    dst.mkdir(parents=True, exist_ok=True)
    n_files = 0
    exclude = tuple(DEFAULT_EXCLUDE.get(name, ())) + tuple(exclude)
    for p in ([] if params_only else sorted(src.rglob('*'))):
        if not p.is_file() or p.suffix in ('.pyc', '.html'):
            continue   # .html: artifact pages live in reports/, not here
        rel = p.relative_to(src)
        if excluded(rel.as_posix(), exclude):
            continue
        ascii_only(rel.as_posix())
        arel = archive_rel(name, rel)
        target = dst / arel
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists() or target.stat().st_size != p.stat().st_size:
            shutil.copy2(p, target)
        if arel != rel and (dst / rel).is_file():
            # the copy archived before the folder prefix applied
            (dst / rel).unlink()
        n_files += 1
    # parameters: the stats tables' identifying / parameter columns
    frames = []
    # the fixed names plus any per-experiment summary / stats / fits table (atlas's
    # experiments write E<n>_summary.csv, window_stats.csv, fits.csv, ...)
    tables = []
    for d in ([res] if res else []) + [src]:
        if name in PARAMS_TABLES:
            tables.extend(d / t for t in PARAMS_TABLES[name])
            continue
        for sname in STATS_NAMES:
            tables.append(d / sname)
        for pat in ('*_summary.csv', 'stats*.csv', 'window_stats.csv', 'fits.csv', 'variants.csv'):
            tables.extend(sorted(d.glob(pat)))
    seen = set()
    for sp in tables:
        if sp in seen:
            continue
        seen.add(sp)
        for _ in (0,):
            if sp.exists():
                df = pd.read_csv(sp)
                keep = [c for c in PARAM_COLS if c in df.columns]
                if keep:
                    d = df[keep].copy(); d.insert(0, 'source_table', sp.name); frames.append(d)
    params = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if len(params):
        params.to_csv(dst / 'params.csv', index=False)
    missing_params = not len(params) and any(dst.rglob('*_spike_counts.parquet'))
    if missing_params:
        # E49 (2026-11-08) was archived with spike counts and no params.csv, so olfdb listed none of its models
        print(f'!!! {name}: spike count files archived but NO params.csv built: none of the tables it is built from '
            f'({", ".join(STATS_NAMES)}, *_summary.csv, stats*.csv, window_stats.csv, fits.csv, variants.csv) exists in '
            f'results/{name}/ or outputs/{name}/. Write params.csv by hand (one row per saved model, tag = the file names\' '
            'variant) or save the run\'s table under one of those names, then re-run.', file=sys.stderr)
    if res:
        # provenance.json, plus the per-fit provenance_<tag>.json the joint-fit scripts write
        for pf in sorted(res.glob('provenance*.json')):
            shutil.copy2(pf, dst / pf.name)
    if params_only:
        return dict(name=name, host=host, n_files=len(list(dst.rglob('*'))), n_variants=params['tag'].nunique() if 'tag' in params.columns else len(params), conditions='', params_only=True,
            missing_params=missing_params)
    parsed = [split_variant_condition(p.name) for p in list(dst.rglob('*_spike_counts.parquet'))
        + list(dst.rglob('*_responses.parquet'))]
    variants = sorted({v for v, _ in parsed})
    conds = sorted({c for _, c in parsed if c})
    lines = [f'# {name}: model outputs', '',
        f'Copied {datetime.now():%Y-%m-%d %H:%M} from host {host} (`outputs/{name}/` of the model_experiment_queue checkout). '
        f'{n_files} files; {len(variants)} variants; conditions: {", ".join(conds) if conds else "(see file names)"}.', '',
        'Files: `<variant>_<condition>_spike_counts.parquet` (KCs x odors, total spikes per KC per odor), '
        '`<variant>_<condition>_responses.parquet` (same shape, 1 = at least 2 spikes); rows indexed by hemibrain '
        'KC body ID (and type), columns the odor at its concentration. `params.csv`: identifying and parameter '
        'columns of the stats tables (tag = variant name used in the file names); `provenance.json`: script, '
        'commit hashes, base kwargs and sweep settings; run logs as written. Regenerate with the script named in '
        'provenance.json at the recorded commit (`queue_commit`, of ejhonglab/model_experiment_queue; in provenance '
        'written before 2026-11-11, `al_analysis_commit`, of al_analysis, where the scripts were), with the base kwargs '
        'there and the tuned parameters in params.csv.', '']
    if name in FOLDER_PREFIXES:
        lines += ['In ' + ', '.join(f'`{pat}`' for pat, _ in FOLDER_PREFIXES[name]) + ' each file\'s name begins with its folder\'s '
                  'step (' + ', '.join(f'`{pre}`' for _, pre in FOLDER_PREFIXES[name]) + ', the folder\'s path parts from 0), '
                  f'which `outputs/{name}/` does not have: the same names repeat in every step\'s folder there.', '']
    if variants:
        lines += ['Variants:', ''] + [f'- `{v}`' for v in variants] + ['']
    (dst / 'README.md').write_text('\n'.join(lines))
    return dict(name=name, host=host, n_files=n_files, n_variants=len(variants), conditions=', '.join(conds), missing_params=missing_params)


def save_inputs() -> None:
    from experiments.E3_drive_dependent_gain import BOUTON, BASES
    from experiments.E4_substrate import CLAWS
    from pyolfsysm.mb_model import fit_mb_model
    dst = ARCHIVE / 'inputs'; dst.mkdir(parents=True, exist_ok=True)
    orn = H.apl16_orn_deltas()
    orn.to_csv(dst / 'orn_deltas_apl16.csv')
    bases = {'claw_base': dict(CLAWS), 'claw_dynamics_base': dict(BASES['claw-dynamics']), 'bouton_base': dict(BOUTON),
        'hemibrain_wd20_base': dict(BASES['hemibrain-wd20'])}
    meta = {}
    for bname, kws in bases.items():
        _, _, wPNKC, pd0 = fit_mb_model(orn, fixed_thr=100.0, wAPLKC=0.1, **kws, make_plots=False, silent=True)
        wPNKC.to_parquet(dst / f'wPNKC_{bname}.parquet')
        meta[bname] = dict(kwargs={k: (v if isinstance(v, (int, float, str, bool)) else str(v)) for k, v in kws.items()},
            shape=list(wPNKC.shape), index_levels=list(wPNKC.index.names), column_levels=list(wPNKC.columns.names))
    (dst / 'bases.json').write_text(json.dumps(meta, indent=1))
    (dst / 'README.md').write_text('# Inputs shared by all experiments\n\n'
        '- `orn_deltas_apl16.csv`: the ORN input matrix (glomeruli x the 16 APL16 odors; odor-evoked rate changes in Hz '
        'relative to the Hallem spontaneous rates olfsysm adds), as `kcq.harness.apl16_orn_deltas()` returns it. '
        'paa, euc and -aPine are at substituted ORN concentrations (see the harness).\n'
        '- `wPNKC_<base>.parquet`: the PN->KC claw weight matrix of each standard base (rows: claws with KC id / type / '
        'coordinates; columns: glomeruli, or glomerulus / PN / bouton for the bouton base); `bases.json`: the fit_mb_model '
        'kwargs that build each. Experiments that modified the wiring saved their own matrices (substrate experiment E4; '
        'per-glomerulus gains of E20 are (median peak / peak)^alpha on the bouton base columns, peaks in its glomerulus_peaks.csv).\n')
    print('inputs saved:', sorted(p.name for p in dst.iterdir()))


def update_index(rows: list) -> None:
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    idx = ARCHIVE / 'README.md'
    head = ('# Model outputs of the KC-remapping experiments\n\n'
        'One folder per experiment, `E<n>_<short-description>` (the same name as its `results/` folder and its entry in '
        '`experiment_priorities.md`), regardless of the host that ran it. Each holds the KC x odor spike-count and response '
        'matrices of every model variant and condition, the variants\' parameters (`params.csv`, `tuned_params.json` where '
        'written), `provenance.json` (script, commit hashes, base kwargs, sweep settings), the run logs and a README listing '
        'the variants. `inputs/` holds what all experiments share: the ORN input matrix and the standard PN->KC wiring matrices. '
        'File and folder names are ASCII only. Add or refresh a folder with `archive_tools/archive_model_outputs.py` from the repository root.\n\n'
        + REBUILD_README + '\n\n'
        '| experiment | host | files | variants | conditions | copied |\n|---|---|---|---|---|---|\n')
    existing = {}
    if idx.exists():
        for line in idx.read_text().splitlines():
            if line.startswith('| E') or line.startswith('| inputs'):
                existing[line.split('|')[1].strip()] = line
    for r in rows:
        existing[r['name']] = f"| {r['name']} | {r['host']} | {r['n_files']} | {r['n_variants']} | {r['conditions']} | {datetime.now():%Y-%m-%d %H:%M} |"
    idx.write_text(head + '\n'.join(existing[k] for k in sorted(existing, key=lambda s: (len(s.split('_')[0]), s))) + '\n'
        + rebuild_table())


# the rebuild archives (archive_tools/rebuild_archive.py, pyolfsysm.rebuild), described in the index's head
REBUILD_README = (
    '`<experiment>/rebuild/<tag>/` (some fits; listed at the end): a rebuild archive of one fitted model, enough to run '
    'it on any ORN input without the model_experiment_queue repository or the rest of this directory, with olfsysm\'s '
    '`pyolfsysm.rebuild` (`load_archive`; `simulate_at` or `simulate_archive` to run it as fitted, returning '
    '`fit_mb_model`\'s outputs; `make_simulator` with `tune_threshold` or `pyolfsysm.tuning.hold_rates` to re-tune it). '
    'Files: `rebuild.json` (format `pyolfsysm.rebuild/1`; `params`: the fitted threshold and APL / KC-KC scales as '
    'absolute values; `model_kws`: the scalar `fit_mb_model` kwargs; `apl_off`: the acute APL-off condition; '
    '`response_criterion`: two or more spikes; `files`: each file with its sha256; `wiring` (hemibrain, or the shuffle '
    'and seed, or the re-drawn matrix it came from), `gain` (gain and glomeruli), `breadth` (exponent, base parameters, '
    'definition), `aplpn_cap`, `operating_point`, `fit` (its fits.csv row and provenance, with the olfsysm commit it ran '
    'at), `built_by`, `checks`), `wPNKC.parquet` (the claw x bouton wiring as fitted, any gain applied), '
    '`wAPLKC.parquet`, `wKCAPL.parquet`, `wAPLPN.parquet`, `wPNAPL.parquet` (the unscaled APL weights as fitted, APL->KC '
    'breadth-weighted; absolute weight = the column times `params[\'<name>_scale\']`), `orn_deltas.parquet` (the fit\'s '
    'ORN input), `breadth.parquet`, and `expected_<responses|spike_counts>_<on|acute_off>.parquet` (the fit\'s saved '
    'outputs). Written by `archive_tools/rebuild_archive.py`, which keeps an archive only if its own rebuild reproduces '
    'every saved spike count and response with APL on and off (`--check` re-checks; each check is appended to '
    '`checks`).')


def rebuild_table() -> str:
    """The index's list of rebuild archives (`<experiment>/rebuild/<tag>/rebuild.json`), with their last check."""
    rows = []
    for js in sorted(ARCHIVE.glob('E*/rebuild/*/rebuild.json')):
        try:
            rec = json.loads(js.read_text())
        except (OSError, ValueError):
            continue
        ch = (rec.get('checks') or [{}])[-1]
        diffs = [ch.get(r, {}).get(c, {}).get('spike_count_diffs') for r in ('fit_mb_model', 'simulator') for c in ('on', 'off')]
        rows.append(f"| `{js.parent.relative_to(ARCHIVE).as_posix()}` | {'yes' if ch.get('ok') else 'NO'} | "
            f"{' / '.join('-' if d is None else str(d) for d in diffs)} | {str(ch.get('time', ''))[:16].replace('T', ' ')} |")
    if not rows:
        return ''
    return ('\n## Rebuild archives\n\n| archive | reproduces the fit | differing spike counts (fit_mb_model on / off, '
        'simulator on / off) | last checked (UTC) |\n|---|---|---|---|\n' + '\n'.join(rows) + '\n')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('experiments', nargs='*')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--inputs', action='store_true')
    ap.add_argument('--params-only', action='store_true',
        help='only rebuild params.csv from the shared results tables (works for any host\'s experiments)')
    ap.add_argument('--exclude', default='',
        help='comma-separated paths under outputs/<name>/ not to copy (e.g. tmp,cache: pytest temp folders, caches)')
    ap.add_argument('--force-hand-built', action='store_true',
        help=f'also re-archive the folders whose params.csv was rebuilt by hand ({", ".join(HAND_BUILT)}), overwriting it')
    args = ap.parse_args()
    host = socket.gethostname()
    if args.inputs:
        save_inputs()
    names = args.experiments
    if args.all:
        where = ARCHIVE if args.params_only else (H.HERE / 'outputs')
        names = sorted(p.name for p in where.glob('E*') if p.is_dir())
    rows = []
    for n in names:
        r = archive_experiment(n, host, params_only=args.params_only, force_hand_built=args.force_hand_built,
                               exclude=tuple(x.strip() for x in args.exclude.split(',') if x.strip()))
        if r.get('skipped'):
            continue
        rows.append(r)
        print(f"{n}: {r['n_files']} files, {r['n_variants']} variants -> {ARCHIVE / n}" + (' (params.csv only)' if args.params_only else ''))
    if rows and not args.params_only:
        update_index(rows)
    missing = [r['name'] for r in rows if r.get('missing_params')]
    if missing:
        sys.exit(f'!!! no params.csv for {", ".join(missing)} (see above); olfdb will list none of their models')


if __name__ == '__main__':
    main()
