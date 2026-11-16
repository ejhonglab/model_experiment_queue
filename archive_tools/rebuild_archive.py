#!/usr/bin/env python3
"""Writes rebuild archives of fitted models (pyolfsysm.rebuild), so that a fit can be run on new ORN input
without this repository or the shared directory (olfdb runs them this way).

A fit of E10_joint_fit.py's family (a row of an experiment's fits.csv: the bouton base on the real or a re-drawn
wiring, breadth-weighted APL->KC weights, optionally a per-claw threshold and a glomerulus gain) is rebuilt here
as the fit built it (`E10_joint_fit.bouton_base`, `breadth_weights`, `base_kws`), and its pieces are written to
<shared>/data/model_outputs/<experiment>/rebuild/<tag>/ with `pyolfsysm.rebuild.write_archive`: the wiring and the
unscaled APL weights as fitted (parquet), the lab's ORN input, the breadth, the fit's saved APL-on and acute APL-off
outputs, and rebuild.json (the absolute parameters, the scalar model kwargs, the response criterion, and here the
wiring, gain, breadth, APL->bouton cap, operating point and the fit's provenance; the layout is described in the
archive's README.md, from archive_model_outputs.py, whose index this refreshes). The archive is then checked from
its files alone (`pyolfsysm.rebuild.check_archive`: `simulate_at` with APL on and off, and the fixed-parameter
simulator at the same two points, against the saved spike counts and responses) and moved into place only if every
count and response is identical; otherwise it stays in rebuild/.<tag>.partial/ and the script exits with status 1.

Fits with options the rebuild does not cover (another ORN input, spike-driven APL, an input-scaled threshold, a
claw transfer, the flies' pulse, pooled or divisive inhibition, a PN normalization, a fitted flat APL->KC
component, a KC refractory period or adaptation, a base other than the bouton base's) are refused.

Usage (from the repository root):
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u archive_tools/rebuild_archive.py C2
    ... archive_tools/rebuild_archive.py C2-wirings E44          # presets (PRESETS)
    ... archive_tools/rebuild_archive.py E67_apl_effect_fits:rate0.068_aplboth_k2_cap3   # any fits.csv row
    ... archive_tools/rebuild_archive.py --check C2 E44           # re-check written archives (recorded in rebuild.json)
    --out-root DIR writes under DIR/<experiment>/rebuild/ instead of the shared archive.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

ARCHIVE = H.shared_dir() / 'data' / 'model_outputs'
NAME = 'rebuild_archive'
C2 = ('E67_apl_effect_fits', 'rate0.068_aplboth_k2_cap3')
E44 = ('E44_operating_point', 'thr40_gain1.5_rate0.068_k2_cap3')
# fits.csv columns of options the rebuild does not cover: each must be missing, empty, False or the default
UNSUPPORTED = ('orn_input', 'base_params', 'kc_thr_input_exp', 'claw_thr_subgain', 'claw_thr_width', 'time_stim_end',
               'time_end', 'apl_taum', 'claw_inh_pool', 'apl_inh_divisive', 'pn_inhadd', 'pn_inhsc', 'aplkc_flat',
               'aplkc_flat_kind', 'kc_refractory_period', 'kc_adapt_inc', 'kc_adapt_tau')
PARAM_COLS = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')
ROW_COLS = ('variant', 'k', 'claw_thr', 'gain', 'gain_glomeruli', 'wiring', 'wiring_seed', 'wiring_file', 'wiring_tag',
            'apl_drive', 'loss_variant', 'apl_readout', 'rate_target', 'rate_target_spec', 'aplpn_max', 'init_from',
            'init_from_experiment', 'method', 'maxfev', 'n_evals', 'loss', 'rate_on', 'rate_off', 'rate_ratio_off_over_on',
            'thr_factor', 'apl_factor_vs_e3', 'aplpn_factor_vs_e3', 'kc_kc_inh', 'kc_kc_tau')


def presets() -> Dict[str, List[Tuple[str, str]]]:
    """{name: [(experiment folder, tag)]}: 'C2' (the breadth-2 APL-effect fit on the real wiring, E67), 'C2-wirings'
    (its fits on E72's seven re-drawn wirings, as E83 uses them), 'E44' (E44's 40 Hz operating-point refit)."""
    from experiments import E83_cross_wiring_fit as E83
    m = E83.members()
    return {'C2': [C2], 'C2-wirings': [m[w] for w in E83.ALL_WIRINGS if w != 'real'], 'E44': [E44]}


def resolve(specs: List[str]) -> List[Tuple[str, str]]:
    out, pre = [], None
    for s in specs:
        if ':' in s:
            out.append(tuple(s.split(':', 1)))
        else:
            pre = presets() if pre is None else pre
            if s not in pre:
                raise SystemExit(f'unknown preset {s!r} (presets: {sorted(pre)}; or EXPERIMENT:TAG)')
            out.extend(pre[s])
    return out


def _na(v) -> bool:
    return v is None or (isinstance(v, float) and np.isnan(v)) or (isinstance(v, str) and v.strip() == '')


def _num(v) -> Optional[float]:
    return None if _na(v) else float(v)


def fits_row(folder: str, tag: str) -> pd.Series:
    f = pd.read_csv(H.results_dir(folder) / 'fits.csv')
    r = f[f['tag'] == tag]
    if len(r) != 1:
        raise ValueError(f'{folder}: {len(r)} rows with tag {tag!r} in fits.csv')
    return r.iloc[0]


def check_supported(row: pd.Series) -> None:
    bad = {c: row[c] for c in UNSUPPORTED if c in row.index and not _na(row[c]) and row[c] is not False
           and str(row[c]) != 'False'}
    if not _na(row.get('apl_drive')) and row['apl_drive'] != 'claws':
        bad['apl_drive'] = row['apl_drive']
    if bad:
        raise NotImplementedError(f'{row["tag"]}: options the rebuild does not cover: {bad}')


def fit_provenance(folder: str, tag: str) -> dict:
    """The fit's own provenance (provenance_<tag>.json in results/, else in the archive), the fields a rebuild needs."""
    for d in (H.shared_dir() / 'results' / folder, ARCHIVE / folder):
        p = d / f'provenance_{tag}.json'
        if p.exists():
            pr = json.loads(p.read_text())
            build = (pr.get('model') or {}).get('olfsysm', {}).get('build', {})
            keep = ('generated', 'script', 'host', 'olfsysm_commit', 'al_analysis_commit', 'queue_commit',
                    'pyolfsysm_commit', 'hong2p_commit')
            return dict({k: pr[k] for k in keep if k in pr}, file=str(p.relative_to(H.shared_dir())),
                        olfsysm_source_sha256=build.get('source_sha256'), olfsysm_build_git_hash=build.get('git_hash'))
    return {}


def build_pieces(folder: str, tag: str, out_dir: Path) -> dict:
    """The fit's pieces, built as E10_joint_fit.py built them: 'params' (absolute), 'wiring' (claw x bouton, any gain
    applied), 'weights' (unscaled APL weights), 'model_kws', 'orn', 'breadth', 'outputs' (the saved fit), 'meta'."""
    from experiments import E10_joint_fit as J
    from experiments import E72_wiring_population as E72
    from experiments.E27_glomerulus_gain import ORACLE
    row = fits_row(folder, tag)
    check_supported(row)
    orn = H.apl16_orn_deltas()
    claw_thr = _num(row.get('claw_thr'))
    gain = _num(row.get('gain'))
    gloms = None
    if gain is not None:
        gloms = row['gain_glomeruli'].split(',') if not _na(row.get('gain_glomeruli')) else list(ORACLE)
    wkw: dict = {}
    wiring_meta = dict(connectome='hemibrain', claws='one row per claw (prat_claws), boutons (prat_boutons)',
                       shuffle=None, seed=None, matrix_file=None, matrix_tag=None)
    if not _na(row.get('wiring')):
        wkw = dict(wiring=str(row['wiring']), wiring_seed=int(row['wiring_seed']))
        wiring_meta.update(shuffle=wkw['wiring'], seed=wkw['wiring_seed'],
                           how='pyolfsysm.wpnkc_shuffles.shuffle_wPNKC of the hemibrain claw x bouton wiring '
                               '(NAMED_SHUFFLES[shuffle]) at seed')
    elif not _na(row.get('wiring_file')):
        p = Path(row['wiring_file'])
        p = p if p.is_absolute() else H.shared_dir() / p
        wkw = dict(wiring_matrix=pd.read_parquet(p), wiring_tag=str(row['wiring_tag']))
        rel = str(p.relative_to(H.shared_dir())) if p.is_relative_to(H.shared_dir()) else str(p)
        wiring_meta.update(matrix_file=rel, matrix_tag=wkw['wiring_tag'],
                           how='a saved claw x bouton matrix (re-drawn from the hemibrain wiring) replacing it')
    k = float(row['k'])
    initial, wnorm, W, breadth = J.bouton_base(orn, out_dir, claw_thr, gain=gain, gain_gloms=gloms, **wkw)
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], breadth, k)
    params = {c: float(row[c]) for c in PARAM_COLS}
    params['wKCKC_scale'] = float(row['kc_kc_inh'])
    tau = _num(row.get('kc_kc_tau'))
    model_kws = dict(J.base_kws(claw_thr), kc_kc_tau=float(J.TAU if tau is None else tau))
    weights = dict(wAPLKC=wAPLKC, wKCAPL=wnorm['wKCAPL'], wAPLPN=wnorm['wAPLPN'], wPNAPL=wnorm['wPNAPL'])
    outputs = E72.load_fit(folder, tag)
    if outputs is None:
        raise FileNotFoundError(f'{folder}:{tag}: no saved outputs (outputs/ or the archive)')
    cap = _num(row.get('aplpn_max'))
    meta = dict(
        experiment=folder, tag=tag,
        description=("A fitted model of the KC-remapping experiment queue (E10_joint_fit.py's family): the bouton "
                     "base with connectome APL weights, the APL->KC weights breadth-weighted, connectome-weighted "
                     "KC-KC inhibition" + ('' if claw_thr is None else f', a per-claw threshold of {claw_thr:g} Hz')
                     + ('' if gain is None else f', a gain of {gain:g} on the PN->claw weights of {len(gloms)} glomeruli')
                     + '. Absolute weights: wPNKC as stored; APL weights the stored column times '
                     "params['<name>_scale']; KC-KC weights from the connectome (fit_mb_model) times wKCKC_scale."),
        wiring=wiring_meta,
        gain=None if gain is None else dict(gain=gain, glomeruli=gloms,
                                            how='pyolfsysm.rebuild.glomerulus_gain on the wiring (already applied '
                                                'in wPNKC.parquet)'),
        breadth=dict(exponent=k, base_params=initial,
                     definition="each KC's total spike count with APL off in the bouton base at base_params (no KC-KC "
                                "inhibition), on this wiring with the gain, with the lab's 16-odor input; the APL->KC "
                                "weights are the connectome's normalized ones times ((breadth + 1) / mean(breadth + 1)) "
                                "** exponent (pyolfsysm.rebuild.breadth_weights; already applied in wAPLKC.parquet)"),
        aplpn_cap=None if cap is None else dict(aplpn_max=cap, base_wAPLPN_scale=float(initial['wAPLPN_scale']),
                                                wAPLPN_scale_max=cap * float(initial['wAPLPN_scale']),
                                                note='a bound of the fit only (wAPLPN_scale <= aplpn_max x the base\'s); '
                                                     'not needed to rebuild'),
        operating_point=dict(rate_target=_num(row.get('rate_target')), rate_on=_num(row.get('rate_on')),
                             rate_off=_num(row.get('rate_off')),
                             note='response rates (two or more spikes) of the fit with APL on and acute APL off, on the '
                                  "lab's 16-odor input"),
        orn_input="the lab's 16-odor APL-panel ORN input the fit used (kcq.harness.apl16_orn_deltas: 36 glomeruli x "
                  "(panel, odor); paa, -aPine and euc within 1 log10 unit of Yang's concentrations)",
        fit=dict(fits_csv=f'results/{folder}/fits.csv', row={c: (None if _na(row.get(c)) else
                 (row[c].item() if hasattr(row[c], 'item') else row[c])) for c in ROW_COLS if c in row.index},
                 provenance=fit_provenance(folder, tag), saved_outputs=f'data/model_outputs/{folder}/{tag}_<on|acute_off>_'
                 '<responses|spike_counts>.parquet'),
        built_by=dict(script='archive_tools/rebuild_archive.py', **H.code_commits()),
    )
    return dict(params=params, wiring=W, weights=weights, model_kws=model_kws, orn=orn, breadth=breadth,
                outputs=outputs, meta=meta)


def archive_dir(root: Path, folder: str, tag: str) -> Path:
    return root / folder / 'rebuild' / tag


def write_one(folder: str, tag: str, root: Path) -> dict:
    from pyolfsysm import rebuild as R
    t0 = time.time()
    final = archive_dir(root, folder, tag)
    partial = final.with_name(f'.{tag}.partial')
    if partial.exists():
        shutil.rmtree(partial)
    P = build_pieces(folder, tag, H.outputs_dir(NAME))
    R.write_archive(partial, params=P['params'], wiring=P['wiring'], weights=P['weights'], model_kws=P['model_kws'],
                    orn_deltas=P['orn'], outputs=P['outputs'], breadth=P['breadth'], meta=P['meta'])
    del P
    res = R.check_archive(partial, record=True)
    print(json.dumps({f'{folder}:{tag}': res}, indent=1), flush=True)
    if not res['ok']:
        print(f'!!! {folder}:{tag}: the rebuild differs from the saved fit; archive left in {partial}', file=sys.stderr)
        return dict(res, path=str(partial), seconds=time.time() - t0)
    if final.exists():
        shutil.rmtree(final)
    partial.rename(final)
    print(f'{folder}:{tag}: archive written to {final} ({time.time() - t0:.0f}s)', flush=True)
    return dict(res, path=str(final), seconds=time.time() - t0)


def check_one(folder: str, tag: str, root: Path) -> dict:
    from pyolfsysm import rebuild as R
    d = archive_dir(root, folder, tag)
    res = R.check_archive(d, record=True)
    print(json.dumps({f'{folder}:{tag}': res}, indent=1), flush=True)
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('variants', nargs='+', help='presets (C2, C2-wirings, E44) or EXPERIMENT:TAG')
    ap.add_argument('--check', action='store_true', help='only re-check written archives')
    ap.add_argument('--out-root', type=Path, default=ARCHIVE, help='default: <shared>/data/model_outputs')
    args = ap.parse_args()
    ok = True
    for folder, tag in resolve(args.variants):
        res = (check_one if args.check else write_one)(folder, tag, args.out_root)
        ok = ok and res['ok']
    if args.out_root == ARCHIVE:
        # the index's description and list of the rebuild archives
        from archive_tools.archive_model_outputs import update_index
        update_index([])
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
